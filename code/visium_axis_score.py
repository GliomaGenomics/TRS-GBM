"""JBS axis scores of Visium spots in 11 glioblastomas (de Jong 2025).

Spots and section pseudobulks scored as log2(FPKM + 1) · axis loadings, for the
standard-Visium tumours the config lists; sections with fewer than
MIN_MEDIAN_COUNTS median counts per spot are left out. Panels:

    visium_axis_score_maps               spot scores over each section's H&E, coloured
                                         from its 5th to 95th percentile
    visium_section_axis_score            section scores by tumour and site, with a
                                         permutation p per tumour
    visium_axis_score_dispersion         within- against between-tumour site distances
    visium_axis_score_by_region          mean score per IvyGAP region and section
    visium_axis_score_region_difference  each region minus cellular tumour, Wilcoxon
                                         signed-rank, Benjamini-Hochberg

Example:
    pixi run visium_axis_score
"""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq
import yaml
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy import sparse
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

from common import (cohort, config, pair_distances, patient_clustering, resolve, save_panel,
                    set_paper_theme, within_patient_p)

DATASET = "deJong_2025"
CMAP = "jet_r"
MAP_SIZE_INCH = 2
MAP_DPI = 300
# The grid of maps is 23 inches wide, about three times its printed width:
# 150 dpi still gives over 450 dpi across a full page (183 mm).
GRID_DPI = 150
# Margin around a section's spots in its map, as a fraction of their extent.
MAP_PADDING = 0.10
# Each map's colour scale runs between these percentiles of its section's spot
# scores, so that a few extreme spots do not compress the rest into one colour.
COLOUR_PERCENTILES = (5, 95)
# Fewer counts than this in a section's median spot marks a failed section: one
# section, about 1% of its tumour's typical depth; the next is 344.
MIN_MEDIAN_COUNTS = 200
# Each tumour's block in an 8 x 11 grid of maps: top row, left column, rows,
# columns, and whether its sections, in order, fill the block along rows (a
# site's two sections side by side) or down columns.
GRID_ROWS, GRID_COLUMNS = 8, 11
TILING = {
    "AT3": (0, 0, 4, 2, "rows"), "AT4": (0, 2, 4, 2, "rows"), "AT5": (0, 4, 4, 2, "rows"),
    "AT6": (0, 6, 4, 2, "rows"), "AT9": (0, 8, 4, 2, "rows"), "AT14": (0, 10, 4, 1, "rows"),
    "AT7": (4, 0, 2, 5, "columns"), "AT11": (6, 0, 2, 5, "columns"),
    "AT10": (4, 5, 4, 2, "rows"), "AT12": (4, 7, 4, 2, "rows"), "AT13": (4, 9, 4, 2, "rows"),
}
# Width of the strip right of the grid that holds the colour bar, as a fraction of a map.
COLOUR_BAR_STRIP = 0.7
# White space around each tumour's block, as a fraction of a map.
BLOCK_MARGIN = 0.15
# Every map is drawn the same size: 90% of a cell of a two-column block, which
# leaves white space between neighbouring maps.
MAP_SIDE = 0.9 * (2 - 2 * BLOCK_MARGIN) / 2
# Gap between a tumour's maps and the box drawn around them, as a fraction of a map.
BOX_PAD = 0.03
# Half the width of a tumour's slot that its sites, or its regions, are spread across.
SITE_SPREAD = 0.3
# Half the width of the spread of a tumour's sections within each region's place.
SECTION_SPREAD = 0.04
MIN_REGION_OVERLAP = 0.5
MIN_REGION_SPOTS = 10
# Fewest sites a tumour needs for its own clustering test.
MIN_SITES = 3
N_PERMUTATIONS = 10_000
SEED = 11
DISPERSION_COLOUR = "silver"


def palette(group: str) -> dict:
    with open(resolve("config/palette.yaml")) as handle:
        return yaml.safe_load(handle)[group]


def section_paths(entry: dict) -> list[Path]:
    """Sections of the tumours the config lists."""
    return sorted(path for path in resolve(entry["directory"]).glob("*.h5ad")
                  if tumour(path.stem) in entry["tumours"])


def tumour(section: str) -> str:
    return section.split("-")[0]


def tumour_order(names) -> list[str]:
    """Tumours in numerical order: AT3 before AT10."""
    return sorted(set(names), key=lambda name: int(name.removeprefix("AT")))


def site(section: str) -> str:
    """Sampling site from a section name, e.g. "4" in AT10-BRA-5-FO-4_1."""
    return section.split("FO")[1].split("-")[1].split("_")[0]


def axis_genes(entry: dict, symbols: dict) -> pd.DataFrame:
    """Axis genes by symbol, with their loading and the length of their own Ensembl ID.

    Genes with no symbol or no recorded length cannot be matched or normalised
    and are left out.
    """
    loadings = pd.read_csv(resolve(entry["axis_loadings"]), index_col=0).iloc[:, 0]
    annotation = pd.read_csv(resolve(symbols["path"]), index_col=0)
    return pd.DataFrame({
        "symbol": annotation[symbols["column"]].reindex(loadings.index),
        "loading": loadings,
        "length_kb": annotation[symbols["length_column"]].reindex(loadings.index) / 1e3,
    }).dropna()


def measured_symbols(symbols: dict) -> pd.Index:
    """Symbols with a recorded length: the genes a spot's library size is taken over."""
    annotation = pd.read_csv(resolve(symbols["path"]), index_col=0)
    return pd.Index(annotation.dropna(subset=[symbols["length_column"]])[symbols["column"]].dropna().unique())


def read_section(path: Path) -> ad.AnnData:
    data = ad.read_h5ad(path)
    # Ten gene symbols appear twice in the counts; the first of each is kept by name.
    data.var_names_make_unique()
    return data


def features(data: ad.AnnData, entry: dict, feature_type: str) -> ad.AnnData:
    """One feature set of a section (expression counts, histology overlaps), image and coordinates kept."""
    return data[:, data.var[entry["feature_type_column"]] == feature_type].copy()


def spot_scores(data: ad.AnnData, genes: pd.DataFrame, measured: pd.Index, pseudocount: float) -> pd.Series:
    """Place every spot of a section on the JBS axis.

    The score is summed over the stored counts alone and offset by the score of
    a spot of nothing but zeros, which is the same for every spot.
    """
    counts = sparse.csr_matrix(data.X, dtype=np.float64)
    names = data.var_names
    library_sizes = np.asarray(counts[:, names.isin(measured)].sum(axis=1)).ravel() / 1e6

    present = genes[genes["symbol"].isin(names)]
    columns = pd.Series(np.arange(len(names)), index=names)[present["symbol"]].to_numpy()
    logged = counts[:, columns].tocoo()
    fpkm = logged.data / (library_sizes[logged.row] * present["length_kb"].to_numpy()[logged.col])
    logged.data = np.log2(fpkm + pseudocount) - np.log2(pseudocount)

    weights = present["loading"].to_numpy()
    zeros = np.log2(pseudocount) * weights.sum()
    return pd.Series(logged.tocsr() @ weights + zeros, index=data.obs_names)


def section_score(data: ad.AnnData, genes: pd.DataFrame, measured: pd.Index, pseudocount: float) -> float:
    """Place a section's pseudobulk, the sum of its spots' counts, on the JBS axis."""
    pseudobulk = ad.AnnData(sparse.csr_matrix(data.X.sum(axis=0)), var=data.var)
    return spot_scores(pseudobulk, genes, measured, pseudocount).iloc[0]


def spot_regions(overlaps: ad.AnnData, regions: dict) -> pd.Series:
    """Each spot's histological region: the one it overlaps by at least half.

    A spot is left unassigned when it overlaps no region that much, or more than
    one, and so is every spot of an unannotated section.
    """
    overlap = pd.DataFrame(sparse.csr_matrix(overlaps.X).toarray(), index=overlaps.obs_names,
                           columns=overlaps.var_names)
    by_region = pd.DataFrame({region: overlap.reindex(columns=names, fill_value=0).sum(axis=1)
                              for region, names in regions.items()})
    assigned = (by_region >= MIN_REGION_OVERLAP).sum(axis=1) == 1
    return by_region.idxmax(axis=1).where(assigned)


def section_map(data: ad.AnnData) -> np.ndarray:
    """A section's spots over its H&E image, cropped square around the spots, as RGBA pixels."""
    library = next(iter(data.uns["spatial"]))
    scale = float(data.uns["spatial"][library]["scalefactors"]["tissue_hires_scalef"])

    xy = np.asarray(data.obsm["spatial"])[:, :2] * scale
    (x_min, y_min), (x_max, y_max) = xy.min(axis=0), xy.max(axis=0)
    half = max(x_max - x_min, y_max - y_min) * (0.5 + MAP_PADDING)
    centre_x, centre_y = (x_min + x_max) / 2, (y_min + y_max) / 2
    height, width = data.uns["spatial"][library]["images"]["hires"].shape[:2]
    left, right = max(0.0, centre_x - half), min(float(width - 1), centre_x + half)
    top, bottom = max(0.0, centre_y - half), min(float(height - 1), centre_y + half)

    fig, ax = plt.subplots(figsize=(MAP_SIZE_INCH, MAP_SIZE_INCH), dpi=MAP_DPI)
    fig.patch.set_alpha(0)
    sq.pl.spatial_scatter(
        data, color=["JBS axis score"], cmap=CMAP, colorbar=False, shape="hex",
        norm=Normalize(*np.percentile(data.obs["JBS axis score"], COLOUR_PERCENTILES)),
        img=True, img_res_key="hires", img_alpha=1, size=1.2, alpha=0.85, frameon=False,
        title=[""], linewidth=0, fig=fig, ax=ax, library_id=library,
    )
    fig.set_layout_engine("none")
    main = max(fig.axes, key=lambda a: a.get_position().width)
    main.set_xlim(left, right)
    # Image coordinates run downwards from the top-left corner.
    main.set_ylim(bottom, top)
    main.set_aspect("equal", adjustable="box")

    fig.canvas.draw()
    pixels = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    # Trimmed to the drawn tissue and spots; the figure around them is transparent.
    drawn = pixels[..., 3] > 0
    rows, columns = np.flatnonzero(drawn.any(axis=1)), np.flatnonzero(drawn.any(axis=0))
    return pixels[rows[0]:rows[-1] + 1, columns[0]:columns[-1] + 1]


def block_cells(sections: list[str], top: int, left: int, rows: int, columns: int,
                fill: str) -> dict[str, tuple[int, int]]:
    """Grid cell (row, column) of each of a tumour's sections within its block."""
    order = sorted(sections)
    if fill == "rows":
        return {s: (top + i // columns, left + i % columns) for i, s in enumerate(order)}
    return {s: (top + i % rows, left + i // rows) for i, s in enumerate(order)}


def map_grid(maps: dict[str, np.ndarray], colours: dict):
    """Section maps tiled by tumour, each tumour's maps boxed and named in its colour on the box."""
    figure_width = GRID_COLUMNS + COLOUR_BAR_STRIP
    fig = plt.figure(figsize=(figure_width * MAP_SIZE_INCH, GRID_ROWS * MAP_SIZE_INCH))
    ax = fig.add_axes((0, 0, 1, 1))
    for name, (top, left, rows, columns, fill) in TILING.items():
        own = [section for section in maps if tumour(section) == name]
        extents = []
        for section, (row, column) in block_cells(own, top, left, rows, columns, fill).items():
            # Maps shrink a little towards the block's centre, leaving its margin white.
            x = left + BLOCK_MARGIN + (column - left) * (columns - 2 * BLOCK_MARGIN) / columns
            y = top + BLOCK_MARGIN + (row - top) * (rows - 2 * BLOCK_MARGIN) / rows
            width = (columns - 2 * BLOCK_MARGIN) / columns
            height = (rows - 2 * BLOCK_MARGIN) / rows
            pixels = maps[section]
            scale = MAP_SIDE / max(pixels.shape[:2])
            map_width, map_height = pixels.shape[1] * scale, pixels.shape[0] * scale
            x0, y0 = x + (width - map_width) / 2, y + (height - map_height) / 2
            ax.imshow(pixels, extent=(x0, x0 + map_width, y0 + map_height, y0))
            extents.append((x0, x0 + map_width, y0, y0 + map_height))
        x0, x1 = min(e[0] for e in extents) - BOX_PAD, max(e[1] for e in extents) + BOX_PAD
        y0, y1 = min(e[2] for e in extents) - BOX_PAD, max(e[3] for e in extents) + BOX_PAD
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor="black",
                                   linewidth=2))
        ax.text(x0 + 0.08, y0, name, ha="left", va="bottom", fontsize=16,
                color=colours[name], fontweight="bold")
    ax.set_xlim(0, figure_width)
    ax.set_ylim(GRID_ROWS, 0)
    ax.set_axis_off()

    # Three maps tall, centred on the grid's height.
    colour_bar = ax.inset_axes(((GRID_COLUMNS + 0.05) / figure_width, (GRID_ROWS - 3) / 2 / GRID_ROWS,
                                0.12 / figure_width, 3 / GRID_ROWS))
    fig.colorbar(ScalarMappable(cmap=CMAP), cax=colour_bar, ticks=[0, 1])
    colour_bar.set_yticklabels(["Low", "High"])
    colour_bar.set_ylabel("JBS Axis Score", fontsize=14)
    colour_bar.tick_params(labelsize=12)
    return fig


def tumour_ticks(ax, labels: list[str]) -> None:
    """Tumour names under their slots, with tick marks between the slots."""
    ax.set_xticks(np.arange(1, len(labels)) - 0.5, labels=[])
    ax.set_xticks(range(len(labels)), labels, minor=True)
    ax.tick_params(axis="x", which="minor", length=0)
    ax.set_xlim(-0.5, len(labels) - 0.5)


def section_scores_panel(sections: pd.DataFrame, p: pd.Series, colours: dict):
    fig, ax = plt.subplots(figsize=(8, 3))
    order = tumour_order(sections["Patient"])
    for x, name in enumerate(order):
        own = sections[sections["Patient"] == name]
        sites = own.groupby("Site")["JBS axis score"].mean().sort_values().index
        for offset, site_id in zip(np.linspace(-SITE_SPREAD, SITE_SPREAD, len(sites)), sites):
            scores = own.loc[own["Site"] == site_id, "JBS axis score"]
            ax.plot([x + offset] * len(scores), scores, color="black", linewidth=0.6, zorder=1)
            ax.scatter([x + offset] * len(scores), scores, s=20, color=colours[name],
                       edgecolor="black", linewidth=0.4, zorder=2)

    scores = sections["JBS axis score"]
    pad = 0.05 * (scores.max() - scores.min())
    highest = sections.groupby("Patient")["JBS axis score"].max()
    for x, name in enumerate(order):
        if name in p.index:
            ax.text(x, highest[name] + pad, f"$p$={p[name]:.2f}", ha="center", va="bottom",
                    fontsize=7)
    ax.set_ylim(top=scores.max() + 4 * pad)
    tumour_ticks(ax, order)
    ax.set_xlabel("Patient")
    ax.set_ylabel("Section pseudobulk\nJBS Axis Score")
    fig.tight_layout()
    return fig


def labelled_pair_distances(sites: pd.DataFrame) -> pd.DataFrame:
    """Score distance between every pair of sites, each site named by tumour and site."""
    distances = pair_distances(sites)
    names = (sites["Patient"] + " site " + sites["Site"]).to_numpy()
    return distances.assign(**{"First site": names[distances["first"]],
                               "Second site": names[distances["second"]]})


def dispersion_panel(distances: pd.DataFrame, p: float, title: str):
    fig, ax = plt.subplots(figsize=(2, 4))
    # Cut at the data: the smoothed tails would otherwise run below zero.
    sns.violinplot(data=distances, x="Type", y="Distance", order=["Intra-patient", "Inter-patient"],
                   color=DISPERSION_COLOUR, cut=0, ax=ax)
    for violin in ax.collections:
        violin.set_edgecolor("none")

    low, high = ax.get_ylim()
    pad = 0.03 * (high - low)
    bracket = distances["Distance"].max() + 4 * pad
    ax.plot([0, 0, 1, 1], [bracket - pad, bracket, bracket, bracket - pad], lw=1, c="black")
    ax.text(0.5, bracket + pad, f"$p$ = {p:.2g}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(low, bracket + 4 * pad)

    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Site JBS Axis Score distances")
    fig.tight_layout()
    return fig


def depth_adjusted(spots: pd.DataFrame) -> pd.Series:
    """Spot scores with their linear dependence on log sequencing depth removed within each section.

    Each section keeps its mean score.
    """
    def adjust(section: pd.DataFrame) -> pd.Series:
        depth = np.log10(section["Counts"] + 1)
        slope = np.polyfit(depth, section["JBS axis score"], 1)[0]
        return section["JBS axis score"] - slope * (depth - depth.mean())

    return spots.groupby("Section", group_keys=False)[["Counts", "JBS axis score"]].apply(adjust)


def region_means(spots: pd.DataFrame, reference: str) -> pd.DataFrame:
    """Mean score of each region of each annotated section, and its difference from the reference region."""
    table = (spots.dropna(subset=["Region"])
             .groupby(["Section", "Region"], sort=False)
             .agg(Patient=("Patient", "first"), Spots=("JBS axis score", "size"),
                  **{"Mean spot JBS axis score": ("JBS axis score", "mean"),
                     "Depth-adjusted mean": ("Depth-adjusted score", "mean")})
             .reset_index())
    table = table[table["Spots"] >= MIN_REGION_SPOTS].reset_index(drop=True)

    baseline = table[table["Region"] == reference].set_index("Section")
    for column, difference in [("Mean spot JBS axis score", "Difference"),
                               ("Depth-adjusted mean", "Depth-adjusted difference")]:
        table[difference] = table[column] - table["Section"].map(baseline[column])
    table.loc[table["Region"] == reference, ["Difference", "Depth-adjusted difference"]] = np.nan
    return table


def region_tests(table: pd.DataFrame, compared: list[str]) -> pd.DataFrame:
    """Wilcoxon signed-rank tests of each region's within-section difference from the reference."""
    rows = []
    for region in compared:
        paired = table[(table["Region"] == region) & table["Difference"].notna()]
        rows.append({
            "Region": region, "Sections": len(paired), "Tumours": paired["Patient"].nunique(),
            "Median difference": paired["Difference"].median(),
            "p": wilcoxon(paired["Difference"]).pvalue,
            "Median depth-adjusted difference": paired["Depth-adjusted difference"].median(),
            "p depth-adjusted": wilcoxon(paired["Depth-adjusted difference"]).pvalue,
        })
    tests = pd.DataFrame(rows)
    tests["p BH"] = multipletests(tests["p"], method="fdr_bh")[1]
    tests["p depth-adjusted BH"] = multipletests(tests["p depth-adjusted"], method="fdr_bh")[1]
    return tests


def region_panel(table: pd.DataFrame, sections: pd.Series, regions: list[str], colours: dict):
    """One point per region per annotated section, regions from periphery to core within each tumour."""
    fig, ax = plt.subplots(figsize=(7, 3))
    order = tumour_order(table["Patient"])
    places = dict(zip(regions, np.linspace(-SITE_SPREAD, SITE_SPREAD, len(regions))))
    for x, name in enumerate(order):
        own = table[table["Patient"] == name]
        names = sorted(own["Section"].unique())
        for nudge, section in zip(np.linspace(-SECTION_SPREAD, SECTION_SPREAD, len(names)), names):
            points = own[own["Section"] == section].set_index("Region").reindex(regions).dropna(
                subset=["Mean spot JBS axis score"])
            positions = [x + places[region] + nudge for region in points.index]
            ax.plot(positions, points["Mean spot JBS axis score"], color="black", linewidth=0.5,
                    alpha=0.5, zorder=1)
            ax.scatter(positions, points["Mean spot JBS axis score"],
                       c=[colours[region] for region in points.index], s=18, edgecolor="black",
                       linewidth=0.4, zorder=2)
    handles = [plt.Line2D([], [], marker="o", linestyle="", markerfacecolor=colours[region],
                          markeredgecolor="black", markeredgewidth=0.4, markersize=5, label=region)
               for region in regions]
    ax.legend(handles=handles, title="Spots' IvyGAP region", fontsize=7, title_fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    tumour_ticks(ax, [f"{name}\n{sections[name]} sections" for name in order])
    ax.set_xlabel("Patient (sections with IvyGAP annotation)")
    ax.set_ylabel("Mean JBS Axis Score of a\nsection's spots in the region")
    fig.tight_layout()
    return fig


def region_difference_panel(table: pd.DataFrame, tests: pd.DataFrame, reference: str,
                            colours: dict):
    fig, ax = plt.subplots(figsize=(3, 3))
    compared = list(tests["Region"])
    sns.stripplot(data=table.dropna(subset=["Difference"]), x="Region", y="Difference",
                  order=compared, hue="Region", palette=colours, legend=False, size=4,
                  edgecolor="black", linewidth=0.5, ax=ax)
    ax.axhline(0, color="grey", linestyle="--", linewidth=1)
    top = table["Difference"].max()
    for x, row in tests.iterrows():
        ax.text(x, top * 1.05, f"$p$={row['p BH']:.1e}", ha="center", va="bottom", fontsize=7)
    ax.set_ylim(top=top * 1.2)
    ax.set_xlabel("")
    ax.set_ylabel(f"Difference from {reference.lower()}")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    return fig


def main() -> None:
    entry = cohort("spatial", DATASET)
    symbols = config()["gene_symbols"]
    genes = axis_genes(entry, symbols)
    measured = measured_symbols(symbols)
    regions = list(entry["histology_regions"])

    set_paper_theme()
    maps, spots, sections, failed = {}, [], [], []
    for path in section_paths(entry):
        section = read_section(path)
        data = features(section, entry, entry["expression_features"])
        counts = np.asarray(data.X.sum(axis=1)).ravel()
        if np.median(counts) < MIN_MEDIAN_COUNTS:
            failed.append(f"{path.stem} ({np.median(counts):.0f})")
            continue
        data.obs["JBS axis score"] = spot_scores(data, genes, measured, entry["pseudocount"])
        maps[path.stem] = section_map(data)
        spots.append(pd.DataFrame({
            "Section": path.stem,
            "Patient": tumour(path.stem),
            "Spot": data.obs_names,
            "x": data.obsm["spatial"][:, 0],
            "y": data.obsm["spatial"][:, 1],
            "JBS axis score": data.obs["JBS axis score"].to_numpy(),
            "Counts": counts,
            "Region": spot_regions(features(section, entry, entry["histology_features"]),
                                   entry["histology_regions"]).to_numpy(),
        }))
        sections.append({
            "Section": path.stem,
            "Patient": tumour(path.stem),
            "Site": site(path.stem),
            "Spots": data.n_obs,
            "JBS axis score": section_score(data, genes, measured, entry["pseudocount"]),
        })
    spots = pd.concat(spots, ignore_index=True)
    spots["Depth-adjusted score"] = depth_adjusted(spots)
    sections = pd.DataFrame(sections)

    sites = sections.groupby(["Patient", "Site"], sort=False)["JBS axis score"].mean().reset_index()
    distances = labelled_pair_distances(sites)
    p_dispersion = within_patient_p(sites, distances, N_PERMUTATIONS, SEED)
    p_tumour = patient_clustering(sites, MIN_SITES, N_PERMUTATIONS, SEED)
    title = f"{sites['Patient'].nunique()} tumours, {len(sites)} sites"

    reference = entry["reference_region"]
    by_region = region_means(spots, reference)
    tests = region_tests(by_region, [region for region in regions if region != reference])
    annotated = spots.dropna(subset=["Region"]).groupby("Patient")["Section"].nunique()
    print(f"Failed sections left out (median counts): {', '.join(failed)}")
    print(f"{len(maps)} sections, {len(spots)} spots; {annotated.sum()} sections from "
          f"{len(annotated)} tumours with histology")
    print(tests.round(4).to_string(index=False))
    pairs = distances.groupby("Type")["Distance"].agg(["size", "mean"])
    print(f"{title}: intra {pairs.loc['Intra-patient', 'mean']:.1f} "
          f"({pairs.loc['Intra-patient', 'size']} pairs), inter {pairs.loc['Inter-patient', 'mean']:.1f} "
          f"({pairs.loc['Inter-patient', 'size']}); p = {p_dispersion:.2e}")
    print("Per-tumour clustering of sites:")
    print(pd.DataFrame({"p": p_tumour, "p (BH)": multipletests(p_tumour, method="fdr_bh")[1]})
          .round(4).to_string())

    tumours, region_palette = palette("deJong_tumour"), palette("histology_region")
    panels = {
        "visium_axis_score_maps": (map_grid(maps, tumours),
                                   spots[["Section", "Patient", "Spot", "x", "y", "JBS axis score"]]),
        "visium_section_axis_score": (section_scores_panel(sections, p_tumour, tumours),
                                      sections.assign(p=sections["Patient"].map(p_tumour))),
        "visium_axis_score_dispersion": (
            dispersion_panel(distances, p_dispersion, title),
            distances.drop(columns=["first", "second"]).assign(p=p_dispersion)),
        "visium_axis_score_by_region": (region_panel(by_region, annotated, regions, region_palette),
                                        by_region.drop(columns=["Difference", "Depth-adjusted difference"])),
        "visium_axis_score_region_difference": (
            region_difference_panel(by_region, tests, reference, region_palette),
            by_region.dropna(subset=["Difference"])),
    }
    for name, (fig, plotted) in panels.items():
        save_panel(fig, plotted, name, dpi=GRID_DPI if name == "visium_axis_score_maps" else 300)
        plt.close(fig)


if __name__ == "__main__":
    main()
