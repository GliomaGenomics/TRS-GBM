"""Panels of the primary-tumour variant comparison between Up and Down responders.

From the tables responder_variants.R writes: mutational burden by responder subtype
(mutational_burden_by_responder_<cohort>), log10 odds ratios of the 20 genes with
the lowest Fisher p, capped at ±3 (variant_responder_association_<cohort>), and the
SRSF11 variants on the protein (srsf11_variants).

Example:
    pixi run responder_variant_panels
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.lines import Line2D

from common import cohort, config, resolve, save_panel, set_paper_theme

COHORTS = ("EORTC", "Discovery_Stead")
SUBTYPES = {"D": "Down", "U": "Up"}

BURDEN_CATEGORIES = {"total_perMB": "Non-synonymous", "total_perMB_HIGH": "High Impact"}
VIOLIN_WIDTH = 0.9
# Violin tails run three bandwidths past the extreme samples.
VIOLIN_CUT = 3
BOX_WIDTH = 0.1
JITTER = 0.2
JITTER_SEED = 20

TOP_GENES = 20
LOG_OR_CAP = 3

VARIANT_LABELS = {"Missense_Mutation": "Missense", "Nonsense_Mutation": "Nonsense"}


def palette(group: str) -> dict:
    with open(resolve("config/palette.yaml")) as handle:
        return yaml.safe_load(handle)[group]


def nrd0_bandwidth(values: np.ndarray) -> float:
    """Silverman's rule of thumb as R's bw.nrd0 computes it."""
    spread = min(values.std(ddof=1), np.subtract(*np.percentile(values, [75, 25])) / 1.34)
    if spread == 0:
        spread = values.std(ddof=1) or abs(values[0]) or 1
    return 0.9 * spread * len(values) ** -0.2


def violin_density(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    bandwidth = nrd0_bandwidth(values)
    grid = np.linspace(values.min() - VIOLIN_CUT * bandwidth,
                       values.max() + VIOLIN_CUT * bandwidth, 512)
    kernels = np.exp(-0.5 * ((grid[:, None] - values[None, :]) / bandwidth) ** 2)
    density = kernels.sum(axis=1) / (len(values) * bandwidth * np.sqrt(2 * np.pi))
    return grid, density


def burden_axis(ax, burden: pd.DataFrame, column: str, colours: dict, rng) -> None:
    """One impact class: violins scaled to equal area, box, and jittered samples."""
    groups = {subtype: np.log10(burden.loc[burden["Responder"] == subtype, column].to_numpy())
              for subtype in SUBTYPES}
    densities = {subtype: violin_density(values) for subtype, values in groups.items()}
    peak = max(density.max() for _, density in densities.values())

    for position, (subtype, values) in enumerate(groups.items()):
        colour = colours[SUBTYPES[subtype]]
        grid, density = densities[subtype]
        half_width = density / peak * VIOLIN_WIDTH / 2
        ax.fill_betweenx(grid, position - half_width, position + half_width,
                         color=colour, alpha=0.4, linewidth=0, zorder=1)
        ax.scatter(position + rng.uniform(-JITTER, JITTER, len(values)), values,
                   s=4, color="black", alpha=0.5, linewidths=0, zorder=2)
        ax.boxplot(values, positions=[position], widths=BOX_WIDTH, whis=1.5,
                   showfliers=False, patch_artist=True, zorder=3,
                   boxprops={"facecolor": colour, "edgecolor": "black", "linewidth": 0.8},
                   medianprops={"color": "black", "linewidth": 1},
                   whiskerprops={"color": "black", "linewidth": 0.8},
                   capprops={"linewidth": 0})

    ax.set_xticks(range(len(SUBTYPES)), list(SUBTYPES))
    ax.set_xlim(-0.6, len(SUBTYPES) - 0.4)
    ax.set_title(BURDEN_CATEGORIES[column], fontsize=9)


def burden_panel(burden: pd.DataFrame, label: str, colours: dict):
    fig, axes = plt.subplots(1, len(BURDEN_CATEGORIES), figsize=(4, 4), sharey=True,
                             layout="constrained")
    rng = np.random.default_rng(JITTER_SEED)
    for ax, column in zip(axes, BURDEN_CATEGORIES):
        burden_axis(ax, burden, column, colours, rng)
    low, high = axes[0].get_ylim()
    exponents = np.arange(np.ceil(low), np.floor(high) + 1).astype(int)
    axes[0].set_yticks(exponents, [f"$10^{{{k}}}$" for k in exponents])
    axes[0].set_ylabel("Mutations per Mb (log scale)")
    fig.supxlabel("Responder Subtype", fontsize=9)
    fig.suptitle(f"{label} mutational burden", fontsize=10)

    table = burden.melt(id_vars=["Tumor_Sample_Barcode", "Responder"],
                        value_vars=list(BURDEN_CATEGORIES),
                        var_name="Category", value_name="Mutations per Mb")
    table["Category"] = table["Category"].map(BURDEN_CATEGORIES)
    return fig, table


def top_associations(results: pd.DataFrame) -> pd.DataFrame:
    shown = results.sort_values("pval", kind="stable").head(TOP_GENES).copy()
    shown["direction"] = np.where(shown["or"] < 1, "Down", "Up")
    with np.errstate(divide="ignore"):
        shown["logOR"] = np.log10(shown["or"])
        shown["logCI_low"] = np.log10(shown["ci.low"])
        shown["logCI_high"] = np.log10(shown["ci.up"])
    shown["logOR_capped"] = shown["logOR"].clip(-LOG_OR_CAP, LOG_OR_CAP)
    shown["ci_low_capped"] = shown["logCI_low"].clip(lower=-LOG_OR_CAP)
    shown["ci_high_capped"] = shown["logCI_high"].clip(upper=LOG_OR_CAP)
    shown["capped"] = shown["logOR"].abs() > LOG_OR_CAP
    return shown


def association_panel(shown: pd.DataFrame, label: str, colours: dict):
    fig, ax = plt.subplots(figsize=(5, 5))
    rows = np.arange(len(shown))[::-1]
    ax.errorbar(shown["logOR_capped"], rows,
                xerr=[shown["logOR_capped"] - shown["ci_low_capped"],
                      shown["ci_high_capped"] - shown["logOR_capped"]],
                fmt="none", ecolor="0.25", elinewidth=0.8, capsize=2.5, zorder=1)
    for capped, marker, size in ((False, "o", 40), (True, "D", 32)):
        drawn = (shown["capped"] == capped).to_numpy()
        ax.scatter(shown["logOR_capped"][drawn], rows[drawn], marker=marker, s=size,
                   color=[colours[direction] for direction in shown["direction"][drawn]],
                   linewidths=0, zorder=2)
    ax.axvline(0, color="0.5", linestyle="--", linewidth=0.8, zorder=0)
    ax.set_yticks(rows, shown["Hugo_Symbol"])
    ax.set_ylim(-0.7, len(shown) - 0.3)
    ax.set_xlabel(r"Log$_{10}$ Odds Ratio")
    ax.set_ylabel(f"Top {TOP_GENES} genes")
    ax.set_title(f"{label} - Odds of being Up given mutation", fontsize=10)
    fig.tight_layout()
    return fig, shown


def protein_track(ax, variants: pd.DataFrame, protein: pd.DataFrame, baseline: float,
                  label: str, border: str, variant_colours: dict) -> None:
    length = protein["aa.length"].iloc[0]
    ax.add_patch(plt.Rectangle((0, baseline - 0.15), length, 0.3, facecolor="0.75",
                               edgecolor=border, linewidth=1.5, zorder=1))
    for domain in protein.itertuples():
        ax.add_patch(plt.Rectangle((domain.Start, baseline - 0.22), domain.End - domain.Start,
                                   0.44, facecolor="white", edgecolor="black", linewidth=0.8,
                                   zorder=2))
        ax.text((domain.Start + domain.End) / 2, baseline, domain.Label, ha="center",
                va="center", fontsize=8, zorder=3)
    ax.text(length - 5, baseline, label, ha="right", va="center", fontsize=8, zorder=3)

    # Variants at the same residue stack on one stem.
    counts = variants.groupby(["Protein_position", "Variant_Classification"]).size()
    for (position, classification), count in counts.items():
        top = baseline + 0.15 + 0.35 * count
        ax.plot([position, position], [baseline + 0.15, top], color="0.4", linewidth=0.8,
                zorder=0)
        ax.scatter(position, top, s=30, color=variant_colours[classification],
                   edgecolors="black", linewidths=0.4, zorder=4)


def gene_panel(variants: pd.DataFrame, protein: pd.DataFrame, colours: dict,
               variant_colours: dict):
    fig, ax = plt.subplots(figsize=(6, 3))
    for baseline, subtype in ((1.6, "Up"), (0, "Down")):
        protein_track(ax, variants[variants["Responder"] == subtype], protein, baseline,
                      f"{subtype} responders", colours[subtype], variant_colours)

    length = protein["aa.length"].iloc[0]
    ax.set_xlim(0, length)
    ax.set_xticks([*range(0, length - 49, 100), length])
    ax.set_ylim(-0.4, 2.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Amino acid position")
    ax.set_title(protein["Gene"].iloc[0], fontsize=10, fontstyle="italic")
    shown = [name for name in VARIANT_LABELS if name in set(variants["Variant_Classification"])]
    handles = [Line2D([], [], marker="o", linestyle="", markersize=6,
                      markerfacecolor=variant_colours[name], markeredgecolor="black",
                      markeredgewidth=0.4, label=VARIANT_LABELS[name]) for name in shown]
    ax.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1, 1),
              fontsize=8)
    fig.tight_layout()
    return fig, variants


def main() -> None:
    set_paper_theme()
    colours = palette("responder")
    variant_colours = palette("variant_class")

    for name in COHORTS:
        entry = cohort("bulk", name)
        burden = pd.read_csv(resolve(entry["mutational_burden"]))
        fig, table = burden_panel(burden, entry["label"], colours)
        save_panel(fig, table, f"mutational_burden_by_responder_{name}")
        plt.close(fig)

        results = pd.read_csv(resolve(entry["responder_variant_association"]))
        fig, table = association_panel(top_associations(results), entry["label"], colours)
        save_panel(fig, table, f"variant_responder_association_{name}")
        plt.close(fig)

    settings = config()
    variants = pd.read_csv(resolve(settings["srsf11_variants"]))
    protein = pd.read_csv(resolve(settings["srsf11_protein"]))
    fig, table = gene_panel(variants, protein, colours, variant_colours)
    save_panel(fig, table, "srsf11_variants")
    plt.close(fig)


if __name__ == "__main__":
    main()
