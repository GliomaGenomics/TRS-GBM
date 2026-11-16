import re
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parent.parent
DATASETS = REPO_ROOT / "config" / "datasets.yaml"
FIGURES = REPO_ROOT / "figures"


def config() -> dict:
    with open(DATASETS) as handle:
        return yaml.safe_load(handle) or {}


def resolve(path: str) -> Path:
    path = Path(path)
    return path if path.is_absolute() else REPO_ROOT / path


def set_nature_style() -> None:
    """Nature figure style, as in common.R.

    A panel is drawn at its final printed size, with Arial text of 5 to 7 pt
    (titles and axis titles 7, everything else 6) and lines of at least 0.25 pt.
    Text stays editable in the PDF.
    """
    mpl.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "font.size": 6,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.5,
        "axes.labelpad": 2,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "xtick.major.pad": 1.5,
        "ytick.major.pad": 1.5,
        "lines.linewidth": 0.75,
        "patch.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def cohort(modality: str, name: str) -> dict:
    return config()["datasets"][modality][name]


def natural_order(names) -> list:
    """Identifiers sorted with their numbers compared as numbers (P2 before P10)."""
    return sorted(names, key=lambda name: [int(part) if part.isdigit() else part
                                           for part in re.split(r"(\d+)", str(name))])


def selected(obs: pd.DataFrame, filters: dict) -> pd.Series:
    mask = pd.Series(True, index=obs.index)
    for column, accepted in filters.items():
        mask &= obs[column].astype(object).isin(accepted)
    return mask


def bulk_counts(entry: dict, data) -> pd.DataFrame:
    """A bulk h5ad as the (genes x samples) count matrix jbs_axis_score expects."""
    identifiers = data.var[entry["gene_id_column"]].astype(object)
    named = identifiers.notna().to_numpy()
    return pd.DataFrame(
        data.X[:, named].T,
        # Ensembl IDs are truncated so versioned and unversioned IDs match.
        index=pd.Index(identifiers[named].str.split(".").str[0], name="gene"),
        columns=data.obs_names,
    )


def trs_table(name: str) -> pd.DataFrame:
    """One row per patient of a bulk cohort: its TRS and its primary tumour's metadata.

    TRS is read from the sample metadata, where trs_column.py writes it. The
    clinical record is carried on the primary sample, which is also the row the
    purity threshold and the inclusion flag are read from.
    """
    entry = cohort("bulk", name)
    obs = ad.read_h5ad(resolve(entry["path"]), backed="r").obs

    patients = obs[entry["patient_column"]].astype(object)
    included = (
        selected(obs, entry["filters"])
        & obs[entry["include_column"]].to_numpy()
        & (obs[entry["pair_purity_column"]] >= entry["min_purity"])
    )

    primary = obs[included].set_index(patients[included]).rename_axis("Patient")
    # The file's own Cohort column names the pooled cohort; the label names the
    # entry the samples were selected by, which is what the figures group on.
    primary["Cohort"] = entry["label"]
    return primary


def stage_samples(names: list[str], stage: str) -> pd.DataFrame:
    """Metadata of one stage's samples ("primary" or "recurrent") of the included
    pairs of bulk cohorts sharing one file, in the order the file holds them."""
    entries = [cohort("bulk", name) for name in names]
    obs = ad.read_h5ad(resolve(entries[0]["path"]), backed="r").obs
    keep = pd.Series(False, index=obs.index)
    for name, entry in zip(names, entries):
        keep |= (
            selected(obs, entry["filters"])
            & obs[entry["patient_column"]].isin(trs_table(name).index)
            & (obs[entry["stage_column"]] == entry[f"{stage}_value"])
        )
    return obs[keep]


def fpkm_axis_scores(entry: dict) -> pd.Series:
    """Every sample of the pipeline FPKM matrix placed on the JBS axis.

    Genes the pipeline left unquantified count as unexpressed rather than missing.
    """
    loadings = pd.read_csv(resolve(entry["axis_loadings"]), index_col=0).iloc[:, 0]
    expression = pd.read_csv(resolve(entry["fpkm"]), index_col=0)
    genes = expression.columns.intersection(loadings.index)
    logged = np.log2(expression.loc[:, genes].fillna(0) + entry["pseudocount"])
    return logged.dot(loadings.loc[genes])


PNMES_COUNTS_PER_SAMPLE = 10_000


def pnmes_scores(counts: pd.DataFrame, loadings: pd.Series) -> pd.Series:
    """Place every sample (row) of a count matrix on the PN:MES component.

    Counts are scaled to PNMES_COUNTS_PER_SAMPLE and log10(1 + x) transformed.
    Library sizes are taken over every gene in the matrix, unquantified genes
    counting as zero, before the matrix is cut to the genes with a loading.
    """
    counts = counts.fillna(0)
    normalised = np.log10(1 + counts.div(counts.sum(axis=1), axis=0) * PNMES_COUNTS_PER_SAMPLE)
    genes = normalised.columns.intersection(loadings.index)
    return normalised.loc[:, genes].dot(loadings.loc[genes])


def pseudobulk_axis_scores(data: ad.AnnData, pseudocount: float) -> pd.Series:
    """Place every pseudobulk on the JBS axis.

    The axis loadings and gene lengths travel in .var of the single-cell data, so
    the pseudobulks are scored against those rather than the gene annotation the
    bulk cohorts are matched to. Genes with no recorded length cannot be FPKM
    normalised and are dropped first, which is also what fixes the library size.
    """
    measured = data.var["gene_length"].notna().to_numpy()
    counts = np.asarray(data.X)[:, measured].astype(np.float32)
    lengths_kb = data.var["gene_length"].to_numpy(dtype=np.float32)[measured] / 1e3

    library_sizes = counts.sum(axis=1) / 1e6
    log2_fpkm = np.log2(counts / library_sizes[:, None] / lengths_kb[None, :] + pseudocount)

    loadings = data.var["JRS loading"].to_numpy(dtype=np.float32)[measured]
    scored = ~np.isnan(loadings)
    return pd.Series(log2_fpkm[:, scored].dot(loadings[scored]), index=data.obs_names)


def cell_axis_scores(data: ad.AnnData, pseudocount: float) -> pd.Series:
    """Place every cell of a single-cell dataset on the JBS axis.

    The axis loadings and gene lengths travel in .var of the single-cell data, so
    the cells are scored against those rather than the gene annotation the bulk
    cohorts are matched to. Genes with no recorded length cannot be FPKM
    normalised and are dropped first, which is also what fixes each cell's
    library size.

    A cell by gene matrix is too large to log transform densely, so the score is
    summed over the stored counts alone and offset by the score a cell of nothing
    but zeros would take, which is the same for every cell.
    """
    measured = data.var["gene_length"].notna().to_numpy()
    counts = sparse.csr_matrix(data.X)[:, measured]
    loadings = data.var["JRS loading"].to_numpy(dtype=np.float32)[measured]

    scored = ~np.isnan(loadings)
    weights = loadings[scored]
    lengths_kb = data.var["gene_length"].to_numpy(dtype=np.float32)[measured][scored] / 1e3
    library_sizes = np.asarray(counts.sum(axis=1)).ravel() / 1e6

    logged = counts[:, scored].tocoo()
    fpkm = logged.data / (library_sizes[logged.row] * lengths_kb[logged.col])
    logged.data = np.log2(fpkm + pseudocount) - np.log2(pseudocount)

    zeros = np.log2(pseudocount) * weights.sum()
    return pd.Series(logged @ weights + zeros, index=data.obs_names)


def compositions(obs: pd.DataFrame, entry: dict) -> pd.DataFrame:
    """Turn the cell counts a pseudobulk was summed from into proportions.

    Three sets of proportions are added, each summing to one within a pseudobulk:
    comp_<cell type> and comp_pop_<population> over all of its cells, and
    malignant_comp_<state> over its malignant cells alone.
    """
    obs = obs.copy()
    cell_types = [c for c in obs.columns
                  if c.startswith("n_") and not c.startswith("n_malignant_")]
    states = [c for c in obs.columns if c.startswith("n_malignant_")]

    for column in states:
        state = column.removeprefix("n_malignant_")
        obs[f"malignant_comp_{state}"] = (obs[column] / obs["n_Malignant"]).fillna(0)

    for column in cell_types:
        population = entry["cell_populations"][column.removeprefix("n_")]
        pooled = f"n_pop_{population}"
        obs[pooled] = obs[column] if pooled not in obs.columns else obs[pooled] + obs[column]

    counted = [c for c in obs.columns if c.startswith("n_")]
    proportions = {c.replace("n_", "comp_", 1): (obs[c] / obs["total_cells"]).fillna(0)
                   for c in counted}

    obs["Purity"] = (obs["n_Malignant"] / obs["total_cells"]).fillna(0)
    return obs.assign(**proportions)


def pseudobulk_table(name: str) -> pd.DataFrame:
    """One row per pseudobulk of a single-cell dataset: axis score and composition."""
    entry = cohort("single_cell", name)
    data = ad.read_h5ad(resolve(entry["pseudobulk"]))
    obs = compositions(data.obs, entry)
    obs["Axis Score"] = pseudobulk_axis_scores(data, entry["pseudocount"])
    return obs


def pair_distances(table: pd.DataFrame) -> pd.DataFrame:
    """Axis score distance between every pair of samples, and whether they share a patient."""
    pairs = np.array(list(combinations(range(len(table)), 2)))
    patients = table["Patient"].to_numpy()
    scores = table["JBS axis score"].to_numpy(float)
    within = patients[pairs[:, 0]] == patients[pairs[:, 1]]
    return pd.DataFrame({
        "first": pairs[:, 0],
        "second": pairs[:, 1],
        "Type": np.where(within, "Intra-patient", "Inter-patient"),
        "Distance": np.abs(scores[pairs[:, 0]] - scores[pairs[:, 1]]),
    })


def within_patient_p(table: pd.DataFrame, distances: pd.DataFrame, n_permutations: int,
                     seed: int) -> float:
    """One-sided permutation p for within-patient distances being smaller than between.

    Patient labels are shuffled across samples, so each permutation keeps the
    distances themselves and their dependence on shared samples.
    """
    first, second = distances["first"].to_numpy(), distances["second"].to_numpy()
    values = distances["Distance"].to_numpy()
    patients = table["Patient"].to_numpy()

    def statistic(labels: np.ndarray) -> float | None:
        within = labels[first] == labels[second]
        if not within.any() or within.all():
            return None
        return values[within].mean() - values[~within].mean()

    observed = statistic(patients)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_permutations):
        permuted = statistic(rng.permutation(patients))
        if permuted is not None:
            count += permuted <= observed
    return (1 + count) / (n_permutations + 1)


def mean_pairwise_distance(scores: np.ndarray) -> float:
    return np.fromiter((abs(a - b) for a, b in combinations(scores, 2)), dtype=float).mean()


def patient_clustering(table: pd.DataFrame, min_samples: int, n_permutations: int,
                       seed: int) -> pd.Series:
    """Per-patient permutation p: are the patient's samples closer together than chance?

    Only patients with at least min_samples samples are tested. The null draws as
    many samples as the patient has, without replacement, from the samples of
    every patient tested (its own included).
    """
    groups = table.groupby("Patient")["JBS axis score"].apply(lambda s: s.to_numpy(float))
    groups = groups[groups.map(len) >= min_samples]
    pool = np.concatenate(groups.tolist())

    rng = np.random.default_rng(seed)
    p = {}
    for patient, scores in groups.items():
        observed = mean_pairwise_distance(scores)
        null = np.array([mean_pairwise_distance(rng.choice(pool, size=len(scores), replace=False))
                         for _ in range(n_permutations)])
        p[patient] = (1 + np.sum(null <= observed)) / (n_permutations + 1)
    return pd.Series(p, name="p")


PANELS = FIGURES / "panels"
SOURCE_DATA = FIGURES / "source_data"

PALETTE = {
    "Proneural": "#ef466f",
    "Mesenchymal": "#f3c962",
    "Classical": "#248ab2",
    "Up": "#FFBA49",
    "Down": "#63CCCA",
    "Discovery": "#50C878",
    "EORTC": "#FE7F9C",
    "GLASS": "#0197f6",
    "Primary": "#29AB87",
    "Recurrent": "#E52B50",
}


# Cell types and malignant states, coloured by the compartment they belong to.
CELL_COLOURS = {
    "Malignant": "#e07271",
    "Brain": "#449ddf",
    "Immune": "#006d2c",
    "Other": "lightgrey",
    "Unknown": "lightgrey",

    "Myeloid": "#f3c962",
    "Lymphocyte": "#006d2c",
    "TAM": "#78c679",

    "Oligodendrocyte": "#8c6bb1",
    "Neuron": "#31b2d2",
    "Excitatory neuron": "#31b2d2",
    "Inhibitory neuron": "#31b2d2",
    "Astrocyte": "#88419d",
    "OPC": "#388ae0",

    "Endothel": "goldenrod",
    "Pericyte": "brown",

    "MES1": "#f3cc0b",
    "MES2": "#FFBA49",
    "NPC1": "#5f5feb",
    "NPC2": "#450dbb",
}

# The Nomura malignant states are drawn from a rainbow ramp, in this order, with
# the unresolved cells left grey.
RAINBOW_STATES = ["Stress", "Hypoxia", "MES", "Unresolved", "GPC", "AC", "Cilia",
                  "Neuron", "OPC", "NPC"]


def cell_colours() -> dict:
    palette = dict(CELL_COLOURS)
    ramp = sns.color_palette("rainbow_r", len(RAINBOW_STATES), desat=0.8)
    palette.update(dict(zip(RAINBOW_STATES, ramp)))
    palette["Unresolved"] = "lightgrey"
    return palette


def set_paper_theme() -> None:
    """The seaborn theme the single-cell panels were drawn with."""
    sns.set_theme(
        style="ticks",
        context="paper",
        font_scale=1,
        font="sans-serif",
        rc={
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "savefig.transparent": False,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        },
    )


def set_correlate_style() -> None:
    """Larger type for the TRS correlate panels (clinical variables, cell types)."""
    mpl.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })


def correlation_label(name: str, r: float, p: float) -> str:
    """Legend entry for a regression panel, emboldening a significant p-value."""
    formatted = rf"$\mathbf{{p={p:.3f}}}$" if p < 0.05 else f"p={p:.3f}"
    return f"{name}: R={r:.2f}, " + formatted


def save_panel(fig, table: pd.DataFrame, name: str, bbox_inches: str | None = "tight",
               dpi: int = 300) -> None:
    """Write a panel and the tidy table it was plotted from.

    A panel drawn at its final size passes bbox_inches=None, so the files keep
    that size rather than being cropped to their contents. A panel drawn several
    times larger than it will be printed may pass a lower dpi.
    """
    (PANELS / name).parent.mkdir(parents=True, exist_ok=True)
    (SOURCE_DATA / name).parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(PANELS / f"{name}.{suffix}", dpi=dpi, bbox_inches=bbox_inches)
    table.to_csv(SOURCE_DATA / f"{name}.csv", index=False)
    print(f"wrote {name}: {len(table)} rows")
