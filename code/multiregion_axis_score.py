"""JBS axis scores across the regions of single tumours (Mathur 2025).

94 multi-region samples of ten glioblastomas, scored as log2(TPM + 0.01) · axis
loadings, matched by gene symbol. Panels:

    multiregion_axis_score                     scores per patient by contrast
                                               enhancement, with a permutation p for
                                               the patient's samples lying closer
                                               together than random samples
    multiregion_axis_score_dispersion          within- against between-patient
                                               score distances
    multiregion_axis_score_vs_centre_distance  score against relative distance from
                                               the tumour centroid, Pearson's r

Example:
    pixi run multiregion_axis_score
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

from common import (cohort, config, pair_distances, patient_clustering, resolve, save_panel,
                    set_paper_theme, within_patient_p)

DATASET = "mathur_2025"
N_PERMUTATIONS = 10_000
SEED = 11
# Fewest samples a patient needs for its own clustering test.
MIN_SAMPLES = 3
# Leaves room above the highest points for the p-value labels and the legend.
SCORE_AXIS_TOP = 183
DISPERSION_COLOURS = {"All samples": "slategrey", "CEL samples": "#73c05a"}
CEL_LABEL_COLOUR = "forestgreen"


def contrast_enhancement_colours() -> dict:
    with open(resolve("config/palette.yaml")) as handle:
        return yaml.safe_load(handle)["contrast_enhancement"]


def read_tpm(entry: dict) -> pd.DataFrame:
    """Expression as TPM, samples x genes."""
    logged = pd.read_excel(resolve(entry["expression"]), sheet_name=entry["expression_sheet"],
                           index_col=0)
    return (np.power(2, logged) - 1).T


def read_samples(entry: dict) -> pd.DataFrame:
    return pd.read_excel(resolve(entry["metadata"]), sheet_name=entry["sample_sheet"],
                         skiprows=1).set_index(entry["sample_column"])


def symbol_loadings(path: str, symbols: dict) -> pd.Series:
    """JBS axis loadings indexed by gene symbol; genes with no symbol are dropped."""
    loadings = pd.read_csv(resolve(path), index_col=0).iloc[:, 0]
    names = pd.read_csv(resolve(symbols["path"]), index_col=0)[symbols["column"]]
    loadings.index = names.loc[loadings.index].to_numpy()
    return loadings[loadings.index.notna()]


def axis_scores(tpm: pd.DataFrame, loadings: pd.Series, pseudocount: float) -> pd.Series:
    genes = loadings.index.intersection(tpm.columns)
    return np.log2(tpm.loc[:, genes] + pseudocount).dot(loadings.loc[genes])


def sample_table(entry: dict, loadings: pd.Series) -> pd.DataFrame:
    """One row per sample in both supplementary tables, in the expression table's order."""
    tpm = read_tpm(entry)
    samples = read_samples(entry)
    ids = tpm.index.intersection(samples.index)
    return pd.DataFrame({
        "Sample": ids,
        "Patient": samples.loc[ids, entry["patient_column"]].to_numpy(),
        "CEL": samples.loc[ids, entry["contrast_enhancement_column"]].to_numpy(),
        "Relative distance from centroid": samples.loc[ids, entry["centre_distance_column"]].to_numpy(),
        "JBS axis score": axis_scores(tpm.loc[ids], loadings, entry["pseudocount"]).to_numpy(),
    })


def patient_panel(table: pd.DataFrame, p_all: pd.Series, p_cel: pd.Series, colours: dict):
    fig, ax = plt.subplots(figsize=(6, 4))
    order = sorted(table["Patient"].unique())
    sns.stripplot(data=table, x="Patient", y="JBS axis score", order=order, hue="CEL",
                  hue_order=list(colours), palette=colours, edgecolor="black", linewidth=.2,
                  alpha=.8, ax=ax)

    scores = table["JBS axis score"]
    pad = 0.05 * (scores.max() - scores.min())
    highest = table.groupby("Patient")["JBS axis score"].max()
    for x, patient in enumerate(order):
        y = highest[patient] + pad
        if patient in p_all.index:
            ax.text(x, y, f"$p$={p_all[patient]:.2f}", ha="center", va="bottom", fontsize=7,
                    color="black")
        if patient in p_cel.index:
            ax.text(x, y + pad, f"$p_{{CEL}}$={p_cel[patient]:.2f}", ha="center", va="bottom",
                    fontsize=7, color=CEL_LABEL_COLOUR)

    ax.set_ylim(top=SCORE_AXIS_TOP)
    ax.legend(title="Contrast Enhancing", loc="upper center", ncol=3)
    ax.set_ylabel("JBS Axis Score")
    fig.tight_layout()
    return fig


def dispersion_panel(distances: pd.DataFrame, p: dict):
    fig, axes = plt.subplots(1, 2, figsize=(4, 4))
    for ax, (subset, colour) in zip(axes, DISPERSION_COLOURS.items()):
        data = distances[distances["Samples"] == subset]
        sns.violinplot(data=data, x="Type", y="Distance", order=["Intra-patient", "Inter-patient"],
                       color=colour, ax=ax)

        low, high = ax.get_ylim()
        pad = 0.03 * (high - low)
        bracket = data["Distance"].max() + 4 * pad
        ax.plot([0, 0, 1, 1], [bracket - pad, bracket, bracket, bracket - pad], lw=1, c="black")
        ax.text(0.5, bracket + pad, f"$p$ = {p[subset]:.2e}", ha="center", va="bottom", fontsize=8)
        ax.set_ylim(low, bracket + 4 * pad)

        ax.set_title(subset)
        ax.set_xlabel("")
        ax.set_ylabel("Sample JBS Axis Score distances")
    fig.tight_layout()
    return fig


def centre_distance_panel(data: pd.DataFrame, r: float, p: float):
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.regplot(data=data, x="JBS axis score", y="Relative distance from centroid",
                line_kws={"color": "blue"}, scatter_kws={"alpha": 0.5, "color": "black"}, ax=ax)
    ax.text(0.6, 0.12, f"Pearson r = {r:.2f}\np-value = {p:.2e}", transform=ax.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=None, alpha=0.7))
    ax.set_xlabel("JBS Axis Score")
    ax.set_ylabel("Relative distance from centroid")
    return fig


def main() -> None:
    entry = cohort("bulk", DATASET)
    loadings = symbol_loadings(entry["axis_loadings"], config()["gene_symbols"])
    table = sample_table(entry, loadings)
    subsets = {"All samples": table, "CEL samples": table[table["CEL"].eq("Yes")]}

    p_all = patient_clustering(subsets["All samples"], MIN_SAMPLES, N_PERMUTATIONS, SEED)
    p_cel = patient_clustering(subsets["CEL samples"], MIN_SAMPLES, N_PERMUTATIONS, SEED)
    for name, p in {"all samples": p_all, "CEL samples": p_cel}.items():
        adjusted = multipletests(p, method="fdr_bh")[1]
        print(f"Per-patient clustering, {name}:")
        print(pd.DataFrame({"p": p, "p (BH)": adjusted}).round(4).to_string())

    distances = pd.concat([pair_distances(data).assign(Samples=name)
                           for name, data in subsets.items()], ignore_index=True)
    p_dispersion = {name: within_patient_p(data, distances[distances["Samples"] == name],
                                           N_PERMUTATIONS, SEED)
                    for name, data in subsets.items()}
    print(f"Within- vs between-patient distances: {p_dispersion}")

    located = table.dropna(subset=["Relative distance from centroid"])
    r, p = pearsonr(located["JBS axis score"], located["Relative distance from centroid"])
    print(f"Distance from centroid: n = {len(located)}, r = {r:.4f}, p = {p:.3e}")

    set_paper_theme()
    panels = {
        "multiregion_axis_score": (
            patient_panel(table, p_all, p_cel, contrast_enhancement_colours()),
            table.assign(**{"p": table["Patient"].map(p_all), "p CEL": table["Patient"].map(p_cel)}),
        ),
        "multiregion_axis_score_dispersion": (
            dispersion_panel(distances, p_dispersion),
            distances.drop(columns=["first", "second"]).assign(p=distances["Samples"].map(p_dispersion)),
        ),
        "multiregion_axis_score_vs_centre_distance": (
            centre_distance_panel(located, r, p),
            located[["Sample", "Patient", "JBS axis score", "Relative distance from centroid"]],
        ),
    }
    for name, (fig, plotted) in panels.items():
        save_panel(fig, plotted, name)
        plt.close(fig)


if __name__ == "__main__":
    main()
