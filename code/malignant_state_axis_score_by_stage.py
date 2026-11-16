"""Per-cell JBS axis scores of each malignant state, primary against recurrent.

Malignant cells of patients with a primary and a first recurrence, scored with the
loadings in the single-cell .var. Per state, a Wilcoxon signed-rank test of the
patients' paired median scores (at least MIN_CELLS cells at both stages),
Benjamini-Hochberg across states. Repeated for +TRS and -TRS patients (sign of the
whole-tumour pseudobulk TRS), and each state's shift correlated with TRS
(Spearman, Benjamini-Hochberg).

Writes malignant_state_axis_score_by_stage_<dataset>, its _by_trs_sign_ version and
malignant_state_axis_score_shift_vs_trs_<dataset>, each with a _tests table.

Example:
    pixi run malignant_state_axis_score_by_stage nomura_2025
"""

import argparse

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

from common import (PALETTE, SOURCE_DATA, cell_axis_scores, cell_colours, cohort,
                    pseudobulk_axis_scores, resolve, save_panel, set_paper_theme)

PANEL_WIDTH = 4
ROW_HEIGHT = 0.5

STAGES = ("Primary", "Recurrent")
STAGE_COLOURS = {"Primary": "#0197f6", "Recurrent": PALETTE["Recurrent"]}
TRS_SIGNS = ("+TRS", "−TRS")

# Cells of a state a patient needs at both stages for their medians to enter that
# state's test; a median over fewer cells is too noisy to pair.
MIN_CELLS = 10


def scored_cells(data: ad.AnnData, entry: dict) -> pd.DataFrame:
    """Every cell of a dataset, with its patient, stage, malignant state and axis score.

    The state is blank for every cell outside the malignant population.
    """
    return pd.DataFrame({
        "Malignant_state": data.obs[entry["malignant_state_column"]],
        "Patient": data.obs[entry["patient_column"]].astype(str),
        "Stage": data.obs[entry["stage_column"]].astype(object),
        "Axis Score": cell_axis_scores(data, entry["pseudocount"]),
    })


def mean_score_order(cells: pd.DataFrame) -> np.ndarray:
    """States ordered by their mean score over every cell, the unpaired ones included."""
    means = cells.groupby("Malignant_state", observed=False)["Axis Score"].mean()
    return means.sort_values().index.values


def paired_cells(cells: pd.DataFrame, entry: dict) -> pd.DataFrame:
    """The cells of the primary and recurrent samples of patients that have both.

    A stage outside the pair, such as a second recurrence, is not part of the
    comparison, and a patient missing one of the two cannot show a shift. A
    patient qualifies on any of their cells, not only their malignant ones.
    """
    stages = dict(zip((entry["primary_value"], entry["recurrent_value"]), STAGES))
    paired = cells.assign(Stage=cells["Stage"].map(stages))
    paired = paired.dropna(subset=["Axis Score", "Stage", "Patient"])

    present = paired.groupby("Patient", observed=True)["Stage"].nunique()
    return paired[paired["Patient"].isin(present.index[present == len(STAGES)])]


def whole_tumour_trs(data: ad.AnnData, entry: dict) -> pd.Series:
    """Each paired patient's TRS: recurrent minus primary whole-tumour pseudobulk score.

    A sample's whole-tumour pseudobulk is the gene-wise sum of the raw counts of
    all of its cells, whatever their type.
    """
    samples = data.obs[entry["sample_column"]].astype(str)
    names = samples.unique()
    summed = np.vstack([np.asarray(data.X[(samples == sample).to_numpy()].sum(axis=0)).ravel()
                        for sample in names])
    scores = pseudobulk_axis_scores(ad.AnnData(summed, var=data.var), entry["pseudocount"])

    first_cells = data.obs.assign(Sample=samples).drop_duplicates("Sample").set_index("Sample")
    stages = dict(zip((entry["primary_value"], entry["recurrent_value"]), STAGES))
    table = pd.DataFrame({
        "Patient": first_cells[entry["patient_column"]].astype(str),
        "Stage": first_cells[entry["stage_column"]].astype(object).map(stages),
        "Axis Score": scores.to_numpy(),
    }, index=names).dropna(subset=["Stage"])
    by_stage = table.pivot(index="Patient", columns="Stage", values="Axis Score")
    return (by_stage["Recurrent"] - by_stage["Primary"]).dropna()


def state_shifts(cells: pd.DataFrame) -> pd.Series:
    """Each patient's shift in each state: recurrent minus primary median cell score.

    Only patients with MIN_CELLS of the state's cells at both stages have one.
    """
    per_patient = (cells.groupby(["Malignant_state", "Patient", "Stage"], observed=True)
                   ["Axis Score"].agg(["median", "size"]).unstack("Stage"))
    enough = (per_patient["size"] >= MIN_CELLS).all(axis=1)
    shifts = per_patient["median", "Recurrent"] - per_patient["median", "Primary"]
    return shifts[enough].rename("Shift")


def stage_tests(cells: pd.DataFrame, order: np.ndarray) -> pd.DataFrame:
    """Each state's shift at recurrence, tested with the patient as the unit.

    A state with no patient that has a shift is not tested.
    """
    by_state = state_shifts(cells).groupby(level="Malignant_state", observed=True)
    tests = pd.DataFrame({
        "Patients": by_state.size(),
        "Median shift": by_state.median(),
        "p": by_state.apply(lambda shift: wilcoxon(shift).pvalue),
    }).reindex(order)
    tests["Patients"] = tests["Patients"].fillna(0).astype(int)

    tested = tests["p"].notna()
    tests.loc[tested, "q"] = multipletests(tests.loc[tested, "p"], method="fdr_bh")[1]
    return tests.rename_axis("Malignant state").reset_index()


def trs_correlations(shifts: pd.DataFrame, order: np.ndarray) -> pd.DataFrame:
    """Spearman correlation of each state's shift with TRS, the patient as the unit.

    A state with fewer than three patients has no p-value and is not tested.
    """
    rows = []
    for state in order:
        patients = shifts[shifts["Malignant_state"] == state]
        result = spearmanr(patients["TRS"], patients["Shift"])
        rows.append({"Malignant state": state, "Patients": len(patients),
                     "rho": result.statistic, "p": result.pvalue})
    tests = pd.DataFrame(rows)

    tested = tests["p"].notna()
    tests.loc[tested, "q"] = multipletests(tests.loc[tested, "p"], method="fdr_bh")[1]
    return tests


def draw_shifts(ax, cells: pd.DataFrame, order: np.ndarray, tests: pd.DataFrame,
                title: str) -> None:
    sns.boxplot(
        data=cells,
        x="Axis Score",
        y="Malignant_state",
        hue="Stage",
        palette=STAGE_COLOURS,
        order=order,
        hue_order=list(STAGES),
        fliersize=1,
        ax=ax,
    )

    ax.set_ylabel("Malignant state")
    ax.set_xlabel("JBS Axis Score")
    ax.legend(title="")
    ax.set_title(title)

    for row, (q, n) in enumerate(zip(tests["q"], tests["Patients"])):
        label = f"q = {q:.2g} (n = {n})" if n else "n = 0"
        ax.text(1.02, row, label, transform=ax.get_yaxis_transform(), va="center",
                fontsize="small")


def shift_panel(cells: pd.DataFrame, order: np.ndarray, tests: pd.DataFrame, title: str):
    fig, ax = plt.subplots(figsize=(PANEL_WIDTH, ROW_HEIGHT * len(order)))
    draw_shifts(ax, cells, order, tests, title)
    return fig


def trs_sign_panel(cells: pd.DataFrame, order: np.ndarray, tests: pd.DataFrame, label: str):
    fig, axes = plt.subplots(1, 2, figsize=(2 * PANEL_WIDTH + 2, ROW_HEIGHT * len(order)),
                             sharex=True, sharey=True, gridspec_kw={"wspace": 0.5})
    for ax, sign in zip(axes, TRS_SIGNS):
        draw_shifts(ax, cells[cells["TRS sign"] == sign], order,
                    tests[tests["TRS sign"] == sign], f"{label}, {sign}")
    axes[1].set_ylabel("")
    axes[1].get_legend().remove()
    return fig


def shift_vs_trs_panel(shifts: pd.DataFrame, order: np.ndarray, tests: pd.DataFrame,
                       label: str):
    n_columns = 5
    n_rows = -(-len(order) // n_columns)
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(2 * n_columns, 2.2 * n_rows),
                             sharex=True, sharey=True)
    colours = cell_colours()
    for ax, state, n, rho, q in zip(axes.flat, tests["Malignant state"], tests["Patients"],
                                    tests["rho"], tests["q"]):
        patients = shifts[shifts["Malignant_state"] == state]
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvline(0, color="grey", linewidth=0.5)
        ax.scatter(patients["TRS"], patients["Shift"], s=12, color=colours[state],
                   edgecolor="black", linewidth=0.3)
        stats = f"ρ = {rho:.2f}, q = {q:.2g} (n = {n})" if pd.notna(q) else f"n = {n}"
        ax.set_title(f"{state}\n{stats}", fontsize="small")
    for ax in axes.flat[len(order):]:
        ax.set_visible(False)

    fig.supxlabel("TRS")
    fig.supylabel("Change in median JBS Axis Score")
    fig.suptitle(label)
    fig.tight_layout()
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "dataset", help="single-cell dataset in config/datasets.yaml, e.g. nomura_2025"
    )
    name = parser.parse_args().dataset
    entry = cohort("single_cell", name)
    data = ad.read_h5ad(resolve(entry["path"]))

    set_paper_theme()
    cells = scored_cells(data, entry)
    order = mean_score_order(cells)

    paired = paired_cells(cells, entry)
    drawn = paired[paired["Malignant_state"].notna()]
    print(f"drawing {len(drawn)} malignant cells from "
          f"{drawn['Patient'].nunique()} paired patients")

    tests = stage_tests(drawn, order)
    fig = shift_panel(drawn, order, tests, entry["label"])
    panel = f"malignant_state_axis_score_by_stage_{name}"
    save_panel(fig, drawn.rename_axis("Cell").reset_index(), panel)
    tests.to_csv(SOURCE_DATA / f"{panel}_tests.csv", index=False)
    plt.close(fig)

    trs = whole_tumour_trs(data, entry)
    signs = pd.Series(np.where(trs > 0, TRS_SIGNS[0], TRS_SIGNS[1]), index=trs.index)
    split = drawn.assign(TRS=drawn["Patient"].map(trs), **{"TRS sign": drawn["Patient"].map(signs)})
    print("paired patients by TRS sign:",
          split.groupby("TRS sign")["Patient"].nunique().to_dict())

    split_tests = pd.concat(
        [stage_tests(split[split["TRS sign"] == sign], order).assign(**{"TRS sign": sign})
         for sign in TRS_SIGNS])
    fig = trs_sign_panel(split, order, split_tests, entry["label"])
    panel = f"malignant_state_axis_score_by_stage_by_trs_sign_{name}"
    save_panel(fig, split.rename_axis("Cell").reset_index(), panel)
    split_tests.to_csv(SOURCE_DATA / f"{panel}_tests.csv", index=False)
    plt.close(fig)

    shifts = state_shifts(drawn).reset_index()
    shifts["TRS"] = shifts["Patient"].map(trs)
    shifts = shifts.dropna(subset=["TRS"])
    correlations = trs_correlations(shifts, order)
    fig = shift_vs_trs_panel(shifts, order, correlations, entry["label"])
    panel = f"malignant_state_axis_score_shift_vs_trs_{name}"
    save_panel(fig, shifts, panel)
    correlations.to_csv(SOURCE_DATA / f"{panel}_tests.csv", index=False)
    plt.close(fig)


if __name__ == "__main__":
    main()
