"""JBS axis score regressed on the cell composition of a sample.

Ordinary least squares on centred log-ratio transformed cell fractions. Panels:
observed against predicted score, and each component's coefficient with its 95%
confidence interval.

    --component celltype   cell types, malignant cells pooled
    --component state      malignant cells split into their states

Example:
    pixi run axis_score_composition_regression nomura_2025 --component celltype
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from skbio.stats.composition import clr, multi_replace
from sklearn.metrics import r2_score

from common import cell_colours, cohort, pseudobulk_table, save_panel, set_paper_theme

FIGSIZE = (10, 4)


def composition_matrix(name: str, component: str) -> pd.DataFrame:
    """Every sample's composition, one column per component, as proportions of its cells."""
    entry = cohort("single_cell", name)
    obs = pseudobulk_table(name)
    obs = obs[obs["Population"] == "All"]

    cell_types = [t for t in entry["cell_populations"] if f"comp_{t}" in obs.columns]
    if component == "celltype":
        columns = [f"comp_{t}" for t in cell_types]
    else:
        states = [c for c in obs.columns if c.startswith("comp_malignant_")]
        columns = [f"comp_{t}" for t in cell_types if t != "Malignant"] + states

    matrix = obs[columns]
    matrix.columns = [c.removeprefix("comp_") for c in columns]
    return matrix.set_index(obs["Sample"]).assign(**{"Axis Score": obs["Axis Score"].to_numpy()})


def fit_composition(matrix: pd.DataFrame):
    """CLR transform the compositions, then regress the axis score on them."""
    parts = matrix.drop(columns="Axis Score").to_numpy()
    if not np.allclose(parts.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("Compositions do not sum to one; the components do not partition the cells.")

    # Zero proportions have no logarithm, so they are replaced by a small value
    # before the centred log-ratio transform.
    scores = matrix["Axis Score"].to_numpy()
    model = sm.OLS(scores, sm.add_constant(clr(multi_replace(parts)))).fit()

    intervals = np.asarray(model.conf_int(alpha=0.05))
    coefficients = pd.DataFrame({
        "Component": ["Intercept"] + list(matrix.columns.drop("Axis Score")),
        "β": np.asarray(model.params),
        "ci_lower": intervals[:, 0],
        "ci_upper": intervals[:, 1],
    }).query("Component != 'Intercept'")

    return model.fittedvalues, coefficients


def regression_panel(matrix: pd.DataFrame, predicted: np.ndarray,
                     coefficients: pd.DataFrame, title: str):
    observed = matrix["Axis Score"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, gridspec_kw={"width_ratios": [1.2, 1]})

    sns.regplot(x=observed, y=predicted, ax=axes[0],
                scatter_kws={"alpha": 0.7, "color": "slategrey"},
                line_kws={"color": "black"})
    axes[0].set_xlabel("JBS Axis Score (observed)")
    axes[0].set_ylabel("JBS Axis Score (predicted)")
    axes[0].axline((0, 0), slope=1, linestyle="--", color="black")
    axes[0].text(0.05, 0.95, f"$R^2 = {r2_score(observed, predicted):.3f}$",
                 transform=axes[0].transAxes, ha="left", va="top", fontsize=10,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.5))

    ranked = coefficients.sort_values("β", ascending=False).reset_index(drop=True)
    positions = np.arange(len(ranked))
    axes[1].errorbar(
        x=ranked["β"].to_numpy(), y=positions,
        xerr=[ranked["β"] - ranked["ci_lower"], ranked["ci_upper"] - ranked["β"]],
        fmt="none", ecolor="black", alpha=0.7, capsize=3,
    )
    colours = pd.Series([c.split("_")[-1] for c in ranked["Component"]]).map(cell_colours())
    axes[1].scatter(ranked["β"].to_numpy(), positions, color=colours, s=40,
                    zorder=3, edgecolor="black")
    axes[1].axvline(0, color="lightgrey", linewidth=1.25, linestyle="--", zorder=0)
    axes[1].set_yticks(positions)
    axes[1].set_yticklabels(ranked["Component"])
    axes[1].set_xlabel("Coefficient (β)")

    fig.suptitle(title)
    fig.tight_layout()
    return fig


def source_table(matrix: pd.DataFrame, predicted: np.ndarray,
                 coefficients: pd.DataFrame) -> pd.DataFrame:
    """The values behind both panels: one row per sample, then one per coefficient."""
    fit = pd.DataFrame({
        "Quantity": "Axis score",
        "Key": matrix.index,
        "Observed": matrix["Axis Score"].to_numpy(),
        "Predicted": predicted,
    })
    weights = pd.DataFrame({
        "Quantity": "Coefficient",
        "Key": coefficients["Component"],
        "β": coefficients["β"],
        "ci_lower": coefficients["ci_lower"],
        "ci_upper": coefficients["ci_upper"],
    })
    return pd.concat([fit, weights], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "dataset", help="single-cell dataset in config/datasets.yaml, e.g. nomura_2025"
    )
    parser.add_argument("--component", required=True, choices=("celltype", "state"))
    arguments = parser.parse_args()

    set_paper_theme()
    matrix = composition_matrix(arguments.dataset, arguments.component)
    predicted, coefficients = fit_composition(matrix)

    # The dataset label brackets its year; this title supplies its own brackets.
    label = cohort("single_cell", arguments.dataset)["label"].replace("(", "").replace(")", "")
    fig = regression_panel(matrix, predicted, coefficients,
                           f"Tumour Composition Regression of Axis Score ({label})")
    save_panel(fig, source_table(matrix, predicted, coefficients),
               f"axis_score_composition_regression_{arguments.dataset}_{arguments.component}")
    plt.close(fig)


if __name__ == "__main__":
    main()
