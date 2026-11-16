"""TRS against the change in JBS axis score between a patient's two tumours.

Both are linear in log2 FPKM, so patients fall on y = x, except those of the Stead
centre whose FPKM matrix lacks ten axis genes. Included pairs of both cohorts
(`Include` and pair purity).

Example:
    pixi run trs_vs_axis_score_difference
"""

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

from common import (cohort, config, fpkm_axis_scores, resolve, save_panel, selected,
                    set_paper_theme, trs_table)

FIGSIZE = (3.2, 3.0)


def score_difference(name: str) -> pd.DataFrame:
    """TRS and Δ JBS axis score, one row per included patient of a bulk cohort."""
    entry = cohort("bulk", name)
    patients = trs_table(name)

    obs = ad.read_h5ad(resolve(entry["path"]), backed="r").obs
    samples = obs[selected(obs, entry["filters"])
                  & obs[entry["patient_column"]].isin(patients.index)]
    by_stage = pd.DataFrame({
        "patient": samples[entry["patient_column"]].astype(object),
        "stage": samples[entry["stage_column"]].astype(object),
        "score": fpkm_axis_scores(entry).reindex(samples.index),
    }).pivot(index="patient", columns="stage", values="score")

    return pd.DataFrame({
        "Cohort": patients["Cohort"],
        "TRS": patients["TRS"],
        "Δ JBS axis score": by_stage[entry["recurrent_value"]] - by_stage[entry["primary_value"]],
    }).rename_axis("Patient").reset_index()


def scatter(table: pd.DataFrame):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.scatter(table["TRS"], table["Δ JBS axis score"], color="black", alpha=0.3)
    ax.axline((0, 0), slope=1, linestyle="--", color="grey", linewidth=1)
    ax.set_xlabel("TRS")
    ax.set_ylabel("Δ JBS Axis Score")
    fig.tight_layout()
    return fig


def main() -> None:
    cohorts = [name for name, entry in config()["datasets"]["bulk"].items() if "log2fc" in entry]
    table = pd.concat([score_difference(name) for name in cohorts], ignore_index=True)

    offset = (table["Δ JBS axis score"] - table["TRS"]).abs()
    print(f"{len(table)} patients ({table['Cohort'].value_counts().to_dict()}); "
          f"{int((offset > 1e-9).sum())} off the line, by at most {offset.max():.3f}")

    set_paper_theme()
    fig = scatter(table)
    save_panel(fig, table, "trs_vs_axis_score_difference")
    plt.close(fig)


if __name__ == "__main__":
    main()
