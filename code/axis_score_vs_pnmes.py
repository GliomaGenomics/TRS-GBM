"""Proneural–mesenchymal (PN:MES) score against JBS axis score in bulk tumours.

Both tumours of the included pairs of EORTC and Discovery (`Include` and pair
purity). PN:MES score: raw counts scaled to 10,000 per sample, log10(1 + x), dotted
with the PC_1 loadings distributed with Wang (2022). Per cohort, a least-squares
line and Pearson's r. Drawn at 55 × 55 mm.

Example:
    pixi run axis_score_vs_pnmes
"""

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import yaml
from scipy.stats import pearsonr

from common import (cohort, config, fpkm_axis_scores, pnmes_scores, resolve, save_panel,
                    selected, set_nature_style, trs_table)

COHORTS = ["EORTC", "Discovery"]
# A third of a double-column (183 mm) figure.
SIZE_MM = (55, 55)
# Marker area in pt².
MARKER_SIZE = 12


def cohort_colours() -> dict:
    with open(resolve("config/palette.yaml")) as handle:
        return yaml.safe_load(handle)["cohort"]


def cohort_scores(name: str, loadings: pd.Series) -> pd.DataFrame:
    """Both scores for every sample of an included pair in a bulk cohort."""
    entry = cohort("bulk", name)
    patients = trs_table(name).index

    data = ad.read_h5ad(resolve(entry["path"]), backed="r")
    chosen = selected(data.obs, entry["filters"]) & data.obs[entry["patient_column"]].isin(patients)
    samples = data[chosen.to_numpy()].to_memory()

    return pd.DataFrame({
        "Sample": samples.obs_names,
        "Patient": samples.obs[entry["patient_column"]].astype(object),
        "Stage": samples.obs[entry["stage_column"]].astype(object),
        "Cohort": entry["label"],
        "JBS axis score": fpkm_axis_scores(entry).reindex(samples.obs_names),
        "PN:MES score": pnmes_scores(samples.to_df(), loadings),
    }).reset_index(drop=True)


def correlations(table: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlation of the two scores within each cohort."""
    rows = []
    for name, group in table.groupby("Cohort", sort=False):
        r, p = pearsonr(group["PN:MES score"], group["JBS axis score"])
        rows.append({"Cohort": name, "n": len(group), "r": r, "p": p})
    return pd.DataFrame(rows)


def scatter(table: pd.DataFrame, stats: pd.DataFrame, colours: dict):
    fig, ax = plt.subplots(figsize=[mm / 25.4 for mm in SIZE_MM], layout="constrained")
    sns.scatterplot(data=table, x="JBS axis score", y="PN:MES score", hue="Cohort",
                    palette=colours, s=MARKER_SIZE, alpha=0.7, ax=ax, legend=False)

    for name, group in table.groupby("Cohort", sort=False):
        sns.regplot(data=group, x="JBS axis score", y="PN:MES score", scatter=False,
                    ci=None, line_kws={"color": colours[name], "linewidth": 1}, ax=ax)

    for i, row in stats.iterrows():
        ax.text(0.02, 0.17 - i * 0.09, f"{row['Cohort']}: r={row['r']:.2f}, p={row['p']:.1e}",
                transform=ax.transAxes, va="top", color="black",
                bbox=dict(boxstyle="round,pad=0.4", fc=colours[row["Cohort"]], ec="none", alpha=0.1))

    ax.set_xlabel("JBS Axis Score")
    ax.set_ylabel("PN:MES")
    return fig


def main() -> None:
    settings = config()["pn_mes_loadings"]
    loadings = pd.read_csv(resolve(settings["path"]), index_col=0)[settings["column"]]

    table = pd.concat([cohort_scores(name, loadings) for name in COHORTS], ignore_index=True)
    stats = correlations(table)
    print(stats.to_string(index=False))

    set_nature_style()
    fig = scatter(table, stats, cohort_colours())
    save_panel(fig, table, "axis_score_vs_pnmes", bbox_inches=None)
    plt.close(fig)


if __name__ == "__main__":
    main()
