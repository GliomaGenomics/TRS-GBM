"""TRS against the clinical profile of the primary tumour.

TRS and its magnitude against continuous variables (age, surgical interval,
survival, MGMT expression, purity; Pearson's r per cohort) and discrete ones
(gender, location, MGMT promoter methylation, subtype; boxplots per cohort).
Gender is recorded only for Discovery, MGMT methylation only for EORTC.

Example:
    pixi run trs_clinical_correlates
"""

from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pingouin as pg
import seaborn as sns
from scipy.stats import pearsonr
from statannotations.Annotator import Annotator

from common import (PALETTE, cohort, correlation_label, save_panel,
                    set_correlate_style, trs_table)

COHORTS = ["EORTC", "Discovery_Stead"]

# Locations outside these four are pooled; the rest are too sparse to test.
LOCATIONS = ["temporal", "frontal", "parietal", "occipital"]

RENAMES = {
    "MGMT": {"M": "Methylated", "U": "Unmethylated"},
    "Gender": {"M": "Male", "F": "Female"},
}


def clinical_table(names: list[str]) -> pd.DataFrame:
    table = pd.concat([trs_table(name) for name in names])
    for column, renames in RENAMES.items():
        table[column] = table[column].astype(object).replace(renames)
    table["Location.Primary"] = table["Location.Primary"].astype(object)
    table.loc[~table["Location.Primary"].isin(LOCATIONS), "Location.Primary"] = "other"
    table["abs(TRS)"] = table["TRS"].abs()
    return table


def continuous_panel(metadata: pd.DataFrame, variables: dict, score: str):
    n_cols = 5
    n_vars = len(variables)
    n_rows = int(np.ceil(n_vars / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3),
                             squeeze=False, sharey=True)
    rows = []
    for idx, (column, label) in enumerate(variables.items()):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]

        for name in metadata["Cohort"].unique():
            data = metadata[metadata["Cohort"] == name]
            mask = np.isfinite(data[column].values)
            y = data[score].to_numpy()[mask]
            x = data[column].values[mask]
            r, p = pearsonr(x, y)
            sns.regplot(x=x, y=y, ax=ax, color=PALETTE[name],
                        label=correlation_label(name, r, p))
            rows.append(pd.DataFrame({
                "Patient": data.index[mask], "Cohort": name, "Variable": label,
                "Value": x, score: y, "R": r, "p": p,
            }))

        ax.set_ylabel(score if idx == 0 else "")
        ax.set_xlabel(label)
        ax.legend(loc="best", frameon=True, fontsize=10,
                  handletextpad=0.4, handlelength=0.3, labelspacing=0.3)

    for idx in range(n_vars, n_rows * n_cols):
        fig.delaxes(axes[idx // n_cols][idx % n_cols])

    sns.despine()
    fig.tight_layout()
    return fig, pd.concat(rows, ignore_index=True)


def annotate_anova(ax, data, value, hue, y_pad=10, bar_height=0.05):
    """Mark each cohort with a non-significant Welch ANOVA; report post-hoc tests otherwise."""
    for i, name in enumerate(data["Cohort"].unique()):
        subset = data[data["Cohort"] == name]
        if subset[hue].dropna().nunique() < 2:
            continue

        p = pg.welch_anova(data=subset, dv=value, between=hue)["p_unc"].iloc[0]
        print(f"  {name} {hue}: Welch ANOVA p={p:.3f}")

        if p >= 0.05:
            start, end = i - 0.4, i + 0.4
            y = subset[value].max() + y_pad
            ax.plot([start, end], [y, y], lw=1.5, c="k")
            ax.plot([start, start], [y, y - bar_height], lw=1, c="k")
            ax.plot([end, end], [y, y - bar_height], lw=1, c="k")
            ax.text((start + end) / 2, y + y_pad * 0.2, "ns",
                    ha="center", va="bottom", fontsize=10, color="k")
        else:
            print(pg.pairwise_tests(data=subset, dv=value, between=hue, padjust="bonf"))


def annotate_pairs(ax, data, value, hue, levels):
    pairs = [
        ((name, levels[0]), (name, levels[1]))
        for name in data["Cohort"].unique()
        if set(levels).issubset(set(data[data["Cohort"] == name][hue]))
    ]
    annotator = Annotator(ax, pairs, data=data, x="Cohort", y=value, hue=hue)
    annotator.configure(test="t-test_ind", text_format="star", loc="inside", verbose=0)
    annotator.apply_and_annotate()


def annotate_subtype(ax, data, value):
    """Welch ANOVA per cohort; where significant, Classical is tested against the rest."""
    for name in data["Cohort"].unique():
        subset = data[data["Cohort"] == name]
        groups = subset["Subtype"].dropna().unique()
        if len(groups) < 2:
            continue

        p = pg.welch_anova(data=subset, dv=value, between="Subtype")["p_unc"].iloc[0]
        print(f"  {name} Subtype: Welch ANOVA p={p:.3f}")
        if p >= 0.05:
            annotate_anova(ax, data, value, "Subtype")
            continue

        pairs = [((name, a), (name, b)) for a, b in combinations(groups, 2)
                 if "Classical" in (a, b)]
        annotator = Annotator(ax, pairs=pairs, data=data, x="Cohort", y=value, hue="Subtype")
        annotator.configure(test="t-test_welch", text_format="star", loc="inside",
                            comparisons_correction="bonferroni", verbose=1,
                            line_offset=100, line_height=0.0, text_offset=-3.3)
        annotator.apply_and_annotate()


def discrete_panel(metadata: pd.DataFrame, variables: dict, score: str, height: float,
                   subtype_annotator):
    fig, axes = plt.subplots(1, 4, figsize=(15, height), sharey=True)
    rows = []

    for ax, (column, title) in zip(axes, variables.items()):
        if column == "Gender":
            data = metadata[metadata["Gender"] != "Unknown"]
            sns.boxplot(data=data, x="Cohort", y=score, hue=column, ax=ax)
            annotate_pairs(ax, data, score, column, ["Male", "Female"])
        elif column == "Location.Primary":
            data = metadata
            sns.boxplot(data=data, x="Cohort", y=score, hue=column, palette="husl", ax=ax)
            annotate_anova(ax, data, score, column)
        elif column == "MGMT":
            data = metadata[metadata["MGMT"].isin(["Methylated", "Unmethylated"])]
            sns.boxplot(data=data, x="Cohort", y=score, hue=column, ax=ax)
            annotate_pairs(ax, data, score, column, ["Methylated", "Unmethylated"])
        else:
            data = metadata
            sns.boxplot(data=data, x="Cohort", y=score, hue=column, palette=PALETTE, ax=ax)
            subtype_annotator(ax, data, score)

        ax.set_title(title)
        ax.set_xlabel("")
        sns.despine(ax=ax)
        rows.append(pd.DataFrame({
            "Patient": data.index, "Cohort": data["Cohort"], "Variable": title,
            "Group": data[column].astype(object), score: data[score],
        }))

    for i, ax in enumerate(axes):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc="center left", bbox_to_anchor=(0.95, 0.5),
                  frameon=False, handletextpad=0.5, handlelength=1.5, labelspacing=0.3)
        if i > 0:
            ax.set_ylabel("")

    fig.tight_layout()
    return fig, pd.concat(rows, ignore_index=True)


def main() -> None:
    set_correlate_style()
    metadata = clinical_table(COHORTS)
    entry = cohort("bulk", COHORTS[0])
    continuous, discrete = entry["continuous_clinical"], entry["discrete_clinical"]

    for score, prefix in (("TRS", "trs"), ("abs(TRS)", "abs_trs")):
        fig, table = continuous_panel(metadata, continuous, score)
        save_panel(fig, table, f"{prefix}_vs_clinical_continuous")
        plt.close(fig)

    # The magnitude panel is a little taller than the signed one, and
    # is the only one that tests subtypes pairwise where the ANOVA is significant.
    panels = (
        ("TRS", "trs", 3.0, lambda ax, data, score: annotate_anova(ax, data, score, "Subtype")),
        ("abs(TRS)", "abs_trs", 3.2, annotate_subtype),
    )
    for score, prefix, height, subtype_annotator in panels:
        print(f"{score} by categorical clinical variable:")
        fig, table = discrete_panel(metadata, discrete, score, height, subtype_annotator)
        save_panel(fig, table, f"{prefix}_vs_clinical_discrete")
        plt.close(fig)


if __name__ == "__main__":
    main()
