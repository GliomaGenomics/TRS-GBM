"""The patients of one single-cell dataset, for the patient metadata table.

One row per patient: published sample metadata, TRS (recurrent minus primary
whole-tumour pseudobulk score) and, per stage (primary_, recurrent_,
second_recurrence_), the sample, cell counts, purity, pseudobulk axis scores and
cell proportions. The metadata comes from the cells' .obs, so run it on the real
cells (hpc/run_single_cell.sh).

Writes figures/tables/patient_table_<dataset>.csv and _columns.csv.

Example:
    pixi run single_cell_patient_table nomura_2025
"""

import argparse
import re

import anndata as ad
import numpy as np
import pandas as pd

from common import FIGURES, cohort, natural_order, pseudobulk_table, resolve

TABLES = FIGURES / "tables"

# Stages outside the primary-recurrent pair, and the prefix their columns take.
OTHER_STAGE_PREFIXES = {"2nd recurrence": "second_recurrence"}

# The pseudobulk each axis score column is taken from.
SCORES = {
    "whole_tumour_score": ("All", "JBS axis score of the pseudobulk of all its cells"),
    "cancer_cell_score": ("Malignant", "JBS axis score of the pseudobulk of its malignant cells"),
    "microenvironment_score": ("Non-malignant",
                               "JBS axis score of the pseudobulk of its non-malignant cells"),
}


def snake(name: str) -> str:
    return re.sub(r"[^0-9a-z]+", "_", name.lower()).strip("_")


def sample_metadata(entry: dict) -> pd.DataFrame:
    """The published metadata of every sample, one row per sample."""
    obs = ad.read_h5ad(resolve(entry["path"]), backed="r").obs
    samples = obs.drop_duplicates(entry["sample_column"]).set_index(entry["sample_column"])
    return samples[list(entry["sample_metadata_columns"])].rename(columns=snake)


def sample_values(obs: pd.DataFrame, entry: dict) -> tuple[pd.DataFrame, dict]:
    """Every sample's cell counts, purity, scores and proportions, one row per
    sample, and a description of each."""
    whole = obs[obs["Population"] == "All"].set_index("Sample")
    values = pd.DataFrame({
        "patient": whole["Patient"].astype(str),
        "stage": whole["PvsR"].astype(str),
        "cells": whole["total_cells"].astype("Int64"),
        "malignant_cells": whole["n_Malignant"].astype("Int64"),
        "purity": whole["Purity"],
    })
    described = {
        "cells": "Cells",
        "malignant_cells": "Malignant cells",
        "purity": "Fraction of cells that are malignant",
    }
    for column, (population, description) in SCORES.items():
        scores = obs[obs["Population"] == population].set_index("Sample")["Axis Score"]
        values[column] = scores.reindex(values.index)
        described[column] = description

    # The malignant population and cell type are the purity column again.
    for population in dict.fromkeys(entry["cell_populations"].values()):
        if population != "Malignant" and f"comp_pop_{population}" in whole.columns:
            column = f"population_{snake(population)}"
            values[column] = whole[f"comp_pop_{population}"]
            described[column] = f"Fraction of cells in the {population} population"
    for cell_type in entry["cell_populations"]:
        if cell_type != "Malignant" and f"comp_{cell_type}" in whole.columns:
            column = f"cell_type_{snake(cell_type)}"
            values[column] = whole[f"comp_{cell_type}"]
            described[column] = f"Fraction of cells that are {cell_type}"
    for state in entry["malignant_state_order"]:
        column = f"malignant_state_{snake(state)}"
        values[column] = whole[f"malignant_comp_{state}"]
        described[column] = f"Fraction of malignant cells in the {state} state"
    return values, described


def patient_table(values: pd.DataFrame, described: dict, metadata: pd.DataFrame,
                  entry: dict) -> tuple[pd.DataFrame, dict]:
    """One row per patient, and a description of every column."""
    published = {snake(c): text for c, text in entry["sample_metadata_columns"].items()}
    samples = values.join(metadata)
    by_patient = samples.groupby("patient")
    shared = [c for c in metadata.columns if (by_patient[c].nunique(dropna=False) <= 1).all()]
    per_stage = ([c for c in metadata.columns if c not in shared]
                 + [c for c in values.columns if c not in ("patient", "stage")])

    table = by_patient[shared].first()
    whole = samples.pivot(index="patient", columns="stage", values="whole_tumour_score")
    trs = whole[entry["recurrent_value"]] - whole[entry["primary_value"]]
    table["trs"] = trs
    table["trs_sign"] = np.where(trs.isna(), None, np.where(trs > 0, "+TRS", "−TRS"))
    descriptions = {
        "patient": f"Patient in {entry['label']}",
        **{c: published[c] for c in shared},
        "trs": "TRS: recurrent minus primary whole-tumour score (blank unless the patient "
               "has both samples)",
        "trs_sign": "+TRS or −TRS, by the sign of TRS",
    }

    stage_prefixes = {entry["primary_value"]: "primary", entry["recurrent_value"]: "recurrent",
                      **OTHER_STAGE_PREFIXES}
    for stage, prefix in stage_prefixes.items():
        at_stage = samples[samples["stage"] == stage]
        if at_stage.empty:
            continue
        at_stage = at_stage.rename_axis("sample").reset_index().set_index("patient")
        for column in ["sample", *per_stage]:
            name = f"{prefix}_{column}"
            table[name] = at_stage[column].reindex(table.index)
            meaning = "Sample" if column == "sample" else {**described, **published}[column]
            descriptions[name] = f"{meaning}: {stage.lower()} sample"
    # Metadata recorded at one stage only, such as the time to recurrence, leaves
    # the other stage's column empty.
    table = table.loc[natural_order(table.index)].dropna(axis=1, how="all")
    return table.reset_index(), descriptions


def main(name: str) -> None:
    entry = cohort("single_cell", name)
    values, described = sample_values(pseudobulk_table(name), entry)
    table, descriptions = patient_table(values, described, sample_metadata(entry), entry)

    TABLES.mkdir(parents=True, exist_ok=True)
    table.to_csv(TABLES / f"patient_table_{name}.csv", index=False)
    pd.DataFrame({
        "Column": table.columns,
        "Description": [descriptions[c] for c in table.columns],
    }).to_csv(TABLES / f"patient_table_{name}_columns.csv", index=False)
    print(f"wrote patient_table_{name}: {len(table)} patients, {table.shape[1]} columns, "
          f"TRS for {int(table['trs'].notna().sum())}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dataset", help="single-cell dataset in config/datasets.yaml")
    main(parser.parse_args().dataset)
