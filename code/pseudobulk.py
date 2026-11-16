"""Sample-level pseudobulks from single-cell counts.

Raw counts summed over a sample's cells: every cell (<sample>_All), the
non-malignant cells, and each cell population. Each carries its number of cells of
each cell type and malignant state. A second file holds pseudobulks of three
interleaved thirds of each sample's cells, as replicates.

Example:
    pixi run pseudobulk nomura_2025
"""

import argparse

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from common import cohort, resolve

DTYPE = "float32"

# Cells are dealt out into this many interleaved subsets, each pseudobulked on its
# own, so every sample gets that many replicates.
SUBSETS = 3


def annotate_cell_populations(data: ad.AnnData, entry: dict) -> ad.AnnData:
    """Label every cell with its population, in place."""
    cell_types = data.obs[entry["cell_type_column"]]
    if "Unknown" not in cell_types.cat.categories:
        cell_types = cell_types.cat.add_categories(["Unknown"])
    cell_types = cell_types.fillna("Unknown")

    populations = cell_types.map(entry["cell_populations"])
    unmapped = sorted(set(cell_types[populations.isna()].unique()))
    if unmapped:
        raise ValueError(
            f"Cell types with no population in the config: {unmapped}.\n"
            "Add them to cell_populations before pseudobulking."
        )

    data.obs[entry["cell_type_column"]] = cell_types
    data.obs["CellPopulation"] = populations
    return data


def check_constant_within_samples(obs: pd.DataFrame, entry: dict) -> None:
    for column in (entry["patient_column"], entry["stage_column"]):
        distinct = obs.groupby(entry["sample_column"], observed=False)[column].nunique()
        if (distinct > 1).any():
            raise ValueError(
                f"Column {column!r} is not constant within samples "
                f"{list(distinct.index[distinct > 1])}.\n"
                "Fix this before pseudobulking."
            )


def summed_counts(cells: ad.AnnData) -> np.ndarray:
    total = cells.X.sum(axis=0)
    if sparse.issparse(total):
        total = total.toarray()
    return np.asarray(total).ravel().astype(DTYPE, copy=False)


def cell_counts(obs: pd.DataFrame, entry: dict) -> dict[str, int]:
    """How many cells of each cell type, and of each malignant state, a pseudobulk holds."""
    counts = {
        f"n_{cell_type}": int(n)
        for cell_type, n in obs[entry["cell_type_column"]].value_counts().items()
    }
    counts.update({
        f"n_malignant_{state}": int(n)
        for state, n in obs[entry["malignant_state_column"]].value_counts().items()
    })
    return counts


def cell_groups(obs: pd.DataFrame, populations: list[str]) -> list[tuple[str, np.ndarray]]:
    """The sets of cells of one sample that each become a pseudobulk, in output order.

    The populations come from the whole dataset rather than from this sample, so a
    sample missing one still gets an empty pseudobulk for it.
    """
    population = obs["CellPopulation"]
    return [
        ("All", np.ones(len(obs), dtype=bool)),
        ("Non-malignant", (population != "Malignant").to_numpy()),
        *((name, (population == name).to_numpy()) for name in populations),
    ]


def pseudobulk_by_sample(data: ad.AnnData, entry: dict) -> ad.AnnData:
    check_constant_within_samples(data.obs, entry)
    samples = data.obs[entry["sample_column"]].unique().tolist()
    populations = data.obs["CellPopulation"].unique().tolist()

    matrix, records = [], []
    for sample in samples:
        cells = data[(data.obs[entry["sample_column"]] == sample).to_numpy()]
        for population, group in cell_groups(cells.obs, populations):
            pseudobulk = cells[group]
            matrix.append(summed_counts(pseudobulk))
            records.append({
                "Pseudobulk": f"{sample}_{population}",
                "Sample": sample,
                "Population": population,
                "Patient": cells.obs[entry["patient_column"]].iloc[0],
                "PvsR": cells.obs[entry["stage_column"]].iloc[0],
                "total_cells": int(group.sum()),
                **cell_counts(pseudobulk.obs, entry),
            })

    obs = pd.DataFrame(records).set_index("Pseudobulk")
    # A cell type absent from a whole subset leaves its count column missing, not zero.
    counted = obs.select_dtypes(include=[np.number]).columns
    obs[counted] = obs[counted].fillna(0)

    return ad.AnnData(
        X=np.vstack(matrix), obs=obs, var=data.var.copy(), uns=data.uns.copy()
    )


def subset_pseudobulks(data: ad.AnnData, entry: dict) -> ad.AnnData:
    """Pseudobulks of interleaved subsets of the cells, as replicates of each sample."""
    subsets = [pseudobulk_by_sample(data[i::SUBSETS], entry) for i in range(SUBSETS)]
    pseudobulks = ad.concat(subsets)
    pseudobulks.obs_names_make_unique()
    pseudobulks.var = data.var.copy()
    return pseudobulks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "cohort", help="single-cell cohort in config/datasets.yaml, e.g. wang_2022"
    )
    name = parser.parse_args().cohort
    entry = cohort("single_cell", name)

    data = annotate_cell_populations(ad.read_h5ad(resolve(entry["path"])), entry)

    for key, build in (
        ("pseudobulk", pseudobulk_by_sample),
        ("subsampled_pseudobulk", subset_pseudobulks),
    ):
        pseudobulks = build(data, entry)
        path = resolve(entry[key])
        path.parent.mkdir(parents=True, exist_ok=True)
        pseudobulks.write_h5ad(path)
        print(f"wrote {path}: {pseudobulks.n_obs} pseudobulks")


if __name__ == "__main__":
    main()
