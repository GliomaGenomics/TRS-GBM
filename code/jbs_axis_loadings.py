"""PC1 gene loadings that define the JBS axis.

Uncentred, unscaled PCA of a cohort's primary-to-recurrent log2 fold changes over
a JBS gene set. Only genes with a fold change in every processing group
(`fold_change_groups`) enter: each group was filtered separately, and a gene its
filter removed is zero. Discovery with jbs_all gives the paper's axis.

Example:
    pixi run jbs_axis_loadings --cohort Discovery --gene-set jbs_all -o data/loadings/jbs_axis_discovery.csv
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from common import cohort, config, resolve, selected


def build_parser() -> argparse.ArgumentParser:
    settings = config()
    cohorts = [name for name, entry in settings["datasets"]["bulk"].items() if "log2fc" in entry]
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cohort", required=True, choices=sorted(cohorts))
    parser.add_argument("--gene-set", required=True, choices=sorted(settings["gene_sets"]))
    parser.add_argument("-o", "--out", type=Path, default=None)
    return parser


def read_log2fc(path) -> pd.DataFrame:
    matrix = pd.read_csv(path, index_col=0)
    # Ensembl IDs are truncated so versioned and unversioned IDs match.
    matrix.index = matrix.index.astype(str).str[:15]
    return matrix


def read_gene_set(path) -> pd.Index:
    with open(path) as handle:
        return pd.Index([line.strip()[:15] for line in handle if line.strip()])


def restrict_to_gene_set(matrix: pd.DataFrame, genes: pd.Index) -> pd.DataFrame:
    return matrix.loc[matrix.index.intersection(genes)]


def drop_incomplete_genes(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix.dropna(axis=0, how="any")


def fold_change_groups(entry: dict) -> list[pd.Index]:
    """The patients of each group the fold changes were built in.

    A group is the union of its filters, each a mapping of metadata column to
    accepted values, all of which must match.
    """
    obs = ad.read_h5ad(resolve(entry["path"]), backed="r").obs
    in_cohort = selected(obs, entry["filters"])
    patients = obs[entry["patient_column"]].astype(object)
    groups = []
    for alternatives in entry["fold_change_groups"]:
        members = pd.concat([selected(obs, filters) for filters in alternatives], axis=1).any(axis=1)
        groups.append(pd.Index(patients[members & in_cohort].unique()))
    return groups


def drop_filtered_genes(matrix: pd.DataFrame, groups: list[pd.Index]) -> pd.DataFrame:
    """Keep the genes with a non-zero fold change for at least one patient of every group."""
    measured = pd.concat(
        [(matrix[matrix.columns.intersection(group)] != 0).any(axis=1) for group in groups], axis=1
    )
    return matrix[measured.all(axis=1)]


def pc1_loadings(matrix: pd.DataFrame) -> pd.Series:
    # prcomp(center = FALSE, scale = FALSE) on the patients x genes matrix; the
    # loadings are the first right singular vector.
    _, _, vt = np.linalg.svd(matrix.T.to_numpy(dtype=float), full_matrices=False)
    return pd.Series(vt[0], index=matrix.index, name="loading").rename_axis("gene")


def main(args: argparse.Namespace) -> None:
    settings = config()
    entry = cohort("bulk", args.cohort)
    genes = read_gene_set(resolve(settings["gene_sets"][args.gene_set]))
    matrix = read_log2fc(resolve(entry["log2fc"]))

    matrix = restrict_to_gene_set(matrix, genes)
    print(f"{len(genes)} genes in {args.gene_set}, {matrix.shape[0]} present in {args.cohort}")
    matrix = drop_incomplete_genes(matrix)
    matrix = drop_filtered_genes(matrix, fold_change_groups(entry))
    print(f"{matrix.shape[0]} of them have a fold change in every group")

    out = args.out or resolve(entry["loadings"][args.gene_set])
    out.parent.mkdir(parents=True, exist_ok=True)
    pc1_loadings(matrix).to_csv(out)
    print(f"PCA on {matrix.shape[1]} patients x {matrix.shape[0]} genes -> {out}")


if __name__ == "__main__":
    main(build_parser().parse_args())
