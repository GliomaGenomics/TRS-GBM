"""Patient TRS, written into the bulk sample metadata.

TRS = a patient's log2 fold changes · the loadings of every axis gene
(`axis_loadings`); for the Discovery patients the axis was fitted on, it is their
PC1 score. Written to both of the patient's samples as `TRS`; patients without fold
changes are left blank.

Example:
    pixi run trs_column
"""

import os

import anndata as ad
import pandas as pd

from common import config, resolve, selected


def axis_loadings(path) -> pd.Series:
    loadings = pd.read_csv(resolve(path), index_col=0).iloc[:, 0]
    # Ensembl IDs are truncated so versioned and unversioned IDs match.
    loadings.index = loadings.index.astype(str).str[:15]
    return loadings


def patient_trs(path, loadings: pd.Series) -> pd.Series:
    log2fc = pd.read_csv(resolve(path), index_col=0)
    log2fc.index = log2fc.index.astype(str).str[:15]
    genes = log2fc.index.intersection(loadings.index)
    print(f"  {len(genes)} of {len(loadings)} axis genes present in {resolve(path).name}")
    return log2fc.loc[genes].T.dot(loadings.loc[genes])


def main() -> None:
    settings = config()
    cohorts = {name: entry for name, entry in settings["datasets"]["bulk"].items()
               if "log2fc" in entry}
    path, = {resolve(entry["path"]) for entry in cohorts.values()}

    data = ad.read_h5ad(path)
    trs = pd.Series(index=data.obs_names, dtype=float)
    for name, entry in cohorts.items():
        print(f"{name}:")
        scores = patient_trs(entry["log2fc"], axis_loadings(entry["axis_loadings"]))
        rows = selected(data.obs, entry["filters"])
        patients = data.obs[entry["patient_column"]].astype(object)[rows]
        trs[rows] = patients.map(scores)
        print(f"  TRS for {patients.isin(scores.index).sum()} of {rows.sum()} samples")

    data.obs["TRS"] = trs
    temporary = path.with_suffix(".h5ad.writing")
    data.write_h5ad(temporary)
    os.replace(temporary, path)
    print(f"wrote TRS for {int(trs.notna().sum())} of {len(trs)} samples to {path}")


if __name__ == "__main__":
    main()
