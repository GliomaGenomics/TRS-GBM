"""Each bulk sample's JBS axis score, written into the bulk sample metadata.

Score = log2(pipeline FPKM + 0.01) · axis loadings, for the cohorts whose config
entry names an FPKM matrix; other samples are left blank. Written to
bulk_counts.h5ad as `JBS axis score`, for the R scripts.

Example:
    pixi run axis_score_column
"""

import os

import anndata as ad
import pandas as pd

from common import config, fpkm_axis_scores, resolve, selected


def main() -> None:
    cohorts = {name: entry for name, entry in config()["datasets"]["bulk"].items()
               if "fpkm" in entry}
    path, = {resolve(entry["path"]) for entry in cohorts.values()}
    column, = {entry["axis_score_column"] for entry in cohorts.values()}

    data = ad.read_h5ad(path)
    scores = pd.Series(index=data.obs_names, dtype=float)
    for name, entry in cohorts.items():
        rows = selected(data.obs, entry["filters"])
        cohort_scores = fpkm_axis_scores(entry).reindex(data.obs_names[rows])
        scores[rows] = cohort_scores
        print(f"{name}: axis scores for {int(cohort_scores.notna().sum())} of {rows.sum()} samples")

    data.obs[column] = scores
    temporary = path.with_suffix(".h5ad.writing")
    data.write_h5ad(temporary)
    os.replace(temporary, path)
    print(f"wrote {column} for {int(scores.notna().sum())} of {len(scores)} samples to {path}")


if __name__ == "__main__":
    main()
