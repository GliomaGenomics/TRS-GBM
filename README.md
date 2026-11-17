# TRS

Code for Thomas et al., *Glioblastoma exhibits bidirectional adaptation to treatment through JARID2-associated developmental plasticity*.

## Scoring samples on the JBS axis

To place your own samples on the JBS axis, run
[code/jbs_axis_score.py](code/jbs_axis_score.py) on an expression matrix:

```sh
pixi install
pixi run jbs_axis_score --mtx counts.tsv --genes data/gene_metadata/gene_metadata_loadings.csv
```

- **Input:** a tab-separated matrix of raw counts, genes in rows and samples in
  columns, with Ensembl IDs in the first column. Use `-T` if samples are in rows,
  `--sep ","` for a CSV, `--no-Ensembl` for gene symbols and `--no-FPKM` for a
  matrix already in FPKM.
- **Output:** `axis_scores.txt` (or `-o <file>`), one JBS axis score per sample.
- **Method:** counts are FPKM-normalised with the gene lengths in the gene metadata,
  log2-transformed after adding a pseudocount, and multiplied by the loadings of the
  3,746 axis genes. A gene missing from a sample counts as unexpressed. The
  pseudocount is 0.01 for bulk data; use `--pseudocount 1` for single-cell,
  pseudobulk or spatial data. Scores made with different pseudocounts are not
  comparable.

A patient's TRS is the score of their recurrent tumour minus that of their primary.
The paper's bulk scores come from the FPKM matrix of the Tanner et al. (2024)
pipeline, which filters genes by expression in primary and recurrent tumours before
normalising.

### Demonstration

[code/run](code/run) scores synthetic counts for three primary/recurrent pairs
([data/demo/counts.tsv](data/demo/counts.tsv)) and writes their axis scores and TRS
to `results/`:

```sh
pixi run bash code/run
```

## Reproducing the paper's analyses

The scripts regenerate the paper's figure panels and supplementary tables from
processed data.

### Environment

The environment is managed with [pixi](https://pixi.sh) (`pixi install`). Every
script has a pixi task of the same name, run as `pixi run <task> [arguments]`.
`responder_variants` runs in the `maftools` environment (`pixi run -e maftools responder_variants`), which is available on linux-64 and osx-64 only.

### Data

Processed data is not included. Place the input files at the paths given in
[config/datasets.yaml](config/datasets.yaml), under `data/`. The gene annotation in
`data/gene_metadata/` is the only data that is tracked.

| Data                                       | Source                                                                                      |
| ------------------------------------------ | ------------------------------------------------------------------------------------------- |
| Bulk RNA-seq, Discovery cohort             | Tanner et al., *Genome Biol* 2024; 10.1186/s13059-024-03172-3                             |
| Bulk RNA-seq, EORTC cohort                 | Hoogstrate et al., *Cancer Cell* 2023; 10.1016/j.ccell.2023.02.019 (EGA, EGAD00001007860) |
| snRNA-seq, `wang_2022`                   | Wang et al., *Nat Cancer* 2022; 10.1038/s43018-022-00475-x                                |
| snRNA-seq, `nomura_2025`                 | Nomura et al., *Nat Genet* 2025; 10.1038/s41588-025-02167-5                               |
| Multi-region bulk RNA-seq, `mathur_2025` | Mathur et al., *Cell* 2024; 10.1016/j.cell.2023.12.013                                    |
| Visium, `deJong_2025`                    | de Jong et al., *Nature* 2026; 10.1038/s41586-026-11234-8                                 |

The bulk cohorts are read from one AnnData object of raw counts
(`data/bulk/bulk_counts.h5ad`), with tumour purity (GBMPurity), cell-type scores
(GBMDeconvoluteR) and subtype calls (GlioVis) in its sample metadata. The FPKM
matrix and the per-patient log2 fold changes come from the Tanner et al. (2024)
pipeline and are stored separately; TRS and the bulk JBS axis scores are computed
from them, not from the raw counts. Variants were called from the RNA-seq reads and
are read as MAF files.

Each single-cell dataset is read from one AnnData object of raw counts
(`data/single-cell/wang_counts.h5ad`, `data/single-cell/nomura_counts.h5ad`), with
the published cell-type and malignant-state annotations in its cell metadata. The
single-cell analyses hold every cell in memory and need an HPC cluster:
[hpc/run_single_cell.sh](hpc/run_single_cell.sh) runs all of one dataset's analyses
as a SLURM job.

### Running

Outputs go to `figures/panels/<analysis>.png` and `.pdf`, the plotted values to
`figures/source_data/<analysis>.csv`, and tables to `figures/tables/`.

Run these first, in order. Later scripts read what they write.

```sh
pixi run jbs_axis_loadings --cohort Discovery --gene-set jbs_all -o data/loadings/jbs_axis_discovery.csv
pixi run trs_column                     # patient TRS, into bulk_counts.h5ad
pixi run axis_score_column              # bulk JBS axis scores, into bulk_counts.h5ad
pixi run pseudobulk nomura_2025         # and wang_2022
```

| Panel                | Command                                                                                                                 |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| 1C                   | `responder_de_overlap`                                                                                                |
| 1D                   | `responder_gsea_overlap` (after 1C)                                                                                   |
| 1E                   | `enriched_gene_set_jaccard` (after 1D)                                                                                |
| 1F                   | `trs_prediction_feature_selection primary`, then `trs_prediction primary`                                           |
| 1G–1I               | `-e maftools responder_variants`, then `responder_variant_panels`                                                   |
| 1J                   | `trs_vs_axis_score_difference`                                                                                        |
| 2A                   | `axis_score_vs_pnmes`                                                                                                 |
| 2B–2D               | `composition_vs_axis_score nomura_2025`                                                                               |
| 2E, 2F               | `axis_score_composition_regression nomura_2025 --component celltype`                                                  |
| 2G–2I               | `multiregion_axis_score`                                                                                              |
| 2J, 2K               | `visium_axis_score`                                                                                                   |
| 3A                   | `axis_score_polarisation Discovery`                                                                                   |
| 3B                   | `celltype_axis_score_distribution nomura_2025`                                                                        |
| 3C                   | `malignant_state_axis_score_by_stage nomura_2025`                                                                     |
| 3F–3J               | scripts in [code/epigenomics/](code/epigenomics/), run with `Rscript` from the folder holding their input tables |
| Supplementary Fig. 1 | `trs_clinical_correlates`                                                                                             |

| Table                                                              | Command                                                                            |
| ------------------------------------------------------------------ | ---------------------------------------------------------------------------------- |
| Clinical metadata                                                  | `patient_table`; single-cell sheets from `single_cell_patient_table <dataset>` |
| Gene sets enriched with TRS in both cohorts, primary and recurrent | `responder_gsea_overlap`                                                         |
| TRS prediction from primary tumours                                | `trs_prediction primary`                                                         |
| TRS prediction from recurrent tumours                              | `trs_prediction_feature_selection recurrent`, then `trs_prediction recurrent`  |
| Gene sets enriched with the JBS axis score                         | `axis_score_gsea` (after `axis_score_column`)                                  |

### Layout

```text
code/                 one script per analysis; common.py and common.R hold shared helpers
code/run              demonstration run, the Code Ocean entry point
code/epigenomics/     histone-mark analyses
config/datasets.yaml  data paths and per-dataset settings
config/palette.yaml   group colours
data/                 inputs (untracked, except data/gene_metadata/ and data/demo/)
figures/              outputs (untracked)
hpc/                  SLURM job for the single-cell analyses
environment/          Code Ocean environment
```

## Notes

Claude Opus 5.5 (Anthropic) was used to collate and document a rather disparate codebase generated over the course of a PhD, previously spread across several local machines and computing clusters.
