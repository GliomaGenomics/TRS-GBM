"""Patient metadata table: one row per paired bulk patient, with the single-cell sheets.

Every EORTC and Discovery pair: clinical record, inclusion and purity flags,
responder type and NES, TRS, and per tumour (primary_, recurrent_) sample,
location, GBMPurity, subtype, JBS axis score, PN:MES score and GBMDeconvoluteR
scores; primary MGMT expression and mutational burden. The workbook adds the
sheets of single_cell_patient_table.py (run it first) and column descriptions.

Writes figures/tables/patient_table_bulk.csv, patient_table_bulk_columns.csv and
patient_table.xlsx.

Example:
    pixi run patient_table
"""

import anndata as ad
import pandas as pd

from common import (FIGURES, cohort, config, natural_order, pnmes_scores, resolve,
                    selected)

TABLES = FIGURES / "tables"

COHORTS = ["EORTC", "Discovery"]

# The clinical record, from the primary sample: (table column, .obs column, description).
CLINICAL = [
    ("centre", "Sample.Source", "Centre (contributing study) of the pair"),
    ("library_type", "LibraryType", "RNA-seq library type"),
    ("gender", "Gender", "Gender (Unknown where not recorded, as for all of EORTC)"),
    ("age_at_diagnosis", "Age", "Age at diagnosis (years)"),
    ("non_surgical_treatment", "Non.Surgical.Treatment", "Non-surgical treatment"),
    ("local_recurrence", "Local_Recurrence..",
     "Recurrence local (TRUE), distant (FALSE) or not recorded (Unknown, TBC)"),
    ("surgical_interval_months", "Surgical.Intervals..months.",
     "Time between the two surgeries (months)"),
    ("overall_survival_months", "survival", "Overall survival (months)"),
    ("vital_status", "Status", "Vital status at last follow-up"),
    ("mgmt_promoter", "MGMT", "MGMT promoter methylation: M methylated, U unmethylated"),
]

# Where each tumour was located, from the primary sample.
LOCATIONS = {"primary": "Location.Primary", "recurrent": "Location.Recurrence"}

# GBMDeconvoluteR cell type scores of each tumour, by compartment: cancer cell states,
# brain resident, immune and vascular cells.
CELL_TYPE_SCORES = ["AC", "MES", "NPC", "OPC",
                    "Astrocyte", "Oligodendrocyte", "Neuron", "Radial glial",
                    "Monocyte", "Macrophage", "Microglia", "DC", "Mast Cells", "T Cells",
                    "NK Cells", "B Cells", "Plasma B",
                    "Mural cell", "Endothelial"]

RESPONDER_TYPES = {"U": "Up", "D": "Down"}


def column_name(cell_type: str) -> str:
    return cell_type.lower().replace(" ", "_")


def mutational_burden() -> pd.DataFrame:
    """Each variant-called primary's mutations per Mb, by patient."""
    tables = []
    for entry in config()["datasets"]["bulk"].values():
        if "mutational_burden" in entry:
            burden = pd.read_csv(resolve(entry["mutational_burden"]))
            patients = burden["Tumor_Sample_Barcode"].str.removesuffix(entry["variant_barcode_suffix"])
            tables.append(burden.set_index(patients)[["total_perMB", "total_perMB_HIGH"]])
    return pd.concat(tables)


def cohort_rows(obs: pd.DataFrame, entry: dict, pn_mes: pd.Series,
                burden: pd.DataFrame) -> pd.DataFrame:
    samples = obs[selected(obs, entry["filters"])]
    samples = samples.assign(patient=samples[entry["patient_column"]].astype(str))
    tumours = {name: samples[samples[entry["stage_column"]] == entry[f"{name}_value"]]
               .set_index("patient") for name in ("primary", "recurrent")}
    primary = tumours["primary"]

    table = pd.DataFrame({"cohort": entry["label"]}, index=primary.index.rename("patient"))
    for column, source, _ in CLINICAL:
        table[column] = primary[source]
    table["meets_clinical_criteria"] = primary[entry["include_column"]].astype(bool)
    table["pair_purity"] = primary[entry["pair_purity_column"]]
    table["passes_pair_purity"] = table["pair_purity"] >= entry["min_purity"]
    table["responder_type"] = primary[entry["responder_column"]].astype(object).map(RESPONDER_TYPES)
    table["responder_nes"] = primary["NES"]
    table["trs"] = primary[entry["trs_column"]]

    for name, rows in tumours.items():
        rows = rows.reindex(table.index)
        table[f"{name}_sample"] = rows["sample"]
        table[f"{name}_location"] = primary[LOCATIONS[name]]
        table[f"{name}_purity"] = rows[entry["purity_column"]]
        table[f"{name}_subtype"] = rows["Subtype"]
        table[f"{name}_jbs_axis_score"] = rows[entry["axis_score_column"]]
        table[f"{name}_pn_mes_score"] = pn_mes.reindex(rows["sample"]).to_numpy()
        if name == "primary":
            table["primary_mgmt_expression"] = primary["MGMT_ExpP"]
            table["primary_mutations_per_mb"] = burden["total_perMB"].reindex(table.index)
            table["primary_high_impact_mutations_per_mb"] = (
                burden["total_perMB_HIGH"].reindex(table.index))
        for cell_type in CELL_TYPE_SCORES:
            table[f"{name}_{column_name(cell_type)}_score"] = rows[cell_type]
    return table.loc[natural_order(table.index)]


def descriptions(min_purity: float) -> dict:
    described = {
        "patient": "Patient",
        "cohort": "Cohort",
        **{column: text for column, _, text in CLINICAL},
        "meets_clinical_criteria": "Matched pair, IDH-wildtype primary GBM, radiotherapy and "
                                   "temozolomide, recurrence not recorded as distant: the pairs "
                                   "fold changes and TRS were computed for",
        "pair_purity": "Lower GBMPurity of the pair's two tumours",
        "passes_pair_purity": f"pair_purity at least {min_purity}, the paired analyses' threshold",
        "responder_type": "Responder subtype, Up or Down",
        "responder_nes": "Responder subtype NES: the normalised enrichment score the subtype "
                         "is assigned from (Up if positive)",
        "trs": "TRS: the patient's log2 fold changes, recurrent over primary, dotted with the "
               "JBS axis loadings",
        "primary_mgmt_expression": "MGMT expression of the primary tumour",
        "primary_mutations_per_mb": "Non-synonymous SNVs per Mb called from the primary's RNA "
                                    "(maftools' nominal 50 Mb; the cohorts were called in "
                                    "different pipeline runs, so compare within a cohort only)",
        "primary_high_impact_mutations_per_mb": "High-impact SNVs per Mb called from the "
                                                "primary's RNA (as above)",
    }
    for name in ("primary", "recurrent"):
        tumour = f"{name} tumour"
        described |= {
            f"{name}_sample": f"Sample of the {tumour}",
            f"{name}_location": f"Location of the {tumour}",
            f"{name}_purity": f"Purity of the {tumour} (GBMPurity)",
            f"{name}_subtype": f"Transcriptional subtype of the {tumour} (GlioVis)",
            f"{name}_jbs_axis_score": f"JBS axis score of the {tumour}: log2(FPKM + 0.01) "
                                      "dotted with the axis loadings",
            f"{name}_pn_mes_score": f"PN:MES score of the {tumour}: log10(1 + counts per "
                                    "10,000) dotted with the Wang (2022) PC_1 loadings",
        }
        described |= {f"{name}_{column_name(c)}_score": f"{c} score of the {tumour} "
                                                          "(GBMDeconvoluteR)"
                      for c in CELL_TYPE_SCORES}
    return described


def bulk_sheet() -> tuple[pd.DataFrame, pd.DataFrame]:
    entries = {name: cohort("bulk", name) for name in COHORTS}
    path, = {resolve(entry["path"]) for entry in entries.values()}
    data = ad.read_h5ad(path, backed="r")
    in_cohorts = pd.concat([selected(data.obs, e["filters"]) for e in entries.values()],
                           axis=1).any(axis=1)
    samples = data[in_cohorts.to_numpy()].to_memory()

    settings = config()["pn_mes_loadings"]
    loadings = pd.read_csv(resolve(settings["path"]), index_col=0)[settings["column"]]
    pn_mes = pnmes_scores(samples.to_df(), loadings)

    obs = samples.obs.assign(sample=samples.obs_names)
    burden = mutational_burden()
    table = pd.concat([cohort_rows(obs, entry, pn_mes, burden) for entry in entries.values()])
    min_purity, = {entry["min_purity"] for entry in entries.values()}
    described = descriptions(min_purity)
    table = table.reset_index()
    columns = pd.DataFrame({"Column": table.columns,
                            "Description": [described[c] for c in table.columns]})
    return table, columns


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    table, columns = bulk_sheet()
    table.to_csv(TABLES / "patient_table_bulk.csv", index=False)
    columns.to_csv(TABLES / "patient_table_bulk_columns.csv", index=False)
    print(f"wrote patient_table_bulk: {len(table)} patients, {table.shape[1]} columns")

    sheets = {"Bulk": table}
    described = [columns.assign(Sheet="Bulk")]
    for name, entry in config()["datasets"]["single_cell"].items():
        sheets[entry["label"]] = pd.read_csv(TABLES / f"patient_table_{name}.csv")
        described.append(pd.read_csv(TABLES / f"patient_table_{name}_columns.csv")
                         .assign(Sheet=entry["label"]))
    with pd.ExcelWriter(TABLES / "patient_table.xlsx") as workbook:
        for sheet, rows in sheets.items():
            rows.to_excel(workbook, sheet_name=sheet, index=False)
        (pd.concat(described)[["Sheet", "Column", "Description"]]
         .to_excel(workbook, sheet_name="Columns", index=False))
    print(f"wrote patient_table.xlsx: sheets {', '.join(sheets)}, Columns")


if __name__ == "__main__":
    main()
