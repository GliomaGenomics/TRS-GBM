# Genes whose expression in the primary tumour, and in the recurrent tumour,
# tracks TRS, and the overlap between cohorts.
#
# Per cohort and stage, the included pairs (inclusion flag, both tumours at least
# 40% pure): DESeq2 on raw counts of protein-coding genes with a count in every
# sample and at least one read, with TRS as covariate (scaled to a maximum
# magnitude of one within the cohort), and centre where a cohort pools centres.
# Significant at BH <= 0.05: Up responder genes rise with TRS, Down responder
# genes fall.
#
# Writes the results, with the ranking metric sign(log2 fold change) *
# -log10(adjusted p), to the config's primary_trs_de and recurrent_trs_de paths,
# and Venn diagrams of the genes the cohorts share, responder_de_overlap_<stage>.
#
# Example:
#     pixi run responder_de_overlap

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(SingleCellExperiment)
  library(yaml)
  library(zellkonverter)
})
source("code/common.R")

COHORTS <- c("Discovery_Stead", "EORTC")
STAGES <- c("primary", "recurrent")
SIGNIFICANCE <- 0.05


# zellkonverter makes .obs column names syntactic (as make.names does), so
# "Min. purity in pair" is read as "Min..purity.in.pair".
obs_column <- function(data, column) {
  colData(data)[[make.names(column)]]
}


stage_tumours <- function(data, entry, stage) {
  in_cohort <- rep(TRUE, ncol(data))
  for (column in names(entry$filters)) {
    in_cohort <- in_cohort & obs_column(data, column) %in% entry$filters[[column]]
  }
  # The inclusion flag sits on the primary sample alone.
  included <- in_cohort &
    obs_column(data, entry$include_column) &
    obs_column(data, entry$pair_purity_column) >= entry$min_purity
  patients <- obs_column(data, entry$patient_column)[which(included)]
  data[, in_cohort &
         obs_column(data, entry$patient_column) %in% patients &
         obs_column(data, entry$stage_column) == entry[[paste0(stage, "_value")]]]
}


with_scaled_trs <- function(data, entry) {
  trs <- obs_column(data, entry$trs_column)
  data$TRS <- trs / max(abs(trs), na.rm = TRUE)
  data
}


tested_genes <- function(data, entry) {
  data <- data[rowData(data)[[entry$gene_type_column]] %in% "protein_coding", ]
  # The file pools cohorts quantified against different gene annotations; a gene
  # absent from a cohort's annotation is missing (NaN) there, not zero.
  data <- data[rowSums(is.na(assay(data, "X"))) == 0, ]
  data[rowSums(assay(data, "X")) > 0, ]
}


fit_trs <- function(data, entry) {
  centre <- make.names(entry$centre_column)
  design <- if (length(unique(data[[centre]])) > 1) {
    as.formula(paste("~", centre, "+ TRS"))
  } else {
    ~ TRS
  }
  dds <- DESeqDataSetFromMatrix(
    countData = assay(data, "X"),
    colData = colData(data),
    design = design
  )
  dds <- DESeq(dds, quiet = TRUE)
  # Independent filtering keeps the genes that maximise the number significant
  # at this adjusted p.
  results(dds, name = "TRS", alpha = SIGNIFICANCE)
}


de_table <- function(results, data) {
  results <- as.data.frame(results)
  results$sig <- ifelse(results$padj <= SIGNIFICANCE, "yes", "no")
  results$expression <- case_when(
    results$log2FoldChange >= 0 & results$padj <= SIGNIFICANCE ~ "Up responders",
    results$log2FoldChange < 0 & results$padj <= SIGNIFICANCE ~ "Down responders",
    TRUE ~ "Unchanged"
  )
  results$ranking_metric <- sign(results$log2FoldChange) * -log10(results$padj)

  metrics <- c("log2FoldChange", "baseMean", "padj", "sig", "expression", "ranking_metric")
  table <- merge(as.data.frame(rowData(data)), results[, metrics],
                 by = "row.names", all.y = TRUE)
  rownames(table) <- table$Row.names
  table$Row.names <- NULL
  table[order(table$padj, -table$log2FoldChange), ]
}


cohort_de <- function(entry, data, stage) {
  data <- stage_tumours(data, entry, stage)
  data <- with_scaled_trs(data, entry)
  data <- tested_genes(data, entry)
  message(entry$label, ": ", ncol(data), " ", stage, " tumours, ", nrow(data), " genes tested")
  de_table(fit_trs(data, entry), data)
}


significant <- function(tables, direction) {
  lapply(tables, function(table) rownames(table)[table$expression == direction])
}


membership <- function(sets, direction) {
  genes <- sort(unique(unlist(sets)))
  in_set <- lapply(sets, function(set) genes %in% set)
  cbind(data.frame(gene = genes, direction = direction),
        as.data.frame(in_set, check.names = FALSE))
}


write_overlap <- function(up, down, name) {
  colours <- RESPONDER_COLOURS
  panel <- venn(up, "Up-regulated", colours$Up, set_size = 4) +
    venn(down, "Down-regulated", colours$Down, set_size = 4) +
    plot_layout(ncol = 2)
  dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
  for (extension in c("png", "pdf")) {
    ggsave(paste0("figures/panels/", name, ".", extension), panel,
           width = 4, height = 3, dpi = 300, bg = "white")
  }

  dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
  write.csv(rbind(membership(up, "Up responders"), membership(down, "Down responders")),
            paste0("figures/source_data/", name, ".csv"), row.names = FALSE)
}


main <- function() {
  datasets <- read_yaml("config/datasets.yaml")$datasets$bulk
  entries <- datasets[COHORTS]
  data <- readH5AD(unique(vapply(entries, function(entry) entry$path, character(1))),
                   reader = "R")

  for (stage in STAGES) {
    tables <- lapply(entries, cohort_de, data = data, stage = stage)
    for (name in COHORTS) {
      path <- entries[[name]][[paste0(stage, "_trs_de")]]
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(tables[[name]], path, row.names = TRUE)
    }

    names(tables) <- vapply(entries, function(entry) entry$label, character(1))
    up <- significant(tables, "Up responders")
    down <- significant(tables, "Down responders")
    for (label in names(tables)) {
      message(label, " ", stage, ": ", length(up[[label]]), " up, ",
              length(down[[label]]), " down")
    }
    write_overlap(up, down, paste0("responder_de_overlap_", stage))
  }
}


main()
