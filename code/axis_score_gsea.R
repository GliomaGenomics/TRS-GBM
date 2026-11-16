# Gene sets enriched at low and high JBS axis scores, in the primary and the
# recurrent tumours of each cohort, and the gene sets on which all four agree.
#
# Per cohort and stage, the included pairs (inclusion flag, both tumours at least
# 40% pure): DESeq2 on raw counts of protein-coding genes with a count in every
# sample and at least one read, with the JBS axis score (axis_score_column.py) as
# covariate, and centre where a cohort pools centres. Genes ranked by
# sign(log2 fold change) * -log10(adjusted p); fgseaMultilevel on sets of 15 to
# 5,000 genes, significant at BH <= 0.05. Shared: significant in the same
# direction in all four tests. Seeds carry over between cohorts, so they run in
# the order listed.
#
# Writes the results to the config's axis_score_de and axis_score_gsea paths and
# figures/tables/axis_score_gsea_shared_gene_sets.csv.
#
# Example:
#     pixi run axis_score_gsea

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(SingleCellExperiment)
  library(yaml)
  library(zellkonverter)
})
source("code/common.R")

COHORTS <- c("Discovery", "EORTC")
STAGES <- c("primary", "recurrent")
SIGNIFICANCE <- 0.05
SEED <- 11


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


with_axis_score <- function(data, entry) {
  data$axis_score <- obs_column(data, entry$axis_score_column)
  data
}


tested_genes <- function(data, entry) {
  data <- data[rowData(data)[[entry$gene_type_column]] %in% "protein_coding", ]
  # The file pools cohorts quantified against different gene annotations; a gene
  # absent from a cohort's annotation is missing (NaN) there, not zero.
  data <- data[rowSums(is.na(assay(data, "X"))) == 0, ]
  data[rowSums(assay(data, "X")) > 0, ]
}


fit_axis_score <- function(data, entry) {
  centre <- make.names(entry$centre_column)
  design <- if (length(unique(data[[centre]])) > 1) {
    as.formula(paste("~", centre, "+ axis_score"))
  } else {
    ~ axis_score
  }
  dds <- DESeqDataSetFromMatrix(
    countData = assay(data, "X"),
    colData = colData(data),
    design = design
  )
  dds <- DESeq(dds, quiet = TRUE)
  # Independent filtering keeps the genes that maximise the number significant
  # at this adjusted p.
  results(dds, name = "axis_score", alpha = SIGNIFICANCE)
}


de_table <- function(results, data) {
  results <- as.data.frame(results)
  results$ranking_metric <- sign(results$log2FoldChange) * -log10(results$padj)
  metrics <- c("log2FoldChange", "baseMean", "padj", "ranking_metric")
  table <- merge(as.data.frame(rowData(data)), results[, metrics],
                 by = "row.names", all.y = TRUE)
  rownames(table) <- table$Row.names
  table$Row.names <- NULL
  table[order(table$padj, -table$log2FoldChange), ]
}


stage_de <- function(data, entry, stage) {
  tumours <- tested_genes(with_axis_score(stage_tumours(data, entry, stage), entry), entry)
  message(entry$label, " ", stage, ": ", ncol(tumours), " tumours, ",
          nrow(tumours), " genes tested")
  de_table(fit_axis_score(tumours, entry), tumours)
}


gene_ranking <- function(de, entry) {
  de <- de[complete.cases(de), ]
  de <- de[order(-de$ranking_metric, -abs(de$log2FoldChange)), ]
  setNames(de$ranking_metric, de[[entry$gene_id_column]])
}


enrichment <- function(ranking, gene_sets) {
  results <- fgseaMultilevel(
    pathways = gene_sets,
    stats = ranking,
    minSize = 15,
    maxSize = 5000,
    eps = 0,
    nPermSimple = 100000,
    gseaParam = 1,
    nproc = 8
  )
  results <- as.data.frame(results)[, c("pathway", "NES", "padj")]
  # fgsea leaves NES and p missing for a gene set when too few permutations
  # give an enrichment score of the same sign to normalise by.
  results <- results[complete.cases(results), ]
  # Most significant positive enrichment first, most significant negative last.
  results[order(log10(results$padj) * sign(results$NES)), ]
}


significant <- function(table, direction) {
  enriched <- table$padj <= SIGNIFICANCE &
    (if (direction == "High") table$NES >= 0 else table$NES < 0)
  table$pathway[enriched]
}


shared_gene_sets <- function(tables) {
  high <- Reduce(intersect, lapply(tables, significant, direction = "High"))
  low <- Reduce(intersect, lapply(tables, significant, direction = "Low"))
  shared <- data.frame(
    `Gene set` = c(high, low),
    Direction = c(rep("High JBS axis score", length(high)),
                  rep("Low JBS axis score", length(low))),
    check.names = FALSE
  )
  for (test in names(tables)) {
    index <- match(shared$`Gene set`, tables[[test]]$pathway)
    shared[[paste0("NES, ", test)]] <- tables[[test]]$NES[index]
    shared[[paste0("Adjusted p, ", test)]] <- tables[[test]]$padj[index]
  }
  first <- shared[[paste0("NES, ", names(tables)[1])]]
  shared[order(shared$Direction, -abs(first)), ]
}


write_table <- function(table, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(table, path, ...)
}


main <- function() {
  config <- read_yaml("config/datasets.yaml")
  entries <- config$datasets$bulk[COHORTS]
  gene_sets <- read_gene_sets(config$gsea_gene_sets)
  data <- readH5AD(unique(vapply(entries, function(entry) entry$path, character(1))),
                   reader = "R")

  tables <- list()
  for (stage in STAGES) {
    de <- lapply(entries, stage_de, data = data, stage = stage)
    set.seed(SEED)
    for (name in COHORTS) {
      entry <- entries[[name]]
      write_table(de[[name]], entry$axis_score_de[[stage]], row.names = TRUE)
      gsea <- enrichment(gene_ranking(de[[name]], entry), gene_sets)
      write_table(gsea, entry$axis_score_gsea[[stage]], row.names = FALSE)
      test <- paste(entry$label, stage)
      message(test, ": ", length(significant(gsea, "High")), " high, ",
              length(significant(gsea, "Low")), " low")
      tables[[test]] <- gsea
    }
  }

  shared <- shared_gene_sets(tables)
  message("Shared: ", sum(shared$Direction == "High JBS axis score"), " high, ",
          sum(shared$Direction == "Low JBS axis score"), " low")
  write_table(shared, "figures/tables/axis_score_gsea_shared_gene_sets.csv", row.names = FALSE)
}


main()
