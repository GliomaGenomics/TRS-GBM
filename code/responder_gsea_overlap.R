# Gene sets enriched with TRS in primary tumours, and in recurrent tumours, and
# the overlap between cohorts.
#
# Per cohort and stage, genes ranked by responder_de_overlap.R's metric;
# fgseaMultilevel on sets of 15 to 5,000 genes, significant at BH <= 0.05: Up
# responder sets have a positive NES, Down responder sets a negative one. Seeds
# carry over between cohorts, so they run in the order listed.
#
# Writes the results to the config's primary_trs_gsea and recurrent_trs_gsea
# paths, the sets significant in the same direction in both cohorts to the
# *_shared_gene_sets paths and figures/tables/responder_gsea_shared_gene_sets.xlsx,
# and Venn diagrams, responder_gsea_overlap_<stage>.
#
# Example:
#     pixi run responder_gsea_overlap

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(patchwork)
  library(writexl)
  library(yaml)
})
source("code/common.R")

COHORTS <- c("Discovery_Stead", "EORTC")
STAGES <- c("primary", "recurrent")
SIGNIFICANCE <- 0.05
SEED <- 11


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


cohort_gsea <- function(entry, gene_sets, stage) {
  de <- read.csv(entry[[paste0(stage, "_trs_de")]], row.names = 1)
  ranking <- gene_ranking(de, entry)
  message(entry$label, " ", stage, ": ", length(ranking), " genes ranked")
  enrichment(ranking, gene_sets)
}


significant <- function(tables, direction) {
  lapply(tables, function(table) {
    enriched <- table$padj <= SIGNIFICANCE &
      (if (direction == "Up responders") table$NES >= 0 else table$NES < 0)
    table$pathway[enriched]
  })
}


shared_gene_sets <- function(up, down, tables) {
  up_counts <- table(unlist(up))
  down_counts <- table(unlist(down))
  shared_up <- names(up_counts[up_counts >= 2])
  shared_down <- names(down_counts[down_counts >= 2])
  shared <- data.frame(
    pathway = c(shared_up, shared_down),
    Direction = c(rep("Up responders", length(shared_up)),
                  rep("Down responders", length(shared_down))),
    Frequency = c(as.vector(up_counts[shared_up]), as.vector(down_counts[shared_down]))
  )
  for (label in names(tables)) {
    scores <- tables[[label]][, c("pathway", "NES", "padj")]
    colnames(scores)[2:3] <- paste0(c("NES_", "padj_"), label)
    shared <- merge(shared, scores, by = "pathway", all.x = TRUE)
  }
  shared[order(abs(shared[[paste0("NES_", names(tables)[1])]]), decreasing = TRUE), ]
}


supplementary_table <- function(shared, labels) {
  table <- data.frame(`Gene set` = shared$pathway, Direction = shared$Direction,
                      check.names = FALSE)
  for (label in labels) {
    table[[paste0("NES, ", label)]] <- shared[[paste0("NES_", label)]]
    table[[paste0("Adjusted p, ", label)]] <- shared[[paste0("padj_", label)]]
  }
  table
}


membership <- function(sets, direction) {
  gene_sets <- sort(unique(unlist(sets)))
  in_set <- lapply(sets, function(set) gene_sets %in% set)
  cbind(data.frame(gene_set = gene_sets, direction = direction),
        as.data.frame(in_set, check.names = FALSE))
}


write_overlap <- function(up, down, name) {
  colours <- RESPONDER_COLOURS
  panel <- venn(up, "Up-regulated", colours$Up) +
    venn(down, "Down-regulated", colours$Down) +
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


stage_gsea <- function(stage, entries, gene_sets, config) {
  set.seed(SEED)
  tables <- lapply(entries, cohort_gsea, gene_sets = gene_sets, stage = stage)
  for (name in COHORTS) {
    path <- entries[[name]][[paste0(stage, "_trs_gsea")]]
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(tables[[name]], path, row.names = FALSE)
  }

  names(tables) <- vapply(entries, function(entry) entry$label, character(1))
  up <- significant(tables, "Up responders")
  down <- significant(tables, "Down responders")
  for (label in names(tables)) {
    message(label, " ", stage, ": ", length(up[[label]]), " up, ",
            length(down[[label]]), " down")
  }

  shared <- shared_gene_sets(up, down, tables)
  path <- config[[paste0(stage, "_trs_shared_gene_sets")]]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(shared, path, row.names = FALSE)
  write_overlap(up, down, paste0("responder_gsea_overlap_", stage))
  supplementary_table(shared, names(tables))
}


main <- function() {
  config <- read_yaml("config/datasets.yaml")
  entries <- config$datasets$bulk[COHORTS]
  gene_sets <- read_gene_sets(config$gsea_gene_sets)

  tables <- lapply(STAGES, stage_gsea, entries = entries, gene_sets = gene_sets,
                   config = config)
  names(tables) <- STAGES
  dir.create("figures/tables", recursive = TRUE, showWarnings = FALSE)
  for (stage in STAGES) {
    write.csv(tables[[stage]],
              paste0("figures/tables/responder_gsea_shared_gene_sets_", stage, ".csv"),
              row.names = FALSE)
  }
  write_xlsx(setNames(tables, tools::toTitleCase(STAGES)),
             "figures/tables/responder_gsea_shared_gene_sets.xlsx")
}


main()
