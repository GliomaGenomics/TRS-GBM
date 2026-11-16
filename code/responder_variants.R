# Somatic variants in primary tumours, compared between Up and Down responders.
#
# Single-nucleotide calls from each primary tumour's RNA (Mutect2 tumour-only,
# VEP) in the included pairs, kept at depth >= 5, at least 2 alternative and 2
# reference reads, and purity-adjusted allele fraction >= 0.05; maftools'
# non-synonymous classes only. Mutational burden per Mb (maftools' nominal 50 Mb);
# per gene, Fisher's exact test of Up against Down (mafCompare, mutated in at
# least 2 samples); genes at p < 0.05 compared between cohorts; SRSF11 variants
# tabulated with the protein's domains.
#
# Writes the tables to the config's paths and the Venn diagram
# variant_responder_overlap; responder_variant_panels.py draws the other panels.
#
# Example:
#     pixi run -e maftools responder_variants

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(maftools)
  library(SingleCellExperiment)
  library(yaml)
  library(zellkonverter)
})
source("code/common.R")

COHORTS <- c("EORTC", "Discovery_Stead")

MIN_DEPTH <- 5
MIN_ALLELE_READS <- 2
MIN_PURITY_ADJUSTED_AF <- 0.05

MIN_MUTATED <- 2
SIGNIFICANCE <- 0.05
GENE <- "SRSF11"

# Most severe first. Classes outside the list (none among these SNVs) sort last.
SEVERITY <- c("Translation_Start_Site", "Splice_Site", "Nonstop_Mutation",
              "Nonsense_Mutation", "Missense_Mutation", "RNA", "3'Flank", "3'UTR",
              "5'Flank", "5'UTR", "IGR", "Intron", "Silent")


# maftools reports its progress on the console; the analysis has no use for it.
quietly <- function(expr) {
  invisible(capture.output(result <- suppressMessages(expr)))
  result
}


# zellkonverter's R reader makes .obs names syntactic, as make.names does.
obs_column <- function(obs, column) {
  obs[[make.names(column)]]
}


primary_records <- function(entry) {
  obs <- as.data.frame(colData(readH5AD(entry$path, reader = "R", use_hdf5 = TRUE)))
  keep <- obs_column(obs, entry$stage_column) == entry$primary_value
  for (column in names(entry$filters)) {
    keep <- keep & obs_column(obs, column) %in% entry$filters[[column]]
  }
  obs <- obs[keep, ]
  data.table(
    Tumor_Sample_Barcode = paste0(obs_column(obs, entry$patient_column),
                                  entry$variant_barcode_suffix),
    included = obs_column(obs, entry$include_column) &
      obs_column(obs, entry$pair_purity_column) >= entry$min_purity,
    purity = obs_column(obs, entry$purity_column),
    Responder = as.character(obs_column(obs, entry$responder_column))
  )
}


filtered_variants <- function(entry) {
  records <- primary_records(entry)
  kept <- records[included == TRUE]

  # Read first with every variant class kept, so the read filters see all calls.
  classes <- fread(entry$variants, select = "Variant_Classification")$Variant_Classification
  classes <- names(sort(table(classes), decreasing = TRUE))
  maf <- quietly(read.maf(entry$variants, clinicalData = records, vc_nonSyn = classes,
                          verbose = FALSE))
  maf <- quietly(subsetMaf(maf, tsb = kept$Tumor_Sample_Barcode))

  calls <- maf@data
  purity <- kept$purity[match(calls$Tumor_Sample_Barcode, kept$Tumor_Sample_Barcode)]
  calls$t_af_purity <- calls$t_alt_count / calls$t_depth / purity
  calls <- calls[t_depth >= MIN_DEPTH &
                   t_af_purity >= MIN_PURITY_ADJUSTED_AF &
                   t_alt_count >= MIN_ALLELE_READS &
                   t_ref_count >= MIN_ALLELE_READS]

  # Read again with maftools' default classes, which sets the synonymous and
  # non-coding calls aside.
  quietly(read.maf(calls, clinicalData = kept, verbose = FALSE))
}


burden_table <- function(maf) {
  responders <- getClinicalData(maf)[, .(Tumor_Sample_Barcode, Responder)]
  burden <- quietly(tmb(maf = maf, plotType = NA))
  high <- quietly(tmb(maf = subsetMaf(maf, query = "IMPACT == 'HIGH'"), plotType = NA))
  burden <- merge(burden, high, by = "Tumor_Sample_Barcode", suffixes = c("", "_HIGH"))
  merge(burden, responders, by = "Tumor_Sample_Barcode")[
    , .(Tumor_Sample_Barcode, Responder, total, total_perMB, total_HIGH, total_perMB_HIGH)]
}


most_severe_per_gene <- function(maf) {
  calls <- maf@data
  calls$Variant_Classification <- factor(calls$Variant_Classification, levels = SEVERITY)
  calls <- calls[order(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification)]
  maf@data <- calls[!duplicated(calls[, .(Tumor_Sample_Barcode, Hugo_Symbol)])]
  maf
}


by_responder <- function(maf) {
  records <- getClinicalData(maf)
  list(
    Up = quietly(subsetMaf(maf, tsb = records[Responder == "U"]$Tumor_Sample_Barcode)),
    Down = quietly(subsetMaf(maf, tsb = records[Responder == "D"]$Tumor_Sample_Barcode))
  )
}


compare_responders <- function(groups) {
  quietly(mafCompare(m1 = groups$Up, m2 = groups$Down, m1Name = "U", m2Name = "D",
                     minMut = MIN_MUTATED))$results
}


gene_variants <- function(groups, labels) {
  rows <- list()
  for (name in names(groups)) {
    for (subtype in c("Up", "Down")) {
      calls <- groups[[name]][[subtype]]@data[Hugo_Symbol == GENE]
      rows[[length(rows) + 1]] <- calls[, .(
        Responder = subtype,
        Cohort = labels[[name]],
        Tumor_Sample_Barcode,
        Variant_Classification = as.character(Variant_Classification),
        HGVSp_Short,
        # The residue is the first number in the protein change, as maftools reads it.
        Protein_position = as.integer(regmatches(HGVSp_Short, regexpr("[0-9]+", HGVSp_Short)))
      )]
    }
  }
  rbindlist(rows)
}


# maftools draws a gene on its longest annotated transcript.
gene_protein <- function() {
  domains <- readRDS(system.file("extdata", "protein_domains.RDs", package = "maftools"))
  domains <- domains[HGNC == GENE]
  domains <- domains[aa.length == max(aa.length)]
  domains[, .(Gene = HGNC, refseq.ID, aa.length, Start, End, Label)]
}


membership <- function(genes) {
  all_genes <- sort(unique(unlist(genes)))
  cbind(data.frame(gene = all_genes),
        as.data.frame(lapply(genes, function(set) all_genes %in% set), check.names = FALSE))
}


write_table <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(table, path)
}


main <- function() {
  config <- read_yaml("config/datasets.yaml")
  datasets <- config$datasets$bulk

  significant <- list()
  groups <- list()
  labels <- list()
  for (name in COHORTS) {
    entry <- datasets[[name]]
    labels[[name]] <- entry$label
    maf <- filtered_variants(entry)
    records <- getClinicalData(maf)
    message(entry$label, ": ", nrow(records), " primary tumours (",
            sum(records$Responder == "U"), " Up, ", sum(records$Responder == "D"), " Down), ",
            nrow(maf@data), " non-synonymous variants")
    write_table(maf@data[, .(Tumor_Sample_Barcode, Hugo_Symbol, Gene,
                             Variant_Classification = as.character(Variant_Classification),
                             IMPACT, HGVSp_Short)],
                entry$filtered_variants)
    write_table(burden_table(maf), entry$mutational_burden)

    groups[[name]] <- by_responder(most_severe_per_gene(maf))
    results <- compare_responders(groups[[name]])
    write_table(results, entry$responder_variant_association)
    significant[[entry$label]] <- results[pval < SIGNIFICANCE]$Hugo_Symbol
    message(entry$label, ": ", nrow(results), " genes tested, ",
            length(significant[[entry$label]]), " at p < ", SIGNIFICANCE,
            ", lowest adjusted p ", signif(min(results$adjPval), 2))
  }

  write_table(gene_variants(groups, labels), config$srsf11_variants)
  write_table(gene_protein(), config$srsf11_protein)

  # The first set is drawn below: Discovery under EORTC, as in the other cohort Venns.
  genes <- significant[c("Discovery", "EORTC")]
  message("shared at p < ", SIGNIFICANCE, ": ", paste(Reduce(intersect, genes), collapse = ", "))
  panel <- venn(genes, paste("Genes at p <", SIGNIFICANCE), "grey75", set_size = 4)
  dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
  for (extension in c("png", "pdf")) {
    ggsave(paste0("figures/panels/variant_responder_overlap.", extension), panel,
           width = 2, height = 3, dpi = 300, bg = "white")
  }
  dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
  write.csv(membership(genes), "figures/source_data/variant_responder_overlap.csv",
            row.names = FALSE)
}


main()
