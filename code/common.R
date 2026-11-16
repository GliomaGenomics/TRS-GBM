# Helpers and figure style shared by the R scripts. Scripts load it with
# source("code/common.R") and are run from the repository root.

# Nature figure style. A panel is drawn at its final printed size, with Arial
# text of 5 to 7 pt and lines of at least 0.25 pt.
FONT_FAMILY <- "Arial"
TITLE_SIZE <- 7
TEXT_SIZE <- 6

# Up and Down responders, and the genes and gene sets that rise or fall with TRS:
# config/palette.yaml's responder colours. R does not read that file, whose
# non-ASCII keys fail to parse outside a UTF-8 locale.
RESPONDER_COLOURS <- list(Up = "#FFBA49", Down = "#63CCCA")

# grid line widths are in units of 1/96 inch; a point is 1/72 inch.
lwd_from_points <- function(points) {
  points * 96 / 72
}


read_gene_sets <- function(path) {
  gene_sets <- fgsea::gmtPathways(path)
  # Four curated sets are named "CAHOY <cell type>"; the others use underscores.
  names(gene_sets) <- gsub("^CAHOY ", "CAHOY_", names(gene_sets))
  # Empty fields in the GMT would otherwise be read as a gene named "".
  lapply(gene_sets, function(genes) genes[genes != ""])
}


venn <- function(sets, title, colour, ...) {
  plot <- ggVennDiagram::ggVennDiagram(sets, label_alpha = 0, label = "count", ...) +
    ggplot2::scale_fill_gradient(low = "white", high = colour, limits = c(0, NA)) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .3)) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(b = 10)))
  # Move the cohort names clear of the circles.
  plot$layers[[3]]$data$X <- plot$layers[[3]]$data$X + c(5, 5)
  plot$layers[[3]]$data$Y <- plot$layers[[3]]$data$Y + c(-5, 5)
  plot
}
