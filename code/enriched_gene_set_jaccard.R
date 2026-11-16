# Gene overlap between the gene sets enriched with TRS in the primary tumours of
# both cohorts.
#
# Jaccard similarity of every pair of the shared gene sets (responder_gsea_overlap.R),
# from their gene lists in the collection, drawn as a heatmap clustered by complete
# linkage on 1 - Jaccard, with each set's direction marked.
#
# Example:
#     pixi run enriched_gene_set_jaccard

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gtable)
  library(pheatmap)
  library(yaml)
})
source("code/common.R")

TITLE <- "Jaccard Similarity of Enriched Gene Sets"
WIDTH_MM <- 55
HEIGHT_MM <- 55


read_shared <- function(path) {
  shared <- read.csv(path)
  # Complete linkage breaks ties in distance by input order, so the gene sets
  # enter the clustering in a fixed order: by name, in the C locale.
  shared[order(shared$pathway, method = "radix"), ]
}


jaccard <- function(gene_sets) {
  n <- length(gene_sets)
  similarity <- matrix(0, n, n, dimnames = list(names(gene_sets), names(gene_sets)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      similarity[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]])) /
        length(union(gene_sets[[i]], gene_sets[[j]]))
    }
  }
  similarity
}


heatmap <- function(similarity, directions, colours) {
  distance <- as.dist(1 - similarity)
  annotation <- data.frame(Direction = factor(directions), row.names = rownames(similarity))
  pheatmap(
    similarity,
    clustering_distance_rows = distance,
    clustering_distance_cols = distance,
    treeheight_row = 10,
    treeheight_col = 10,
    breaks = seq(0, 1, length.out = 100),
    annotation_row = annotation,
    annotation_colors = list(Direction = c("Up responders" = colours$Up,
                                           "Down responders" = colours$Down)),
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    annotation_legend = FALSE,
    annotation_names_row = FALSE,
    angle_col = 0,
    main = TITLE,
    fontsize = TEXT_SIZE,
    silent = TRUE
  )
}


direction_legend <- function(colours) {
  labels <- lapply(c("Up-regulated", "Down-regulated"), textGrob,
                   gp = gpar(fontsize = TEXT_SIZE))
  key <- unit(0.8 * TEXT_SIZE, "points")
  keys <- lapply(c(colours$Up, colours$Down), function(fill) {
    rectGrob(width = key, height = key, gp = gpar(fill = fill, col = NA))
  })
  legend <- gtable(
    widths = unit.c(grobWidth(labels[[1]]), unit(2, "points"), key, unit(8, "points"),
                    grobWidth(labels[[2]]), unit(2, "points"), key),
    heights = unit(1.5 * TEXT_SIZE, "points")
  )
  legend <- gtable_add_grob(legend, c(labels[1], keys[1], labels[2], keys[2]),
                            t = 1, l = c(1, 3, 5, 7))
  legend$vp <- viewport(x = 0, width = sum(legend$widths), just = "left")
  legend
}


raster <- function(fills) {
  rasterGrob(fills, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
}


colour_bar <- function(legend) {
  bar <- legend$children[[1]]
  labels <- legend$children[[2]]$label
  # One colour per interval between breaks, lowest first.
  fills <- bar$gp$fill[seq_along(bar$y)]
  height <- unit(20, "mm")
  width <- unit(2.5, "mm")
  positions <- as.numeric(labels)
  grobTree(
    rasterGrob(matrix(rev(fills), ncol = 1), x = 0, y = unit(1, "npc"),
               width = width, height = height, just = c("left", "top"),
               interpolate = FALSE),
    textGrob(labels, x = width + unit(1, "mm"),
             y = unit(1, "npc") - height * (1 - positions), just = "left",
             gp = gpar(fontsize = TEXT_SIZE))
  )
}


# At print size pheatmap's one rectangle per cell leaves hairline seams, its
# title is bold at 1.3 times the text size and its colour bar is 150 points
# tall. Cells and annotation become rasters, the title and dendrograms are
# restyled, the colour bar is redrawn at 20 mm, and the direction legend goes
# below.
styled <- function(panel, colours) {
  layout <- panel$gtable
  grob <- function(name) which(layout$layout$name == name)
  layout$grobs[[grob("main")]]$gp <- gpar(fontsize = TITLE_SIZE)
  for (tree in c(grob("row_tree"), grob("col_tree"))) {
    layout$grobs[[tree]]$gp <- gpar(lwd = lwd_from_points(0.25))
  }
  layout$grobs[[grob("matrix")]] <- raster(layout$grobs[[grob("matrix")]]$children[[1]]$gp$fill)
  layout$grobs[[grob("row_annotation")]] <- raster(layout$grobs[[grob("row_annotation")]]$gp$fill)
  layout$grobs[[grob("legend")]] <- colour_bar(layout$grobs[[grob("legend")]])

  figure <- gtable(widths = unit(1, "null"),
                   heights = unit.c(unit(1, "null"), unit(2 * TEXT_SIZE, "points")))
  figure <- gtable_add_grob(figure, list(layout, direction_legend(colours)), t = 1:2,
                            l = 1, clip = "off")
  figure$gp <- gpar(fontfamily = FONT_FAMILY)
  figure
}


main <- function() {
  config <- read_yaml("config/datasets.yaml")
  colours <- RESPONDER_COLOURS
  gene_sets <- read_gene_sets(config$gsea_gene_sets)
  shared <- read_shared(config$primary_trs_shared_gene_sets)

  similarity <- jaccard(gene_sets[shared$pathway])
  off_diagonal <- similarity[upper.tri(similarity)]
  message(nrow(similarity), " gene sets; pairwise Jaccard median ",
          signif(median(off_diagonal), 3), ", max ", signif(max(off_diagonal), 3))

  panel <- heatmap(similarity, shared$Direction, colours)
  figure <- styled(panel, colours)
  dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
  ggsave("figures/panels/enriched_gene_set_jaccard.png", figure,
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 300, bg = "white")
  ggsave("figures/panels/enriched_gene_set_jaccard.pdf", figure,
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", device = cairo_pdf)

  plotted <- panel$tree_row$order
  dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
  write.csv(
    cbind(data.frame(gene_set = rownames(similarity)[plotted],
                     direction = shared$Direction[plotted]),
          as.data.frame(similarity[plotted, plotted], check.names = FALSE)),
    "figures/source_data/enriched_gene_set_jaccard.csv",
    row.names = FALSE
  )
}


main()
