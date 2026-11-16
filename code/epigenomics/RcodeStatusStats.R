library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ============================================================
# 1. READ DATA
# ============================================================

readFile <- function(file_path) {

  df <- read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(df)
}

K4  <- readFile("H3K4me3_status.csv")
K27 <- readFile("H3K27me3_status.csv")

datasets <- list(
  H3K27me3 = K27,
  H3K4me3  = K4
)


# ============================================================
# 2. DEFINE ENRICHMENT THRESHOLD
# ============================================================

# A promoter is classified as enriched if at least
# four patients have a significant call for the mark.

threshold <- 4


# ============================================================
# 3. FUNCTION TO COMPUTE STATISTICS
# ============================================================

compute_stats <- function(df, tp_label, threshold) {

  sub <- df %>%
    filter(
      Timepoint == tp_label,
      jbs_tss %in% c(0, 1)
    )

  # Return NA if no data are available
  if (nrow(sub) == 0) {

    return(data.frame(
      Timepoint = tp_label,
      JBS_enriched = NA,
      JBS_not = NA,
      nonJBS_enriched = NA,
      nonJBS_not = NA,
      prop_JBS = NA,
      prop_nonJBS = NA,
      median_JBS = NA,
      median_nonJBS = NA,
      fisher_p = NA,
      odds_ratio = NA,
      wilcox_p = NA
    ))
  }

  # ----------------------------------------------------------
  # Separate JBS and non-JBS promoters
  # ----------------------------------------------------------

  JBS_vals <- sub$n_sig[
    sub$jbs_tss == 1
  ]

  nonJBS_vals <- sub$n_sig[
    sub$jbs_tss == 0
  ]

  # ----------------------------------------------------------
  # Count enriched and non-enriched promoters
  # ----------------------------------------------------------

  JBS_enriched <- sum(
    JBS_vals >= threshold,
    na.rm = TRUE
  )

  JBS_not <- sum(
    JBS_vals < threshold,
    na.rm = TRUE
  )

  nonJBS_enriched <- sum(
    nonJBS_vals >= threshold,
    na.rm = TRUE
  )

  nonJBS_not <- sum(
    nonJBS_vals < threshold,
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # Construct 2x2 contingency table
  # ----------------------------------------------------------

  mat <- matrix(
    c(
      JBS_enriched,
      JBS_not,
      nonJBS_enriched,
      nonJBS_not
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      JBS_status = c("JBS", "Non-JBS"),
      Enrichment = c("Enriched", "Not_enriched")
    )
  )

  # ----------------------------------------------------------
  # Fisher's exact test
  # ----------------------------------------------------------

  fisher_p <- NA_real_
  odds_ratio <- NA_real_

  if (
    all(rowSums(mat) > 0)
  ) {

    fisher_res <- tryCatch(
      fisher.test(mat),
      error = function(e) NULL
    )

    if (!is.null(fisher_res)) {

      fisher_p <- fisher_res$p.value

      odds_ratio <- unname(
        fisher_res$estimate
      )
    }
  }

  # ----------------------------------------------------------
  # Wilcoxon rank-sum test
  # ----------------------------------------------------------

  wilcox_p <- NA_real_

  if (
    length(JBS_vals) > 0 &&
    length(nonJBS_vals) > 0
  ) {

    wilcox_res <- tryCatch(
      wilcox.test(
        JBS_vals,
        nonJBS_vals,
        exact = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(wilcox_res)) {
      wilcox_p <- wilcox_res$p.value
    }
  }

  # ----------------------------------------------------------
  # Calculate proportions and medians
  # ----------------------------------------------------------

  prop_JBS <- if (length(JBS_vals) > 0) {
    mean(JBS_vals >= threshold, na.rm = TRUE)
  } else {
    NA_real_
  }

  prop_nonJBS <- if (length(nonJBS_vals) > 0) {
    mean(nonJBS_vals >= threshold, na.rm = TRUE)
  } else {
    NA_real_
  }

  median_JBS <- if (length(JBS_vals) > 0) {
    median(JBS_vals, na.rm = TRUE)
  } else {
    NA_real_
  }

  median_nonJBS <- if (length(nonJBS_vals) > 0) {
    median(nonJBS_vals, na.rm = TRUE)
  } else {
    NA_real_
  }

  # ----------------------------------------------------------
  # Return statistical summary
  # ----------------------------------------------------------

  data.frame(
    Timepoint = tp_label,

    JBS_enriched = JBS_enriched,
    JBS_not = JBS_not,

    nonJBS_enriched = nonJBS_enriched,
    nonJBS_not = nonJBS_not,

    prop_JBS = prop_JBS,
    prop_nonJBS = prop_nonJBS,

    median_JBS = median_JBS,
    median_nonJBS = median_nonJBS,

    fisher_p = fisher_p,
    odds_ratio = odds_ratio,
    wilcox_p = wilcox_p
  )
}


# ============================================================
# 4. ANALYSE EACH HISTONE MARK
# ============================================================

for (m in names(datasets)) {

  cat("\n====================================\n")
  cat("Processing:", m, "\n")
  cat("====================================\n")

  df <- datasets[[m]]

  mark <- m

  # ----------------------------------------------------------
  # Identify sample columns
  # ----------------------------------------------------------

  # The output files from the preceding pipeline use
  # column names such as:
  #
  # Patient01_Prim H3K27me3
  # Patient01_Rec H3K27me3
  #
  # The regular expression also accepts a full stop
  # between the sample name and mark.

  mark_pattern <- paste0(
    "^(.*)_(Prim|Rec)[ .]",
    mark,
    "$"
  )

  sample_cols <- names(df)[
    str_detect(names(df), mark_pattern)
  ]

  if (length(sample_cols) == 0) {

    warning(
      "No sample columns found for ",
      mark
    )

    next
  }

  # ----------------------------------------------------------
  # Validate required annotation columns
  # ----------------------------------------------------------

  required_cols <- c(
    "TSS",
    "jbs_tss"
  )

  missing_cols <- setdiff(
    required_cols,
    names(df)
  )

  if (length(missing_cols) > 0) {

    stop(
      "Missing columns for ",
      mark,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # Reshape to long format
  # ----------------------------------------------------------

  long_df <- df %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = c(
        "Patient",
        "Timepoint"
      ),
      names_pattern = mark_pattern,
      values_to = "status"
    ) %>%
    mutate(
      status = str_to_lower(
        str_trim(status)
      )
    ) %>%
    filter(
      !is.na(status),
      jbs_tss %in% c(0, 1)
    )

  # ----------------------------------------------------------
  # Validate mark-status values
  # ----------------------------------------------------------

  unexpected_status <- setdiff(
    unique(long_df$status),
    c("sig", "ns")
  )

  if (length(unexpected_status) > 0) {

    stop(
      "Unexpected status values for ",
      mark,
      ": ",
      paste(unexpected_status, collapse = ", ")
    )
  }

  if (nrow(long_df) == 0) {

    warning(
      "No valid observations for ",
      mark
    )

    next
  }

  # ----------------------------------------------------------
  # Collapse duplicate promoter-patient observations
  # ----------------------------------------------------------

  # If duplicate rows occur, a promoter is counted as
  # significant when at least one observation is "sig".
  #
  # This preserves the behaviour of the original script.

  long_df_unique <- long_df %>%
    mutate(
      sig_flag = status == "sig"
    ) %>%
    group_by(
      TSS,
      jbs_tss,
      Patient,
      Timepoint
    ) %>%
    summarise(
      sig_flag = any(sig_flag),
      .groups = "drop"
    )

  # ----------------------------------------------------------
  # Count patients per timepoint
  # ----------------------------------------------------------

  patients_per_tp <- long_df_unique %>%
    distinct(
      Patient,
      Timepoint
    ) %>%
    count(
      Timepoint,
      name = "n_patients"
    )

  # ----------------------------------------------------------
  # Count significant patient calls per TSS
  # ----------------------------------------------------------

  sig_counts <- long_df_unique %>%
    group_by(
      TSS,
      jbs_tss,
      Timepoint
    ) %>%
    summarise(
      n_sig = sum(sig_flag),
      n_informative = n(),
      .groups = "drop"
    ) %>%
    left_join(
      patients_per_tp,
      by = "Timepoint"
    ) %>%
    mutate(
      prop_sig = n_sig / n_informative,
      enriched = n_sig >= threshold
    )

  # ----------------------------------------------------------
  # Build distribution of significant patient calls
  # ----------------------------------------------------------

  plot_df <- sig_counts %>%
    group_by(
      Timepoint,
      jbs_tss,
      n_sig
    ) %>%
    summarise(
      n_tss = n(),
      .groups = "drop"
    ) %>%
    group_by(
      Timepoint,
      jbs_tss
    ) %>%
    mutate(
      percent = 100 * n_tss / sum(n_tss)
    ) %>%
    ungroup()

  if (nrow(plot_df) == 0) {

    warning(
      "No plotting data for ",
      mark
    )

    next
  }

  # ----------------------------------------------------------
  # Run statistical tests
  # ----------------------------------------------------------

  stats_prim <- compute_stats(
    sig_counts,
    "Prim",
    threshold
  )

  stats_rec <- compute_stats(
    sig_counts,
    "Rec",
    threshold
  )

  stats_out <- bind_rows(
    stats_prim,
    stats_rec
  )

  stats_out$mark <- mark

  # ----------------------------------------------------------
  # Print statistical results
  # ----------------------------------------------------------

  cat("\nStatistics for:", mark, "\n")

  print(stats_out)

  # ----------------------------------------------------------
  # Save statistical results
  # ----------------------------------------------------------

  write.csv(
    stats_out,
    paste0(
      "Stats_",
      mark,
      ".csv"
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # Create colour palette
  # ----------------------------------------------------------

  max_sig <- max(
    plot_df$n_sig,
    na.rm = TRUE
  )

  # Vote counts below threshold are grey.
  # Vote counts at or above threshold are gold.

  grey_levels <- 0:min(
    threshold - 1,
    max_sig
  )

  gold_levels <- if (max_sig >= threshold) {
    threshold:max_sig
  } else {
    integer(0)
  }

  grey_vals <- scales::grey_pal()(
    length(grey_levels)
  )

  names(grey_vals) <- as.character(
    grey_levels
  )

  gold_vals <- character(0)

  if (length(gold_levels) > 0) {

    gold_vals <- colorRampPalette(
      c(
        "#FFF5CC",
        "#B8860B"
      )
    )(length(gold_levels))

    names(gold_vals) <- as.character(
      gold_levels
    )
  }

  vals <- c(
    grey_vals,
    gold_vals
  )

  # ----------------------------------------------------------
  # Prepare plotting factors
  # ----------------------------------------------------------

  plot_df <- plot_df %>%
    mutate(
      jbs_label = factor(
        jbs_tss,
        levels = c(0, 1),
        labels = c(
          "Non-JBS",
          "JBS"
        )
      ),

      Timepoint = factor(
        Timepoint,
        levels = c(
          "Prim",
          "Rec"
        )
      ),

      n_sig = factor(
        n_sig,
        levels = 0:max_sig
      )
    )

  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot(
    plot_df,
    aes(
      x = jbs_label,
      y = percent,
      fill = n_sig
    )
  ) +

    geom_col() +

    facet_wrap(
      ~Timepoint
    ) +

    scale_fill_manual(
      values = vals,
      drop = FALSE
    ) +

    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(
        mult = c(0, 0.02)
      )
    ) +

    labs(
      x = "",
      y = "% of promoters",
      fill = "No. significant\npatient calls",

      title = paste0(
        mark,
        " enrichment across patients"
      )
    ) +

    theme_minimal() +

    theme(
      panel.grid = element_blank(),

      axis.line = element_line(
        colour = "black"
      ),

      panel.border = element_blank()
    )

  # ----------------------------------------------------------
  # Print plot
  # ----------------------------------------------------------

  print(p)

  # ----------------------------------------------------------
  # Save plot
  # ----------------------------------------------------------

  ggsave(
    filename = paste0(
      "Plot_",
      mark,
      "_plot.png"
    ),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )

  # ----------------------------------------------------------
  # Print processing summary
  # ----------------------------------------------------------

  cat("\nCompleted:", mark, "\n")

  cat(
    "Number of TSSs analysed:",
    n_distinct(sig_counts$TSS),
    "\n"
  )

  cat(
    "Primary patients:",
    patients_per_tp$n_patients[
      patients_per_tp$Timepoint == "Prim"
    ],
    "\n"
  )

  cat(
    "Recurrent patients:",
    patients_per_tp$n_patients[
      patients_per_tp$Timepoint == "Rec"
    ],
    "\n"
  )

  cat(
    "Enrichment threshold:",
    threshold,
    "significant patient calls\n"
  )
}

cat("\nAll analyses completed.\n")