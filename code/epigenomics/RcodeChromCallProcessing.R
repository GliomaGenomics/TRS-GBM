library(tidyverse)

# ============================================================
# 1. AMEND INPUT TSV FILES
# ============================================================

amend_columns <- function(file_path) {

  # Read file
  df <- read.table(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Check required columns are present
  required_cols <- c(
    "z_nb",
    "pvalue",
    "seqnames",
    "start"
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      file_path,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Assign mark status
  #
  # Missing z_nb or pvalue values remain NA rather than being
  # classified as non-significant.
  df <- df %>%
    mutate(
      significant = case_when(
        is.na(z_nb) | is.na(pvalue) ~ NA_character_,
        z_nb > 0 & pvalue < 0.2     ~ "sig",
        TRUE                        ~ "ns"
      )
    )

  # Ensure seqnames has a "chr" prefix
  df <- df %>%
    mutate(
      seqnames = case_when(
        is.na(seqnames) ~ NA_character_,
        str_detect(seqnames, "^chr") ~ seqnames,
        TRUE ~ paste0("chr", seqnames)
      )
    )

  # Create TSS identifier
  #
  # This assumes each interval starts exactly 1,000 bp upstream
  # of the intended TSS.
  df <- df %>%
    mutate(
      TSS = if_else(
        is.na(seqnames) | is.na(start),
        NA_character_,
        paste(seqnames, start + 1000, sep = ":")
      )
    )

  # Check that each non-missing TSS occurs only once
  duplicated_tss <- df %>%
    filter(!is.na(TSS)) %>%
    count(TSS, name = "n") %>%
    filter(n > 1)

  if (nrow(duplicated_tss) > 0) {
    stop(
      "Duplicate TSS values found in ",
      file_path,
      ". Resolve duplicate rows before merging files."
    )
  }

  # Overwrite the original file
  write.table(
    df,
    file = file_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  return(df)
}


# ============================================================
# 2. AMEND ALL TSV FILES IN THE CURRENT DIRECTORY
# ============================================================

folder <- "."

files <- list.files(
  path = folder,
  pattern = "\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No TSV files were found in the current directory.")
}

for (file in files) {

  message("Amending: ", file)

  amend_columns(file)
}


# ============================================================
# 3. READ AND VALIDATE SAMPLE METADATA
# ============================================================

patqc <- read.csv(
  "PatQC.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_patqc_cols <- c(
  "Prefix",
  "Sample_Name",
  "H3K27me3",
  "H3K4me3"
)

missing_patqc_cols <- setdiff(
  required_patqc_cols,
  names(patqc)
)

if (length(missing_patqc_cols) > 0) {
  stop(
    "PatQC.csv is missing required columns: ",
    paste(missing_patqc_cols, collapse = ", ")
  )
}

# Remove leading and trailing spaces from key metadata columns
patqc <- patqc %>%
  mutate(
    Prefix = str_trim(Prefix),
    Sample_Name = str_trim(Sample_Name),
    H3K27me3 = str_to_upper(str_trim(H3K27me3)),
    H3K4me3 = str_to_upper(str_trim(H3K4me3))
  )

# Sample names must be present
if (any(is.na(patqc$Sample_Name) | patqc$Sample_Name == "")) {
  stop("PatQC.csv contains missing or empty Sample_Name values.")
}

# Prefix values must be present
if (any(is.na(patqc$Prefix) | patqc$Prefix == "")) {
  stop("PatQC.csv contains missing or empty Prefix values.")
}

# Sample names must be unique because they become output column names
duplicated_sample_names <- unique(
  patqc$Sample_Name[duplicated(patqc$Sample_Name)]
)

if (length(duplicated_sample_names) > 0) {
  stop(
    "Duplicated Sample_Name values found in PatQC.csv: ",
    paste(duplicated_sample_names, collapse = ", ")
  )
}

# Check QC values
allowed_qc_values <- c("Y", "N", NA_character_)

unexpected_k27_qc <- setdiff(
  unique(patqc$H3K27me3),
  allowed_qc_values
)

unexpected_k4_qc <- setdiff(
  unique(patqc$H3K4me3),
  allowed_qc_values
)

if (length(unexpected_k27_qc) > 0) {
  stop(
    "Unexpected H3K27me3 QC values in PatQC.csv: ",
    paste(unexpected_k27_qc, collapse = ", ")
  )
}

if (length(unexpected_k4_qc) > 0) {
  stop(
    "Unexpected H3K4me3 QC values in PatQC.csv: ",
    paste(unexpected_k4_qc, collapse = ", ")
  )
}


# ============================================================
# 4. FUNCTION TO READ ONE MARK FILE
# ============================================================

process_mark <- function(prefix, sample_name, mark) {

  file_name <- paste0(
    prefix,
    "_",
    mark,
    ".tsv"
  )

  if (!file.exists(file_name)) {
    warning("File not found: ", file_name)
    return(NULL)
  }

  dat <- read.table(
    file_name,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "TSS",
    "significant"
  )

  missing_cols <- setdiff(
    required_cols,
    names(dat)
  )

  if (length(missing_cols) > 0) {
    warning(
      "Required columns missing from ",
      file_name,
      ": ",
      paste(missing_cols, collapse = ", ")
    )

    return(NULL)
  }

  dat <- dat %>%
    select(
      TSS,
      significant
    ) %>%
    filter(!is.na(TSS))

  # Validate the mark-status classifications
  unexpected_status <- setdiff(
    unique(dat$significant),
    c("sig", "ns", NA_character_)
  )

  if (length(unexpected_status) > 0) {
    stop(
      "Unexpected status values in ",
      file_name,
      ": ",
      paste(unexpected_status, collapse = ", ")
    )
  }

  # Check that each TSS occurs only once
  duplicated_tss <- dat %>%
    count(TSS, name = "n") %>%
    filter(n > 1)

  if (nrow(duplicated_tss) > 0) {
    stop(
      "Duplicate TSS values found in ",
      file_name,
      ". Resolve duplicate rows before merging files."
    )
  }

  # Rename the status column using the sample and mark
  output_col <- paste(
    sample_name,
    mark
  )

  dat <- dat %>%
    rename(
      !!output_col := significant
    )

  return(dat)
}


# ============================================================
# 5. FUNCTION TO BUILD A MARK-SPECIFIC STATUS MATRIX
# ============================================================

build_mark_matrix <- function(
    metadata,
    mark,
    qc_column
) {

  mark_dfs <- vector(
    mode = "list",
    length = nrow(metadata)
  )

  for (i in seq_len(nrow(metadata))) {

    # Only process samples that passed QC for this mark
    if (!isTRUE(metadata[[qc_column]][i] == "Y")) {
      next
    }

    prefix <- metadata$Prefix[i]
    sample_name <- metadata$Sample_Name[i]

    mark_dfs[[i]] <- process_mark(
      prefix = prefix,
      sample_name = sample_name,
      mark = mark
    )
  }

  # Remove NULL entries arising from failed QC or missing files
  mark_dfs <- compact(mark_dfs)

  if (length(mark_dfs) == 0) {
    return(NULL)
  }

  # Merge all samples by TSS
  mark_matrix <- reduce(
    mark_dfs,
    full_join,
    by = "TSS"
  )

  return(mark_matrix)
}


# ============================================================
# 6. BUILD H3K4me3 AND H3K27me3 STATUS MATRICES
# ============================================================

df_h3k4me3 <- build_mark_matrix(
  metadata = patqc,
  mark = "H3K4me3",
  qc_column = "H3K4me3"
)

df_h3k27me3 <- build_mark_matrix(
  metadata = patqc,
  mark = "H3K27me3",
  qc_column = "H3K27me3"
)

if (is.null(df_h3k4me3)) {
  stop(
    paste(
      "No H3K4me3 data were successfully processed.",
      "Check PatQC.csv and the input file names."
    )
  )
}

if (is.null(df_h3k27me3)) {
  stop(
    paste(
      "No H3K27me3 data were successfully processed.",
      "Check PatQC.csv and the input file names."
    )
  )
}


# ============================================================
# 7. IDENTIFY SAMPLES WITH BOTH MARKS
# ============================================================

samples_with_both_marks <- patqc %>%
  filter(
    H3K4me3 == "Y",
    H3K27me3 == "Y"
  ) %>%
  pull(Sample_Name)

if (length(samples_with_both_marks) == 0) {
  stop(
    paste(
      "No samples passed QC for both H3K4me3 and H3K27me3.",
      "Bivalent states cannot be calculated."
    )
  )
}


# ============================================================
# 8. CREATE N/A/R/B CLASSIFICATIONS
# ============================================================

# A separate TSS-matched state table will be generated for each
# sample and then joined into one bivalent-status matrix.

bivalent_sample_dfs <- vector(
  mode = "list",
  length = length(samples_with_both_marks)
)

names(bivalent_sample_dfs) <- samples_with_both_marks

for (sample_name in samples_with_both_marks) {

  col_k27 <- paste(
    sample_name,
    "H3K27me3"
  )

  col_k4 <- paste(
    sample_name,
    "H3K4me3"
  )

  # Check that the expected sample columns were successfully created
  if (!col_k4 %in% names(df_h3k4me3)) {
    warning(
      "H3K4me3 column not found for sample: ",
      sample_name
    )

    next
  }

  if (!col_k27 %in% names(df_h3k27me3)) {
    warning(
      "H3K27me3 column not found for sample: ",
      sample_name
    )

    next
  }

  # Match the two marks explicitly by TSS before assigning a state
  sample_marks <- full_join(
    df_h3k4me3 %>%
      select(
        TSS,
        all_of(col_k4)
      ),
    df_h3k27me3 %>%
      select(
        TSS,
        all_of(col_k27)
      ),
    by = "TSS"
  )

  names(sample_marks)[2:3] <- c(
    "k4",
    "k27"
  )

  sample_marks <- sample_marks %>%
    mutate(
      state = case_when(
        k4 == "ns" & k27 == "ns"   ~ "N",
        k4 == "sig" & k27 == "ns"  ~ "A",
        k4 == "ns" & k27 == "sig"  ~ "R",
        k4 == "sig" & k27 == "sig" ~ "B",
        TRUE                        ~ NA_character_
      )
    )

  # Retain the TSS and final state only, then use the sample name
  # as the output column name.
  sample_state <- sample_marks %>%
    select(
      TSS,
      state
    ) %>%
    rename(
      !!sample_name := state
    )

  bivalent_sample_dfs[[sample_name]] <- sample_state
}

# Remove samples that could not be processed
bivalent_sample_dfs <- compact(
  bivalent_sample_dfs
)

if (length(bivalent_sample_dfs) == 0) {
  stop(
    paste(
      "Bivalent states could not be generated for any samples.",
      "Check the input file names and sample columns."
    )
  )
}

# Join all sample-specific state tables by TSS
df_bivalent <- reduce(
  bivalent_sample_dfs,
  full_join,
  by = "TSS"
)


# ============================================================
# 9. READ AND VALIDATE JBS ANNOTATION
# ============================================================

jbs <- read.csv(
  "jbs_status.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"TSS" %in% names(jbs)) {
  stop("jbs_status.csv does not contain a TSS column.")
}

if (anyDuplicated(jbs$TSS)) {

  duplicated_jbs_tss <- jbs %>%
    count(TSS, name = "n") %>%
    filter(n > 1)

  stop(
    paste(
      "Duplicate TSS values were found in jbs_status.csv.",
      "Resolve these duplicates before adding JBS annotations."
    )
  )
}


# ============================================================
# 10. ADD JBS ANNOTATION
# ============================================================

# Use left_join so that only TSSs present in the experimental
# status matrices are retained.

df_h3k4me3 <- df_h3k4me3 %>%
  left_join(
    jbs,
    by = "TSS"
  )

df_h3k27me3 <- df_h3k27me3 %>%
  left_join(
    jbs,
    by = "TSS"
  )

df_bivalent <- df_bivalent %>%
  left_join(
    jbs,
    by = "TSS"
  )


# ============================================================
# 11. WRITE OUTPUT FILES
# ============================================================

write.csv(
  df_h3k4me3,
  file = "H3K4me3_status.csv",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

write.csv(
  df_h3k27me3,
  file = "H3K27me3_status.csv",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

write.csv(
  df_bivalent,
  file = "Bivalent_status.csv",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)


# ============================================================
# 12. REPORT OUTPUT SUMMARY
# ============================================================

message(
  "Analysis complete."
)

message(
  "H3K4me3 samples: ",
  ncol(df_h3k4me3) - ncol(jbs)
)

message(
  "H3K27me3 samples: ",
  ncol(df_h3k27me3) - ncol(jbs)
)

message(
  "Samples with N/A/R/B classifications: ",
  length(bivalent_sample_dfs)
)

message(
  "TSSs in Bivalent_status.csv: ",
  nrow(df_bivalent)
)