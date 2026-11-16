library(dplyr)

# Read file
df <- read.csv("Bivalent_status.csv",
               stringsAsFactors = FALSE,
               check.names = FALSE)

# Identify sample columns
prim_cols <- grep("_Prim$", names(df), value = TRUE)
rec_cols  <- grep("_Rec$",  names(df), value = TRUE)

# Summarise each row
gene_votes <- df %>%
  rowwise() %>%
  mutate(
    prim_B = sum(c_across(all_of(prim_cols)) == "B", na.rm = TRUE),
    prim_total = sum(!is.na(c_across(all_of(prim_cols)))),
    
    rec_B = sum(c_across(all_of(rec_cols)) == "B", na.rm = TRUE),
    rec_total = sum(!is.na(c_across(all_of(rec_cols))))
  ) %>%
  ungroup()

# Collapse to gene level
gene_summary <- gene_votes %>%
  group_by(gene) %>%
  summarise(
    # use all rows belonging to the gene
    prim_B = sum(prim_B),
    prim_total = sum(prim_total),

    rec_B = sum(rec_B),
    rec_total = sum(rec_total),

    jbs_tss = as.integer(any(jbs_tss == 1)),

    .groups = "drop"
  ) %>%
  mutate(
    prim_majority_B = prim_B > (prim_total / 2),
    rec_majority_B  = rec_B  > (rec_total / 2)
  )

# Genes with majority B in either primary or recurrent samples
majority_B_genes <- gene_summary %>%
  filter(prim_majority_B | rec_majority_B) %>%
  arrange(desc(jbs_tss), gene)

write.csv(
  majority_B_genes,
  "Genes_majority_B.csv",
  row.names = FALSE
)

print(head(majority_B_genes))