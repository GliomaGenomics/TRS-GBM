library(tidyverse)

# ---------------------------
# Load and clean data
# ---------------------------
df <- read_csv("Bivalent_status.csv")

# Remove NA jbs_tss
df <- df %>% 
  filter(!is.na(jbs_tss))

# ---------------------------
# Identify Prim / Rec columns
# ---------------------------
prim_cols <- names(df)[str_detect(names(df), "_Prim$")]
rec_cols  <- names(df)[str_detect(names(df), "_Rec$")]

# ---------------------------
# Majority vote function
# ---------------------------
get_majority <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(tibble(votes = 0))  # no majority
  }

  tab <- table(x)
  max_votes <- max(tab)
  winners <- names(tab)[tab == max_votes]

  if (length(winners) > 1) {
    tibble(votes = 0)  # tie = no majority
  } else {
    tibble(votes = max_votes)
  }
}

# ---------------------------
# Apply to Prim / Rec
# ---------------------------
summarise_subset <- function(data, cols, label) {
  data %>%
    rowwise() %>%
    mutate(res = list(get_majority(c_across(all_of(cols))))) %>%
    unnest(res) %>%
    ungroup() %>%
    mutate(diagnosis = label)
}

prim_df <- summarise_subset(df, prim_cols, "Prim")
rec_df  <- summarise_subset(df, rec_cols,  "Rec")

# Combine
plot_df <- bind_rows(prim_df, rec_df) %>%
  mutate(
    jbs_group = ifelse(jbs_tss == 1, "JBS", "Non-JBS"),
    vote_strength = votes   # 0 = no majority
  )

# ---------------------------
# Summarise to proportions
# ---------------------------
plot_counts <- plot_df %>%
  group_by(diagnosis, jbs_group, vote_strength) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(diagnosis, jbs_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ---------------------------
# Create ordered grouping for plot
# ---------------------------
plot_counts <- plot_counts %>%
  mutate(
    group = paste(diagnosis, jbs_group, sep = "_"),
    group = factor(group, levels = c(
      "Prim_JBS", "Prim_Non-JBS",
      "Rec_JBS",  "Rec_Non-JBS"
    )),
    vote_strength = factor(
      vote_strength,
      levels = sort(unique(vote_strength))
    )
  )

# ---------------------------
# Create colour palette
# ---------------------------
levels_vs <- levels(plot_counts$vote_strength)

gold_cols <- colorRampPalette(c("#FFD700", "#B8860B"))(length(levels_vs) - 1)
fill_cols <- c("grey70", gold_cols)
names(fill_cols) <- levels_vs

# ---------------------------
# Plot
# ---------------------------
bothMarkVotes<-ggplot(plot_counts, aes(
  x = group,
  y = prop,
  fill = vote_strength
)) +
  geom_bar(stat = "identity", position = "fill") +

  # Labels for layout
  scale_x_discrete(labels = c(
    "Prim_JBS" = "Prim\nJBS",
    "Prim_Non-JBS" = "Prim\nNon-JBS",
    "Rec_JBS" = "Rec\nJBS",
    "Rec_Non-JBS" = "Rec\nNon-JBS"
  )) +

  scale_y_continuous(labels = scales::percent_format()) +

  scale_fill_manual(
    values = fill_cols,
    name = "Vote strength\n(# agreeing samples)"
  ) +

  labs(
    x = "",
    y = "Proportion of promoters",
    title = "Majority vote strength distribution\nPrim vs Rec, split by JBS status"
  ) +

  theme_minimal() +

  # Optional visual separator
  geom_vline(xintercept = 2.5, linetype = "dashed", colour = "grey50")

print(bothMarkVotes)
ggsave("Plot_bothMark_VoteDist.png", bothMarkVotes, width = 6, height = 4, dpi = 300)

#modify majority function to also give winning status (NABR)
get_majority_with_status <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(tibble(status = "no_majority"))
  }

  tab <- table(x)
  max_votes <- max(tab)
  winners <- names(tab)[tab == max_votes]

  if (length(winners) > 1) {
    tibble(status = "no_majority")
  } else {
    tibble(status = winners)
  }
}

#recompute Prim/Rec using this function
summarise_subset_status <- function(data, cols, label) {
  data %>%
    rowwise() %>%
    mutate(res = list(get_majority_with_status(c_across(all_of(cols))))) %>%
    unnest(res) %>%
    ungroup() %>%
    mutate(diagnosis = label)
}

prim_df_status <- summarise_subset_status(df, prim_cols, "Prim")
rec_df_status  <- summarise_subset_status(df, rec_cols,  "Rec")

plot_df_status <- bind_rows(prim_df_status, rec_df_status) %>%
  mutate(
    jbs_group = ifelse(jbs_tss == 1, "JBS", "Non-JBS")
  )

# ---------------------------
# Test enrichment of B state in JBS genes
# ---------------------------

run_B_test <- function(dat, tp){

  tmp <- dat %>%
    filter(diagnosis == tp)

  jbs_B <- sum(tmp$jbs_group == "JBS" & tmp$status == "B")
  jbs_notB <- sum(tmp$jbs_group == "JBS" & tmp$status != "B")

  nonJBS_B <- sum(tmp$jbs_group == "Non-JBS" & tmp$status == "B")
  nonJBS_notB <- sum(tmp$jbs_group == "Non-JBS" & tmp$status != "B")

  mat <- matrix(
    c(
      jbs_B,
      jbs_notB,
      nonJBS_B,
      nonJBS_notB
    ),
    nrow = 2,
    byrow = TRUE
  )

  ft <- fisher.test(mat)

  tibble(
    Timepoint = tp,
    JBS_B = jbs_B,
    JBS_notB = jbs_notB,
    nonJBS_B = nonJBS_B,
    nonJBS_notB = nonJBS_notB,
    prop_JBS_B = jbs_B / (jbs_B + jbs_notB),
    prop_nonJBS_B = nonJBS_B / (nonJBS_B + nonJBS_notB),
    fisher_p = ft$p.value,
    odds_ratio = unname(ft$estimate)
  )
}

B_stats <- bind_rows(
  run_B_test(plot_df_status, "Prim"),
  run_B_test(plot_df_status, "Rec")
)

print(B_stats)

write_csv(
  B_stats,
  "Stats_Bivalent_B_Enrichment_JBS_vs_nonJBS.csv"
)

#summarise proportions
plot_counts_status <- plot_df_status %>%
  group_by(diagnosis, jbs_group, status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(diagnosis, jbs_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

#create plotting group
plot_counts_status <- plot_counts_status %>%
  mutate(
    group = paste(diagnosis, jbs_group, sep = "_"),

    # BAR ORDER: Prim Non-JBS, Prim JBS, Rec Non-JBS, Rec JBS
    group = factor(
      group,
      levels = c(
        "Prim_Non-JBS",
        "Prim_JBS",
        "Rec_Non-JBS",
        "Rec_JBS"
      )
    ),

    # STACK ORDER: bottom -> top
    status = factor(
      status,
      levels = c(
        "no_majority",
        "N",
        "A",
        "B",
        "R"
      )
    )
  )

# colours
status_cols <- c(
  "no_majority" = "grey70",
  "N" = "#66a61e",
  "A" = "#1b9e77",
  "B" = "#d95f02",
  "R" = "#7570b3"
)

bothMarkStatus <- ggplot(
  plot_counts_status,
  aes(
    x = group,
    y = prop,
    fill = status
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +

  scale_x_discrete(
    labels = c(
      "Prim_Non-JBS" = "Prim\nNon-JBS",
      "Prim_JBS"     = "Prim\nJBS",
      "Rec_Non-JBS"  = "Rec\nNon-JBS",
      "Rec_JBS"      = "Rec\nJBS"
    )
  ) +

  scale_y_continuous(
    labels = scales::percent_format()
  ) +

  scale_fill_manual(
    values = status_cols,
    breaks = c(
      "no_majority",
      "N",
      "A",
      "B",
      "R"
    ),
    labels = c(
      "Non-majority",
      "N",
      "A",
      "B",
      "R"
    ),
    name = "Majority status"
  ) +

  labs(
    x = "",
    y = "Proportion of promoters",
    title = "Majority-assigned promoter status distribution\n(grey = no majority)"
  ) +

  theme_minimal() +

  geom_vline(
    xintercept = 2.5,
    linetype = "dashed",
    colour = "grey50"
  )

print(bothMarkStatus)

ggsave(
  "Plot_bothMark_Status.png",
  bothMarkStatus,
  width = 6,
  height = 4,
  dpi = 300
)



# ---------------------------
# DO A QUICK RESET AND THEN PLOT THE DISTRIBUTION OF STATUS (MAJORITY CALLS ONLY) WIHTIN THE HIGHER REOLUTION OF JBS GENES I.E. BROKEN DOWN TO LE50/70 ETC, 
# ---------------------------

library(tidyverse)

# ---------------------------
# Load and clean data
# ---------------------------
df <- read_csv("Bivalent_status.csv")

# Remove NA jbs_tss
df <- df %>%
  filter(!is.na(jbs_tss))

# ---------------------------
# Identify Prim / Rec columns
# ---------------------------
prim_cols <- names(df)[str_detect(names(df), "_Prim$")]
rec_cols  <- names(df)[str_detect(names(df), "_Rec$")]

# ---------------------------
# Majority function (returns status)
# ---------------------------
get_majority_status <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(tibble(status = "no_majority"))
  }

  tab <- table(x)
  max_votes <- max(tab)
  winners <- names(tab)[tab == max_votes]

  if (length(winners) > 1) {
    tibble(status = "no_majority")
  } else {
    tibble(status = winners)
  }
}

# ---------------------------
# Apply to Prim
# ---------------------------
prim_df <- df %>%
  rowwise() %>%
  mutate(res = list(get_majority_status(c_across(all_of(prim_cols))))) %>%
  unnest(res) %>%
  ungroup() %>%
  mutate(diagnosis = "Prim")

# ---------------------------
# Apply to Rec
# ---------------------------
rec_df <- df %>%
  rowwise() %>%
  mutate(res = list(get_majority_status(c_across(all_of(rec_cols))))) %>%
  unnest(res) %>%
  ungroup() %>%
  mutate(diagnosis = "Rec")

# ---------------------------
# Combine + filter to majority only
# ---------------------------
combined_df <- bind_rows(prim_df, rec_df) %>%
  filter(status != "no_majority")

# ---------------------------
# Summarise gene_status vs majority status
# ---------------------------
plot_df <- combined_df %>%
  group_by(diagnosis, gene_status, status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(diagnosis, gene_status) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ---------------------------
# Order factors
# ---------------------------
plot_df <- plot_df %>%
  mutate(
    gene_status = factor(gene_status, levels = c("none", "jbs", "le50", "le70")),
    status = factor(status, levels = c("N", "A", "B", "R")),
    diagnosis = factor(diagnosis, levels = c("Prim", "Rec"))
  )

# ---------------------------
# Plot (Prim left, Rec right)
# ---------------------------
statusPerGeneType<-ggplot(plot_df, aes(
  x = gene_status,
  y = prop,
  fill = status
)) +
  geom_bar(stat = "identity", position = "fill") +

  facet_wrap(~ diagnosis, nrow = 1) +   # ← Prim left, Rec right

  scale_y_continuous(labels = scales::percent_format()) +

  scale_fill_manual(
    values = c(
  		"N" = "#E79A95",  # salmon
  		"A" = "#A8C65A",  # green
  		"B" = "#67C5CC",  # cyan
  		"R" = "#C7A1E8"   # purple
    ),
    name = "Majority status"
  ) +

  labs(
    x = "Gene status",
    y = "Proportion of promoters",
    title = "Majority call distribution within gene_status groups\nPrim (left) vs Rec (right)"
  ) +

  theme_minimal()
  
print(statusPerGeneType)
ggsave("Plot_bothMark_StatusPerGeneType.png", statusPerGeneType, width = 6, height = 4, dpi = 300)


