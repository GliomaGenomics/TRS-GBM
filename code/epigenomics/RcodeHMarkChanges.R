library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# -----------------------------
# 1. READ DATA
# -----------------------------
df <- read.csv("Bivalent_status.csv", stringsAsFactors = FALSE)

# -----------------------------
# 2. KEEP VALID PATIENTS + NON-NA JBS
# -----------------------------
all_cols <- colnames(df)

patients <- unique(str_extract(all_cols, "^[^_]+"))

valid_patients <- patients[sapply(patients, function(p) {
  paste0(p, "_Prim") %in% all_cols &&
  paste0(p, "_Rec") %in% all_cols
})]

paired_cols <- unlist(lapply(valid_patients, function(p) {
  c(paste0(p, "_Prim"), paste0(p, "_Rec"))
}))

df_sub <- df %>%
  select(TSS, jbs_tss, all_of(paired_cols)) %>%
  filter(!is.na(jbs_tss))

# -----------------------------
# 3. RESHAPE TO Prim/Rec
# -----------------------------
long_df <- df_sub %>%
  pivot_longer(
    cols = -c(TSS, jbs_tss),
    names_to = c("Patient", "Timepoint"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = value
  )

# -----------------------------
# 4. DEFINE TRANSITIONS
# -----------------------------
is_BR <- function(x) x %in% c("B","R")
is_NA <- function(x) x %in% c("N","A")

long_df <- long_df %>%
  mutate(
    state = case_when(
      is_NA(Prim) & is_BR(Rec) ~ "GAIN",
      is_BR(Prim) & is_NA(Rec) ~ "LOSS",
      is_BR(Prim) & is_BR(Rec) ~ "RETENTION",
      is_NA(Prim) & is_NA(Rec) ~ "NEITHER",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(state)

# -----------------------------
# 5. MAJORITY VOTE PER PROMOTER
# -----------------------------
vote_summary <- long_df %>%
  group_by(TSS, jbs_tss, state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(TSS, jbs_tss) %>%
  mutate(
    max_votes = max(n),
    n_max = sum(n == max_votes)
  ) %>%
  ungroup()

# determine majority
majority_k27 <- vote_summary %>%
  filter(n == max_votes)

majority_k27 <- majority_k27 %>%
  mutate(
    majority_state = ifelse(n_max > 1, "no_majority", state),
    majority_votes = ifelse(n_max > 1, 0, max_votes)
  ) %>%
  select(TSS, jbs_tss, majority_state, majority_votes) %>%
  distinct()

# -----------------------------
# 6. PLOT 1: CONSISTENCY (votes)
# -----------------------------
plot1_df <- majority_k27 %>%
  mutate(
    vote_cat = ifelse(majority_state == "no_majority",
                      "no_majority",
                      as.character(majority_votes))
  ) %>%
  group_by(jbs_tss, vote_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot1_df$jbs_label <- ifelse(plot1_df$jbs_tss == 1, "JBS", "Non-JBS")

# create colors
vote_levels <- sort(unique(plot1_df$vote_cat))
vote_levels <- c("no_majority", setdiff(vote_levels, "no_majority"))

gold_palette <- colorRampPalette(c("#FFF5CC", "#B8860B"))(length(vote_levels) - 1)

cols <- c("grey70", gold_palette)
names(cols) <- vote_levels

p1 <- ggplot(plot1_df,
             aes(x = jbs_label, y = percent, fill = factor(vote_cat, levels = vote_levels))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Majority votes",
    title = "Consistency of promoter transitions"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# -----------------------------
# 7. PLOT 2: TRANSITION TYPES (only majority)
# -----------------------------
plot2_df <- majority_k27 %>%
  filter(majority_state != "no_majority") %>%
  group_by(jbs_tss, majority_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot2_df$jbs_label <- ifelse(plot2_df$jbs_tss == 1, "JBS", "Non-JBS")

p2 <- ggplot(plot2_df,
             aes(x = jbs_label, y = percent, fill = majority_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "GAIN" = "#D55E00",
    "LOSS" = "#0072B2",
    "RETENTION" = "#009E73",
    "NEITHER" = "grey50"
  )) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Transition",
    title = "Majority transition types"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# -----------------------------
# 8. SAVE + PRINT
# -----------------------------
print(p1)
print(p2)

ggsave("Plot_H3K27me3_majority_consistency.png", p1, width = 6, height = 4, dpi = 300)
ggsave("Plot_H3K27me3_majority_transitions.png", p2, width = 6, height = 4, dpi = 300)

# -----------------------------
# FILTER: remove NEITHER majority promoters
# -----------------------------
noNEITHER_df <- majority_k27 %>%
  filter(majority_state != "NEITHER")
  
  plot1_noNEITHER <- noNEITHER_df %>%
  mutate(
    vote_cat = ifelse(majority_state == "no_majority",
                      "no_majority",
                      as.character(majority_votes))
  ) %>%
  group_by(jbs_tss, vote_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot1_noNEITHER$jbs_label <- ifelse(plot1_noNEITHER$jbs_tss == 1, "JBS", "Non-JBS")

# reuse same colour logic as before
vote_levels2 <- sort(unique(plot1_noNEITHER$vote_cat))
vote_levels2 <- c("no_majority", setdiff(vote_levels2, "no_majority"))

gold_palette2 <- colorRampPalette(c("#FFF5CC", "#B8860B"))(length(vote_levels2) - 1)

cols2 <- c("grey70", gold_palette2)
names(cols2) <- vote_levels2

p1_noNEITHER <- ggplot(plot1_noNEITHER,
                   aes(x = jbs_label,
                       y = percent,
                       fill = factor(vote_cat, levels = vote_levels2))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols2) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Majority votes",
    title = "Consistency (excluding NEITHER-majority promoters)"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
  
  plot2_noNEITHER <- noNEITHER_df %>%
  filter(majority_state != "no_majority") %>%
  group_by(jbs_tss, majority_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot2_noNEITHER$jbs_label <- ifelse(plot2_noNEITHER$jbs_tss == 1, "JBS", "Non-JBS")

p2_noNEITHER <- ggplot(plot2_noNEITHER,
                   aes(x = jbs_label,
                       y = percent,
                       fill = majority_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "GAIN" = "#D55E00",
    "LOSS" = "#0072B2",
    "RETENTION" = "#009E73"
  )) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Transition",
    title = "Majority transitions (excluding NEITHER)"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
  
  print(p1_noNEITHER)
print(p2_noNEITHER)

ggsave("Plot_H3K27me3_consistency_noNEITHER.png", p1_noNEITHER, width = 6, height = 4, dpi = 300)
ggsave("Plot_H3K27me3_transitions_noNEITHER.png", p2_noNEITHER, width = 6, height = 4, dpi = 300)

#repeat code block for H3K4me3

# -----------------------------
# 4. DEFINE TRANSITIONS
# -----------------------------
is_BA <- function(x) x %in% c("B","A")
is_NR <- function(x) x %in% c("N","R")

long_df <- long_df %>%
  mutate(
    state = case_when(
      is_NR(Prim) & is_BA(Rec) ~ "GAIN",
      is_BA(Prim) & is_NR(Rec) ~ "LOSS",
      is_BA(Prim) & is_BA(Rec) ~ "RETENTION",
      is_NR(Prim) & is_NR(Rec) ~ "NEITHER",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(state)

# -----------------------------
# 5. MAJORITY VOTE PER PROMOTER
# -----------------------------
vote_summary <- long_df %>%
  group_by(TSS, jbs_tss, state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(TSS, jbs_tss) %>%
  mutate(
    max_votes = max(n),
    n_max = sum(n == max_votes)
  ) %>%
  ungroup()

# determine majority
majority_k4 <- vote_summary %>%
  filter(n == max_votes)

majority_k4 <- majority_k4 %>%
  mutate(
    majority_state = ifelse(n_max > 1, "no_majority", state),
    majority_votes = ifelse(n_max > 1, 0, max_votes)
  ) %>%
  select(TSS, jbs_tss, majority_state, majority_votes) %>%
  distinct()

# -----------------------------
# 6. PLOT 1: CONSISTENCY (votes)
# -----------------------------
plot1_df <- majority_k4 %>%
  mutate(
    vote_cat = ifelse(majority_state == "no_majority",
                      "no_majority",
                      as.character(majority_votes))
  ) %>%
  group_by(jbs_tss, vote_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot1_df$jbs_label <- ifelse(plot1_df$jbs_tss == 1, "JBS", "Non-JBS")

# create colors
vote_levels <- sort(unique(plot1_df$vote_cat))
vote_levels <- c("no_majority", setdiff(vote_levels, "no_majority"))

gold_palette <- colorRampPalette(c("#FFF5CC", "#B8860B"))(length(vote_levels) - 1)

cols <- c("grey70", gold_palette)
names(cols) <- vote_levels

p1 <- ggplot(plot1_df,
             aes(x = jbs_label, y = percent, fill = factor(vote_cat, levels = vote_levels))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Majority votes",
    title = "Consistency of promoter transitions"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# -----------------------------
# 7. PLOT 2: TRANSITION TYPES (only majority)
# -----------------------------
plot2_df <- majority_k4 %>%
  filter(majority_state != "no_majority") %>%
  group_by(jbs_tss, majority_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot2_df$jbs_label <- ifelse(plot2_df$jbs_tss == 1, "JBS", "Non-JBS")

p2 <- ggplot(plot2_df,
             aes(x = jbs_label, y = percent, fill = majority_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "GAIN" = "#D55E00",
    "LOSS" = "#0072B2",
    "RETENTION" = "#009E73",
    "NEITHER" = "grey50"
  )) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Transition",
    title = "Majority transition types"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# -----------------------------
# 8. SAVE + PRINT
# -----------------------------
print(p1)
print(p2)

ggsave("Plot_H3K4me3_majority_consistency.png", p1, width = 6, height = 4, dpi = 300)
ggsave("Plot_H3K4me3_majority_transitions.png", p2, width = 6, height = 4, dpi = 300)

# -----------------------------
# FILTER: remove NEITHER majority promoters
# -----------------------------
noNEITHER_df <- majority_k4 %>%
  filter(majority_state != "NEITHER")
  
  plot1_noNEITHER <- noNEITHER_df %>%
  mutate(
    vote_cat = ifelse(majority_state == "no_majority",
                      "no_majority",
                      as.character(majority_votes))
  ) %>%
  group_by(jbs_tss, vote_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot1_noNEITHER$jbs_label <- ifelse(plot1_noNEITHER$jbs_tss == 1, "JBS", "Non-JBS")

# reuse same colour logic as before
vote_levels2 <- sort(unique(plot1_noNEITHER$vote_cat))
vote_levels2 <- c("no_majority", setdiff(vote_levels2, "no_majority"))

gold_palette2 <- colorRampPalette(c("#FFF5CC", "#B8860B"))(length(vote_levels2) - 1)

cols2 <- c("grey70", gold_palette2)
names(cols2) <- vote_levels2

p1_noNEITHER <- ggplot(plot1_noNEITHER,
                   aes(x = jbs_label,
                       y = percent,
                       fill = factor(vote_cat, levels = vote_levels2))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols2) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Majority votes",
    title = "Consistency (excluding NEITHER-majority promoters)"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
  
  plot2_noNEITHER <- noNEITHER_df %>%
  filter(majority_state != "no_majority") %>%
  group_by(jbs_tss, majority_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(jbs_tss) %>%
  mutate(percent = 100 * n / sum(n))

plot2_noNEITHER$jbs_label <- ifelse(plot2_noNEITHER$jbs_tss == 1, "JBS", "Non-JBS")

p2_noNEITHER <- ggplot(plot2_noNEITHER,
                   aes(x = jbs_label,
                       y = percent,
                       fill = majority_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "GAIN" = "#D55E00",
    "LOSS" = "#0072B2",
    "RETENTION" = "#009E73"
  )) +
  labs(
    x = "",
    y = "% promoters",
    fill = "Transition",
    title = "Majority transitions (excluding NEITHER)"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
  
  print(p1_noNEITHER)
print(p2_noNEITHER)

ggsave("Plot_H3K4me3_consistency_noNEITHER.png", p1_noNEITHER, width = 6, height = 4, dpi = 300)
ggsave("Plot_H3K4me3_transitions_noNEITHER.png", p2_noNEITHER, width = 6, height = 4, dpi = 300)

compute_retention_stats <- function(df, mark_name) {
  
  # Remove uncertain / biologically uninformative
  df_sub <- df %>%
    filter(majority_state != "no_majority")
    # Optional (recommended):
    # %>% filter(majority_state != "NEITHER")
  
  # Binary outcome: RETENTION vs everything else
  df_sub <- df_sub %>%
    mutate(retain = majority_state == "RETENTION")
  
  # Contingency table
  tab <- table(df_sub$jbs_tss, df_sub$retain)
  
  print(paste("Contingency table for", mark_name))
  print(tab)
  
  # Fisher’s exact test
  fisher_res <- fisher.test(tab)
  
  # Logistic regression (effect size)
  glm_fit <- glm(retain ~ jbs_tss, data = df_sub, family = binomial())
  
  # Odds ratio + CI
  OR <- exp(coef(glm_fit)[2])
  CI <- suppressMessages(exp(confint(glm_fit)[2, ]))
  
  # Proportions (for interpretation)
  props <- prop.table(tab, margin = 1)
  
  data.frame(
    mark = mark_name,
    
    JBS_retention = props["1", "TRUE"],
    nonJBS_retention = props["0", "TRUE"],
    
    odds_ratio = OR,
    CI_low = CI[1],
    CI_high = CI[2],
    
    fisher_p = fisher_res$p.value,
    glm_p = summary(glm_fit)$coefficients[2,4]
  )
}

stats_k27 <- compute_retention_stats(majority_k27, "H3K27me3")
stats_k4  <- compute_retention_stats(majority_k4,  "H3K4me3")

stats_all <- bind_rows(stats_k27, stats_k4)

print(stats_all)
write.csv(stats_all, "Stats_RETENTION_H3K4me3_H3K27me3.csv", row.names = FALSE)
