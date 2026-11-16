
library(tidyverse)
library(ggalluvial)
library(patchwork)


# ---------------------------
# Load and clean data
# ---------------------------
df <- read_csv("Bivalent_status.csv") %>%
  filter(!is.na(jbs_tss))
  
  # ---------------------------
# Build all promoter-patient transitions
# ---------------------------
prim_cols <- names(df)[str_detect(names(df), "_Prim$")]
rec_cols  <- names(df)[str_detect(names(df), "_Rec$")]

prim_ids <- str_remove(prim_cols, "_Prim$")
rec_ids  <- str_remove(rec_cols, "_Rec$")

paired_ids <- intersect(prim_ids, rec_ids)

if (length(paired_ids) == 0) {
  stop("No paired primary and recurrent patients were identified.")
}

all_transitions <- map_dfr(paired_ids, function(id) {

  tibble(
    promoter = seq_len(nrow(df)),
    patient  = id,
    Prim     = df[[paste0(id, "_Prim")]],
    Rec      = df[[paste0(id, "_Rec")]],
    jbs_group = ifelse(df$jbs_tss == 1, "JBS", "Non-JBS")
  )

}) %>%

  filter(
    !is.na(Prim),
    !is.na(Rec)
  )
  
  desired_order <- c("N","A","B","R")

all_transitions <- all_transitions %>%
  mutate(
    Prim = factor(Prim, levels = desired_order),
    Rec  = factor(Rec, levels = desired_order)
  )
  
# ---------------------------
# Split JBS vs Non-JBS
# ---------------------------
  
  jbs_df <- all_transitions %>%
  filter(jbs_group == "JBS")

nonjbs_df <- all_transitions %>%
  filter(jbs_group == "Non-JBS")
  
# ---------------------------
# Sankey plotting function
# --------------------------- 
  plot_sankey <- function(data, title_text) {

  sankey_df <- data %>%
    count(Prim, Rec) %>%
    mutate(prop = n / sum(n))

  ggplot(
    sankey_df,
    aes(axis1 = Prim,
        axis2 = Rec,
        y = prop)
  ) +

    geom_alluvium(
      aes(fill = Prim),
      width = 0.25,
      discern = TRUE
    ) +

    geom_stratum(
      width = 0.25,
      fill = "grey80",
      colour = "black"
    ) +

    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum))
    ) +

    scale_x_discrete(limits = c("Prim", "Rec"), expand = c(.1, .1)) +

    scale_fill_manual(
      values = c(
        N = "#67C5CC",
        A = "#E79A95",
        B = "#A8C65A",
        R = "#C7A1E8"
      )
    ) +

    labs(
      title = title_text,
      y = "Proportion"
    ) +

    theme_minimal() +
    theme(legend.position = "none")
}

# ---------------------------
# Plot Sankey diagrams
# ---------------------------
jbs_plot <- plot_sankey(
  jbs_df,
  "JBS promoters"
)

nonjbs_plot <- plot_sankey(
  nonjbs_df,
  "Non-JBS promoters"
)

# Combine into a single figure
combined_sankey <-
  nonjbs_plot +
  jbs_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Promoter state transitions from Primary to Recurrent tumours"
  )

print(combined_sankey)

ggsave(
  "Plot_Sankey_JBS_vs_NonJBS_NOTMAJORITY.png",
  combined_sankey,
  width = 6,
  height = 5,
  dpi = 300
)

bivalency_summary <- all_transitions %>%
  mutate(
    B_transition = case_when(
      Prim == "B" & Rec == "B" ~ "Retained_B",
      Prim == "B" & Rec != "B" ~ "Lost_B",
      Prim != "B" & Rec == "B" ~ "Gained_B",
      TRUE                     ~ "Never_B"
    )
  )
  
  jbs_biv <- all_transitions %>%
  mutate(
    retained_B = ifelse(
      Prim == "B" & Rec == "B",
      "Retained_B",
      "Other"
    )
  )
  
  tab <- table(
  jbs_biv$jbs_group,
  jbs_biv$retained_B
)

tab

fisher.test(tab)

biv_percent <- bivalency_summary %>%
  count(jbs_group, B_transition) %>%
  group_by(jbs_group) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  arrange(jbs_group, B_transition)

print(biv_percent)