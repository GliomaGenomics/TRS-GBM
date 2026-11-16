library(tidyverse)

# ---------------------------
# Data
# ---------------------------
pie_df <- tibble(
  group = c(
    "Genes consistently found to be\nbivalent at both primary\nand recurrence",
    "Genes consistently found to be\nbivalent at both primary\nand recurrence",
    "All genes",
    "All genes"
  ),
  category = c(
    "JBS", "Non-JBS",
    "JBS", "Non-JBS"
  ),
  percentage = c(
    74, 26,
    8, 92
  )
)

# ---------------------------
# Plot
# ---------------------------
pie_plot <- ggplot(
  pie_df,
  aes(
    x = 1,
    y = percentage,
    fill = category
  )
) +
  geom_col(
    width = 1.4,
    colour = "white",
    linewidth = 1
  ) +
  coord_polar(theta = "y") +
  xlim(0.2, 1.8) +   # makes pies larger within facet
  facet_wrap(
    ~group,
    ncol = 1
  ) +
  geom_text(
    aes(label = paste0(percentage, "%")),
    position = position_stack(vjust = 0.5),
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "JBS" = "#E79A95",
      "Non-JBS" = "grey80"
    )
  ) +
  labs(fill = NULL) +
  theme_void() +
  theme(
    strip.text = element_text(
      size = 13,
      face = "bold",
      margin = margin(b = 2)
    ),
    panel.spacing = unit(0.2, "cm"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

print(pie_plot)

# ---------------------------
# Save
# ---------------------------
ggsave(
  "Plot_PieCharts.pdf",
  pie_plot,
  width = 7,
  height = 7,
  units = "in"
)