library(ggplot2)
library(dplyr)
library(forcats)
library(extrafont)
library(svglite) 

wellcome_sirius_output <- read.csv("canopus_formula_summary_edit.csv")

filtered_sirius_output <- wellcome_sirius_output %>%
  filter(NPC.superclass.Probability > 0.700)

total_count <- nrow(filtered_sirius_output)

pie_data <- filtered_sirius_output %>%
  mutate(
    NPC.superclass_grouped = fct_lump_prop(f = NPC.superclass, prop = 0.01, other_level = "Other (<1%)")
  ) %>%
  count(NPC.superclass_grouped) %>%
  arrange(desc(n)) %>%
  mutate(
    percentage = n / total_count,
    percent_label = scales::percent(percentage, accuracy = 0.1),
    legend_label = paste0(NPC.superclass_grouped, " (", percent_label, ")"),
        pie_internal_label = ifelse(percentage > 0.07,
                                paste0(NPC.superclass_grouped, "\n", percent_label),
                                NA)
  ) %>%
  mutate(legend_label = factor(legend_label, levels = legend_label))

pie_data_donut <- pie_data %>%
  arrange(desc(legend_label)) %>% 
  mutate(
    end_angle = cumsum(percentage),
    start_angle = lag(end_angle, default = 0),
    label_angle = (start_angle + end_angle) / 2
  )

inner_radius <- 0.6 
outer_radius <- 1 
line_inner_radius <- inner_radius 
line_outer_radius <- outer_radius 

num_classes <- nrow(pie_data)
color_palette <- hcl.colors(num_classes, palette = "Set 3")

wellcomeres_superclass_pie <- ggplot(pie_data_donut) +
  geom_rect(
    aes(xmin = inner_radius, xmax = outer_radius, 
        ymin = start_angle, ymax = end_angle, 
        fill = legend_label)
  ) +
  geom_segment(
    aes(x = line_inner_radius, xend = line_outer_radius, 
        y = end_angle, yend = end_angle),
    color = "white",
    linewidth = 1.2 
  ) +
  
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Distribution of NPC Superclass Hits",
    subtitle = paste0("Based on ", total_count, " total features (Probability > 0.70)"),
    fill = NULL 
  ) +

  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ) +
  
  guides(fill = guide_legend(ncol = 2))

wellcomeres_superclass_piechart

ggsave("wellcomeres_SIRIUS_piechart.svg", plot = wellcomeres_superclass_piechart, 
       width = 10.5, height = 6, units = "in")
