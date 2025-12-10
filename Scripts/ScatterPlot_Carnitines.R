#Scatter plot (carnitines)
library(readr)
library(dplyr)
library(ggplot2)

base_dir <- "ENTER YOUR PATH HERE"

in_path  <- file.path(base_dir, "Carnitines_scatterplot.csv")

df <- read_csv(in_path, show_col_types = FALSE)

df_sig <- df %>%
  filter(significant == TRUE)

df_plot <- df_sig %>%
  filter(!is.na(carbon_delta)) %>%     
  mutate(
    carbon_plot = carbon_delta,
    
    chain_group = case_when(
      carbon_delta >=  2 & carbon_delta <=  5 ~ "Short",
      carbon_delta >=  6 & carbon_delta <= 12 ~ "Medium",
      carbon_delta >= 13 & carbon_delta <= 20 ~ "Long",
      TRUE ~ "Other"
    ),
    
    chain_group = factor(chain_group, levels = c("Short", "Medium", "Long", "Other"))
  )

real_vals <- sort(unique(df_plot$carbon_plot))
x_breaks  <- real_vals
x_labels  <- as.character(real_vals)

cat("---- plotting data preview ----\n")
print(
  head(df_plot %>% select(Compound_Name, carbon_delta, delta_mass, carbon_plot, chain_group))
)

# ===== scatter plot =====
p <- ggplot(df_plot, aes(x = carbon_plot, y = delta_mass, fill = chain_group, color = chain_group)) +
  geom_point(
    shape = 21,
    size  = 3.5,
    alpha = 0.7,
    stroke = 0.9
  ) +
  
  scale_y_continuous(breaks = seq(
    floor(min(df_plot$delta_mass, na.rm = TRUE) / 50) * 50,
    ceiling(max(df_plot$delta_mass, na.rm = TRUE) / 50) * 50,
    by = 50
  )) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  
  scale_fill_manual(
    name   = "Chain length",
    values = c(
      "Short"  = "#F4A7B9",
      "Medium" = "#A8D8B9",
      "Long"   = "#C3CDE6"
    ),
    labels = c(
      "Short"  = "Short (2–5)",
      "Medium" = "Medium (6–12)",
      "Long"   = "Long (13–20)"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Short"  = "#B21653",
      "Medium" = "#2E8B57",
      "Long"   = "#3B4EAD"
    ),
    guide = "none"
  ) +
  
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        stroke = 0.9,
        color = c(
          "Short"  = "#B21653",
          "Medium" = "#2E8B57",
          "Long"   = "#3B4EAD"
        )
      )
    )
  ) +
  
  labs(
    x = "ΔCarbon from carnitine",
    y = "ΔMass from carnitine",
    title = " "
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    plot.title   = element_text(size = 18, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text  = element_text(size = 12, color = "black")
  )

print(p)

out_svg <- file.path(base_dir, "Carnitines_scatterplot.svg")

ggsave(
  filename = out_svg,
  plot     = p,
  width    = 8,
  height   = 6
)

cat("SVG saved → ", out_svg, "\n")

