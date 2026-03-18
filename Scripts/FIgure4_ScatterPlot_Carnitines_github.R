#Select carnitines with statistical significance
# ===== Folder path =====
folder_path <- "PATH"

# ===== Input files =====
file_cnt  <- file.path(folder_path, "carnitines_hurdle-count.csv")
file_zero <- file.path(folder_path, "carnitines_hurdle-zero.csv")

# ===== Read data =====
df_cnt  <- read.csv(file_cnt,  check.names = FALSE, stringsAsFactors = FALSE)
df_zero <- read.csv(file_zero, check.names = FALSE, stringsAsFactors = FALSE)

# ===== Ensure p_adj is numeric =====
df_cnt[["p_adj"]]  <- suppressWarnings(as.numeric(df_cnt[["p_adj"]]))
df_zero[["p_adj"]] <- suppressWarnings(as.numeric(df_zero[["p_adj"]]))

# ===== Add significant column (TRUE if p_adj < 0.05, else FALSE) =====
df_cnt[["significant"]]  <- !is.na(df_cnt[["p_adj"]])  & df_cnt[["p_adj"]]  < 0.05
df_zero[["significant"]] <- !is.na(df_zero[["p_adj"]]) & df_zero[["p_adj"]] < 0.05

# ===== Save outputs =====
out_cnt  <- file.path(folder_path, "carnitines_hurdle-count_with_significant.csv")
out_zero <- file.path(folder_path, "carnitines_hurdle-zero_with_significant.csv")

write.csv(df_cnt,  out_cnt,  row.names = FALSE)
write.csv(df_zero, out_zero, row.names = FALSE)

cat("Saved:\n", out_cnt, "\n", out_zero, "\n\n")









#Scatter plot (carnitines) - count
library(readr)
library(dplyr)
library(ggplot2)

base_dir <- "PATH"
in_path  <- file.path(base_dir, "carnitines_hurdle-count_with_significant.csv")

df <- read_csv(in_path, show_col_types = FALSE)

df_sig <- df %>%
  filter(significant == TRUE)

df_plot <- df_sig %>%
  filter(!is.na(carbon_delta)) %>%
  mutate(
    carbon_plot = carbon_delta,
    
    chain_group = case_when(
      carbon_delta == 1                      ~ "C1",
      carbon_delta >=  2 & carbon_delta <= 5  ~ "Short",
      carbon_delta >=  6 & carbon_delta <= 12 ~ "Medium",
      carbon_delta >= 13 & carbon_delta <= 21 ~ "Long",
      carbon_delta >= 22                      ~ "Very Long",
      TRUE ~ "Other"
    ),
    
    chain_group = factor(chain_group,
                         levels = c("C1","Short","Medium","Long","Very Long")),
    
    # shape based on estimate
    estimate_sign = case_when(
      estimate < 0  ~ "Negative",
      estimate > 0  ~ "Positive",
      TRUE ~ NA_character_
    )
  )

# ===== scatter plot =====
p <- ggplot(df_plot, aes(x = carbon_plot,
                         y = delta_mass,
                         fill = chain_group,
                         color = chain_group,
                         shape = estimate_sign)) +
  geom_point(
    size  = 3.5,
    alpha = 0.7,
    stroke = 0.9
  ) +
  
  # shape scale
  scale_shape_manual(
    name = "Association",
    values = c(
      "Negative" = 21,  # circle
      "Positive" = 22   # square
    )
  ) +
  
  # x-axis tick every 1
  scale_x_continuous(
    breaks = seq(
      floor(min(df_plot$carbon_plot, na.rm = TRUE)),
      ceiling(max(df_plot$carbon_plot, na.rm = TRUE)),
      by = 1
    )
  ) +
  
  scale_y_continuous(breaks = seq(
    floor(min(df_plot$delta_mass, na.rm = TRUE) / 50) * 50,
    ceiling(max(df_plot$delta_mass, na.rm = TRUE) / 50) * 50,
    by = 50
  )) +
  
  scale_fill_manual(
    name   = "Chain length",
    values = c(
      "C1"        = "#BDBDBD",
      "Short"     = "#F4A7B9",
      "Medium"    = "#A8D8B9",
      "Long"      = "#C3CDE6",
      "Very Long" = "#E0BBFF"
    ),
    labels = c(
      "C1"        = "C1",
      "Short"     = "Short (C2–C5)",
      "Medium"    = "Medium (C6–C12)",
      "Long"      = "Long (C13–C21)",
      "Very Long" = "Very long (≥C22)"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "C1"        = "#6E6E6E",
      "Short"     = "#B21653",
      "Medium"    = "#2E8B57",
      "Long"      = "#3B4EAD",
      "Very Long" = "#B388EB"
    ),
    guide = "none"
  ) +
  
  # keep fill legend colors visible + show shape legend
  guides(
    fill  = guide_legend(
      order = 1,
      override.aes = list(shape = 21, color = "black")
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(fill = "grey70", color = "black")
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
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text  = element_text(size = 12, color = "black")
  )

print(p)

out_svg <- file.path(base_dir, "Carnitines_scatterplot-count_0317.svg")

ggsave(
  filename = out_svg,
  plot     = p,
  width    = 8,
  height   = 6
)

cat("SVG saved → ", out_svg, "\n")








#Scatter plot (carnitines) - zero
library(readr)
library(dplyr)
library(ggplot2)

base_dir <- "PATH"
in_path  <- file.path(base_dir, "carnitines_hurdle-zero_with_significant.csv")

df <- read_csv(in_path, show_col_types = FALSE)

df_sig <- df %>%
  filter(significant == TRUE)

df_plot <- df_sig %>%
  filter(!is.na(carbon_delta)) %>%
  mutate(
    carbon_plot = carbon_delta,
    
    chain_group = case_when(
      carbon_delta == 1                      ~ "C1",
      carbon_delta >=  2 & carbon_delta <= 5  ~ "Short",
      carbon_delta >=  6 & carbon_delta <= 12 ~ "Medium",
      carbon_delta >= 13 & carbon_delta <= 21 ~ "Long",
      carbon_delta >= 22                      ~ "Very Long",
      TRUE ~ "Other"
    ),
    
    chain_group = factor(chain_group,
                         levels = c("C1","Short","Medium","Long","Very Long")),
    
    # shape based on estimate
    estimate_sign = case_when(
      estimate < 0  ~ "Negative",
      estimate > 0  ~ "Positive",
      TRUE ~ NA_character_
    )
  )

# ===== scatter plot =====
p <- ggplot(df_plot, aes(x = carbon_plot,
                         y = delta_mass,
                         fill = chain_group,
                         color = chain_group,
                         shape = estimate_sign)) +
  geom_point(
    size  = 3.5,
    alpha = 0.7,
    stroke = 0.9
  ) +
  
  # shape scale
  scale_shape_manual(
    name = "Association",
    values = c(
      "Negative" = 21,  # circle
      "Positive" = 22   # square
    )
  ) +
  
  # x-axis tick every 1
  scale_x_continuous(
    breaks = seq(
      floor(min(df_plot$carbon_plot, na.rm = TRUE)),
      ceiling(max(df_plot$carbon_plot, na.rm = TRUE)),
      by = 1
    )
  ) +
  
  scale_y_continuous(breaks = seq(
    floor(min(df_plot$delta_mass, na.rm = TRUE) / 50) * 50,
    ceiling(max(df_plot$delta_mass, na.rm = TRUE) / 50) * 50,
    by = 50
  )) +
  
  scale_fill_manual(
    name   = "Chain length",
    values = c(
      "C1"        = "#BDBDBD",
      "Short"     = "#F4A7B9",
      "Medium"    = "#A8D8B9",
      "Long"      = "#C3CDE6",
      "Very Long" = "#E0BBFF"
    ),
    labels = c(
      "C1"        = "C1",
      "Short"     = "Short (C2–C5)",
      "Medium"    = "Medium (C6–C12)",
      "Long"      = "Long (C13–C21)",
      "Very Long" = "Very long (≥C22)"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "C1"        = "#6E6E6E",
      "Short"     = "#B21653",
      "Medium"    = "#2E8B57",
      "Long"      = "#3B4EAD",
      "Very Long" = "#B388EB"
    ),
    guide = "none"
  ) +
  
  # keep fill legend colors visible + show shape legend
  guides(
    fill  = guide_legend(
      order = 1,
      override.aes = list(shape = 21, color = "black")
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(fill = "grey70", color = "black")
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
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text  = element_text(size = 12, color = "black")
  )

print(p)

out_svg <- file.path(base_dir, "Carnitines_scatterplot-zero_0317.svg")

ggsave(
  filename = out_svg,
  plot     = p,
  width    = 8,
  height   = 6
)

cat("SVG saved → ", out_svg, "\n")