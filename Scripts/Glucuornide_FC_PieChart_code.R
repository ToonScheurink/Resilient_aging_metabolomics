#FC plot for potential glucuronide conjugates from MassQL query
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

in_dir  <- "PATH"
in_file <- file.path(in_dir, "2026_newvisits_group_only.csv")

out_dir <- file.path(in_dir, "Log2FC_Lollipop_Styled")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

feature_cols <- c(
  "18404","80624","13002","12878","80654","28338","7251","8377","30763","50784","14627","71530",
  "49824","43351","33563","29896","29124","49797"
)

df <- read_csv(in_file, show_col_types = FALSE) %>%
  mutate(
    cog_group_raw = as.numeric(cogresilience_group),
    CogGroup = case_when(
      cog_group_raw %in% c(1, 2) ~ "Low",
      cog_group_raw %in% c(3, 4) ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CogGroup)) %>%
  mutate(CogGroup = factor(CogGroup, levels = c("Low", "High")))

feature_cols <- intersect(feature_cols, names(df))
stopifnot(length(feature_cols) > 0)

df <- df %>%
  mutate(across(all_of(feature_cols), ~ suppressWarnings(as.numeric(.))))

stats_tbl <- lapply(feature_cols, function(fn) {
  a <- log2(df[[fn]][df$CogGroup == "High"] + 1)
  b <- log2(df[[fn]][df$CogGroup == "Low"]  + 1)
  
  tibble(
    feature = fn,
    log2FC = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)
  )
}) |> bind_rows() |>
  mutate(
    direction = case_when(
      is.na(log2FC) ~ "NA",
      log2FC > 0    ~ "High cognitive resilience",
      log2FC < 0    ~ "Low cognitive resilience",
      TRUE          ~ "No diff"
    )
  )

plt <- stats_tbl %>%
  arrange(log2FC) %>%
  mutate(feature_fac = factor(feature, levels = feature))

p <- ggplot(plt, aes(x = log2FC, y = feature_fac, color = direction)) +
  geom_segment(
    aes(x = 0, xend = log2FC, y = feature_fac, yend = feature_fac),
    linewidth = 0.4, color = "grey80"
  ) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "High cognitive resilience" = "#4B3F8C",
      "Low cognitive resilience"  = "#F05A28",
      "No diff" = "grey60",
      "NA" = "grey70"
    ),
    breaks = c("High cognitive resilience", "Low cognitive resilience")
  ) +
  labs(
    x = "log2 (Fold Change)",
    y = "Feature ID",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "top",
    axis.title.x = element_text(face = "bold", size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 15, margin = margin(r = 10)),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 13, face = "bold")
  )

ggsave(file.path(out_dir, "Lollipop_plot.svg"), p, width = 7, height = 8)
p








#Pie chart for potential glucuronide conjugates
library(readr)
library(dplyr)
library(ggplot2)
library(svglite)  

in_file <- "Path"

df <- read_csv(in_file, show_col_types = FALSE) %>%
  filter(Glucuronidation == "yes") %>%
  count(Direction) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(Direction, "\n", round(percent, 1), "% (n=", n, ")")
  )

print(df)

p <- ggplot(df, aes(x = "", y = n, fill = Direction)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.45),
    size = 4.0,
    fontface = "bold",
    color = "white",
    lineheight = 1.0
  ) +
  scale_fill_manual(values = c("Positive" = "#3586B6", "Negative" = "#E14B5A")) +
  labs(
    title = "Distribution of glucuronide conjugates related to cognitive resilience from MassQL query analysis",
    fill = "Direction"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right"
  )

print(p)

out_dir <- dirname(in_file)
ggsave(
  filename = file.path(out_dir, "piechart.svg"),
  plot = p,
  device = svglite::svglite,  
  width = 6, height = 6       
)
