###Spearman correlation analysis for carnitines
library(readr)
library(dplyr)
library(tidyr)

in_base_dir <- "ENTER YOUR PATH HERE"
in_csv      <- file.path(in_base_dir, "3063_features_with_cognitiveresilience.csv")

out_base_dir <- "ENTER YOUR PATH HERE"


dir.create(out_base_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(in_csv, show_col_types = FALSE)
df$cogresilience <- suppressWarnings(as.numeric(df$cogresilience))


start_col <- 4
end_col   <- min(3066, ncol(df))

res_list <- list()

for (i in start_col:end_col) {
  col_name <- colnames(df)[i]
  MoI <- suppressWarnings(as.numeric(df[[i]]))
  
  tmp <- tibble(cogresilience = df$cogresilience, MoI = MoI) %>%
    drop_na(cogresilience, MoI)
  
  ct <- tryCatch(cor.test(tmp$MoI, tmp$cogresilience, method = "spearman", exact = FALSE),
                 error = function(e) NULL)
  if (is.null(ct)) next
  
  res_list[[length(res_list) + 1]] <- tibble(
    feature_col = col_name,
    n_total     = nrow(tmp),
    rho         = unname(ct$estimate),
    p_value     = ct$p.value
  )
}

results_df <- bind_rows(res_list)
if (nrow(results_df) == 0) stop("No valid features for Spearman correlation.")

results_df <- results_df %>%
  mutate(FDR = p.adjust(p_value, method = "BH"),
         significant = !is.na(FDR) & FDR < 0.05)

csv_path <- file.path(out_base_dir, "spearman_all_features_summary.csv")
write_csv(results_df, csv_path)
message(sprintf("Spearman  %s", csv_path))




library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)


in_file <- "ENTER YOUR PATH HERE/spearman_all_features_summary.csv"
out_svg <- "ENTER YOUR PATH HERE/SpearmanFigure.svg"


df <- read_csv(in_file, show_col_types = FALSE)


df_sub <- df %>%
  filter(str_detect(feature_col, regex("Carnitine", ignore_case = TRUE))) %>%
  mutate(
    FDR_safe = ifelse(!is.na(FDR) & FDR > 0, FDR, NA_real_),
    neglogFDR = -log10(FDR_safe),
    is_sig = !is.na(FDR_safe) & FDR_safe < 0.05,
    sign_rho = case_when(
      rho > 0 ~ "pos",
      rho < 0 ~ "neg",
      TRUE ~ "zero"
    ),
    color_group = case_when(
      !is_sig ~ "ns",
      sign_rho == "pos" ~ "pos",
      sign_rho == "neg" ~ "neg",
      TRUE ~ "ns"
    )
  )


col_vals <- c(pos = "#88CCEE", neg = "#EE8866", ns = "grey80")


thr <- 0.05
y_thr <- -log10(thr)


p <- ggplot(df_sub, aes(x = rho, y = neglogFDR, fill = color_group)) +
  geom_hline(yintercept = y_thr, linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(shape = 21, color = "black", size = 1.8, stroke = 0.3, alpha = 0.9, na.rm = TRUE) +
  scale_fill_manual(values = col_vals, breaks = c("pos", "neg", "ns"),
                    labels = c("rho > 0  FDR < 0.05", "rho < 0  FDR < 0.05", "FDR ≥ 0.05")) +
  labs(
    x = "Correlation (Spearman ρ)",
    y = bquote(bold(-log[10](FDR))),
    fill = NULL,
    title = " "
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14)
  )


print(p)

ggsave(out_svg, p, width = 4, height = 5)

p


