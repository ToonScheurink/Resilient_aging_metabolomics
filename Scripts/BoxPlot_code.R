# Box plots
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

in_file <- "Path"

out_dir <- "Path"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(in_file, show_col_types = FALSE)

# cogresilience_group (3,4 → High / 1,2 → Low)
df <- df %>%
  mutate(
    cog_group = case_when(
      cogresilience_group %in% c(3, 4) ~ "High",
      cogresilience_group %in% c(1, 2) ~ "Low",
      TRUE ~ NA_character_
    )
  )

cols <- c("Low" = "#e377a5", "High" = "#377eb8")

start_col <- #Column number
end_col   <- #Column number

for (i in start_col:end_col) {
  
  col_name <- names(df)[i]
  
  df_temp <- df %>%
    mutate(value = suppressWarnings(as.numeric(.data[[col_name]]))) %>%
    filter(!is.na(cog_group), !is.na(value))
  
  if (length(unique(df_temp$cog_group)) < 2 || nrow(df_temp) < 3) {
    message("Skipped (insufficient groups / insufficient rows): ", col_name)
    next
  }
  
  n_high_nz <- df_temp %>% filter(cog_group == "High", value != 0) %>% nrow()
  n_low_nz  <- df_temp %>% filter(cog_group == "Low",  value != 0) %>% nrow()
  if (n_high_nz < 3 || n_low_nz < 3) {
    message("Skipped (insufficient non-zero count: High=", n_high_nz, ", Low=", n_low_nz, "): ", col_name)
    next
  }
  
  # Mann–Whitney U test
  test_res <- tryCatch(
    wilcox.test(value ~ cog_group, data = df_temp, alternative = "two.sided"),
    error = function(e) list(p.value = NA_real_)
  )
  p_raw <- test_res$p.value
  p_val <- ifelse(is.na(p_raw), NA, signif(p_raw, 3))
  
  p <- ggplot(df_temp, aes(x = cog_group, y = value, fill = cog_group)) +
    geom_boxplot(width = 0.6, fill = NA, aes(color = cog_group), outlier.shape = NA) +
    geom_errorbar(
      aes(
        ymin = after_stat(ymin),
        ymax = after_stat(ymax),
        color = cog_group
      ),
      stat = "boxplot",
      width = 0.25,
      linewidth = 0.6
    ) +
    geom_jitter(aes(color = cog_group), width = 0.15, alpha = 0.7, size = 2) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
    labs(
      title = col_name,
      subtitle = paste0(
        "Mann–Whitney U test p = ", ifelse(is.na(p_val), "NA", p_val),
        "  |  non-zero n (High=", n_high_nz, ", Low=", n_low_nz, ")"
      ),
      x = NULL, y = "Peak area"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y  = element_text(size = 10),
      axis.text.x  = element_text(size = 11, face = "bold", color = "black"),
      legend.position = "none",
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank()
    )
  
  safe_name <- gsub("[/:*?\"<>|]", "_", col_name)
  out_svg <- file.path(out_dir, paste0(safe_name, ".svg"))
  
  tryCatch({
    ggsave(out_svg, p, width = 5, height = 5)
  }, error = function(e) {
    message("SVG save failed: ", col_name)
  })
  
}




  