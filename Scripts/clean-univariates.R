#packages
library(data.table)
library(dplyr)
library(tibble)
library(ggpubr)
library(stringr)
library(tidyr)
library(caret)
library(mixOmics)
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(MuMIn)
library(performance)
library(pscl)
library(purrr)
library(broom)
library(readxl)

#data
df1 <- readRDS("clean_data_wellcome.rds")
carnitines <- read_excel("Carnitine_table.xlsx")
carnitines <- carnitines %>%
  rename(scan_number = SCAN)
carnitines <- carnitines %>%
  mutate(scan_number = as.character(scan_number))


metabs <- readRDS("ranked_metabolites.rds")

df1$host_age<-as.numeric(df1$host_age)
df1$host_body_mass_index<-as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl<-as.numeric(df1$physresilience_pctl)
df1$physresilience<-as.numeric(df1$physresilience)
df1$cogresilience<-as.numeric(df1$cogresilience)
df1$cogresilience_pctl<-as.numeric(df1$cogresilience_pctl)
df1$cogresilience_group <- factor(df1$cogresilience_group)
df1$cogresilience_z <- scale(df1$cogresilience)

df <- df1[, colnames(df1) %in% c(
  "filename",
  "cogresilience_group",
  "cogresilience_z",
  metabs$Metabolite
)]

colnames(df) <- make.names(colnames(df), unique = TRUE)

#prep for univariate
meta_cols <- c("filename", "cogresilience_z", "cogresilience_group")
metabolite_cols <- setdiff(names(df), meta_cols)
length(metabolite_cols)   # should be ~3000+


#Loop for hyrdle model
fit_two_part <- function(metab) {
  
  dat <- df %>%
    dplyr::select(cogresilience_z, dplyr::all_of(metab)) %>%
    dplyr::rename(abundance = dplyr::all_of(metab)) %>%
    dplyr::mutate(present = abundance > 0)
  
  out <- list()
  
  ## Part 1: presence / absence
  m1 <- tryCatch(
    glm(present ~ cogresilience_z, data = dat, family = binomial),
    error = function(e) NULL
  )
  
  if (!is.null(m1)) {
    out$zero <- broom::tidy(m1) |>
      dplyr::filter(term == "cogresilience_z") |>
      dplyr::mutate(component = "zero")
  }
  
  ## Part 2: abundance given presence
  dat_pos <- dat |> dplyr::filter(abundance > 0)
  
  if (nrow(dat_pos) >= 20) {
    m2 <- tryCatch(
      lm(log(abundance) ~ cogresilience_z, data = dat_pos),
      error = function(e) NULL
    )
    
    if (!is.null(m2)) {
      out$count <- broom::tidy(m2) |>
        dplyr::filter(term == "cogresilience_z") |>
        dplyr::mutate(component = "count")
    }
  }
  
  dplyr::bind_rows(out) |> dplyr::mutate(metabolite = metab)
}

# Running loop
results <- purrr::map_dfr(metabolite_cols, fit_two_part)

#Making result tables

count_results <- results %>%
  filter(component == "count", term == "cogresilience_z")

zero_results <- results %>%
  filter(component == "zero", term == "cogresilience_z")

count_results <- count_results %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

zero_results <- zero_results %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

#selecting significant results
sig_count <- count_results %>% filter(p_adj < 0.05)
sig_zero  <- zero_results  %>% filter(p_adj < 0.05)


#Assigning OR and percentages

# Zero part: odds ratios
sig_zero_export <- sig_zero %>%
  mutate(odds_ratio = exp(estimate))

# Count part: percent change
sig_count_export <- sig_count %>%
  mutate(percent_change = (exp(estimate) - 1) * 100)

#Export
#write.csv(sig_zero_export,
#    file = "significant_metabolites_zero_part.csv",
#    row.names = FALSE)

#write.csv(sig_count_export,
#    file = "significant_metabolites_count_part.csv",
#     row.names = FALSE)

#Merging with carnitines for downstream analysis

supp_stats <- all_results %>%
  mutate(scan_number = stringr::str_remove(metabolite, "^X")) %>%
  inner_join(carnitines, by = "scan_number")

supp_stats_wide <- supp_stats %>%
  dplyr::select(
    scan_number,
    Compound_Name,
    Precursor_MZ,
    IonMode,
    Predicted_formula,
    component,
    estimate,
    std.error,
    p.value,
    p_adj
  ) %>%
  tidyr::pivot_wider(
    names_from = component,
    values_from = c(
      estimate,
      std.error,
      p.value,
      p_adj
    )
  )


#Export supplementary table 3
#write.csv(
#  supp_stats_wide,
#  "Supplementary_table_3.csv",
#  row.names = FALSE
#)

#Counting carnitine stats
sig_carnitines <- supp_stats %>%
  filter(p_adj < 0.05) %>%
  mutate(
    direction = case_when(
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "zero"
    )
  )

sig_summary_carnitines <- sig_carnitines %>%
  count(component, direction)

sig_summary_carnitines


carbon <- as.numeric(sub("^C([0-9]+).*", "\\1", supp_stats_wide$Predicted_formula))
sum(carbon >= 6 & carbon <= 12, na.rm = TRUE)
sum(carbon >= 13 & carbon <= 20, na.rm = TRUE)
sum(carbon >= 2 & carbon <= 5, na.rm = TRUE)

# Plot 4 a and c
volcano_dat <- supp_stats %>%
  mutate(
    neg_log10_fdr = -log10(p_adj),
    direction = case_when(
      estimate > 0 ~ "Positive",
      estimate < 0 ~ "Negative"
    )
  )
figure4ac<-ggplot(volcano_dat,
                 aes(x = estimate,
                     y = neg_log10_fdr,
                     color = direction)) +
  
  geom_point(alpha = 0.7, size = 2) +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  
  facet_wrap(~ component, scales = "free_x") +
  
  labs(
    x = "Effect size (model-specific)",
    y = expression(-log[10]("FDR")),
    color = "Association",
    title = "Hurdle model results for carnitines",
    subtitle = "Zero part: detection probability; Count part: abundance among detected"
  ) +
  
  scale_color_manual(values = c(
    Positive = "#583F8F",
    Negative = "#E1642F"
  )) +
  
  theme_classic(base_size = 14)

figure4ac

#ggsave("figure4ac.svg", plot = figure4ac, dpi =300, width = 10, height = 6)

#Creating plot loop per metabolite
plot_hurdle <- function(df, metabolite_col, log_scale = TRUE) {
  
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  # Prepare data
  plot_dat <- df %>%
    dplyr::select(cogresilience_z, all_of(metabolite_col)) %>%
    dplyr::rename(abundance = all_of(metabolite_col)) %>%
    dplyr::mutate(
      present = abundance > 0,
      cogresilience_z = as.numeric(cogresilience_z)  # ensure numeric
    )
  
  # Non-zero data for regression (count part)
  reg_dat <- plot_dat %>% dplyr::filter(abundance > 0)
  
  # Fit log10-linear regression on non-zero abundances
  m_count <- lm(log10(abundance) ~ cogresilience_z, data = reg_dat)
  
  # Prediction dataset
  pred_dat <- data.frame(
    cogresilience_z = seq(
      min(reg_dat$cogresilience_z, na.rm = TRUE),
      max(reg_dat$cogresilience_z, na.rm = TRUE),
      length.out = 100
    )
  )
  
  # Predict with standard errors
  pred <- predict(m_count, newdata = pred_dat, se.fit = TRUE)
  pred_dat <- pred_dat %>%
    mutate(
      fit = 10^pred$fit,
      lwr = 10^(pred$fit - 1.96 * pred$se.fit),
      upr = 10^(pred$fit + 1.96 * pred$se.fit)
    )
  
  # Optional: logistic part for probability of detection (hurdle model)
  m_zero <- glm(present ~ cogresilience_z, data = plot_dat, family = binomial)
  pred_dat$prob_present <- predict(m_zero, newdata = pred_dat, type = "response")
  
  # Base plot
  p <- ggplot(plot_dat, aes(x = cogresilience_z, y = abundance)) +
    
    # Points
    geom_point(aes(color = present), alpha = 0.6, size = 2) +
    
    # Regression line (count part)
    geom_line(data = pred_dat, aes(x = cogresilience_z, y = fit),
              color = "black", linewidth = 1) +
    
    # Confidence ribbon
    geom_ribbon(data = pred_dat, aes(x = cogresilience_z, ymin = lwr, ymax = upr),
                inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +
    
    # Labels
    labs(
      x = "Cognitive resilience (z-score)",
      y = "Metabolite abundance",
      color = "Detected",
      title = paste0("Scan ", gsub("^X", "", metabolite_col)),
      subtitle = "Count-part regression on non-zero values; red dashed = probability of detection"
    ) +
    
    scale_color_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "#D95F02")) +
    
    theme_classic(base_size = 14)
  
  # Optional log10 scale for y-axis
  if (log_scale) {
    p <- p + scale_y_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  }
  
  return(p)
}


#Plot 4h
figure4h<-plot_hurdle(df, "X79204")

figure4h
#ggsave("figure4h.svg", plot = figure4h, dpi =300, width = 10, height = 7)



