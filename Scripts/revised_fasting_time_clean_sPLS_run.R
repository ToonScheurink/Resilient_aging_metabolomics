### Packages
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
library(reshape2)
library(ppcor)
library(lme4)
library(lmerTest)       
library(broom.mixed)   
library(future.apply) 
library(patchwork)
library(broom)
library(forcats)
library(flextable)
library(gtsummary)
library(gt)      
library(zCompositions)
library(compositions) 
library(pls)
library(pheatmap)
library(MASS)
library(viridis)
library(ggvegan)
library(ggfortify)
library(ggrepel)
library(gridExtra)
library(kableExtra)
library(lubridate)

### Loading data
df1 <- readRDS("clean_data_wellcome.rds")
df1$host_age<-as.numeric(df1$host_age)
df1$host_body_mass_index<-as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl<-as.numeric(df1$physresilience_pctl)
df1$physresilience<-as.numeric(df1$physresilience)
df1$cogresilience<-as.numeric(df1$cogresilience)
df1$cogresilience_pctl<-as.numeric(df1$cogresilience_pctl)

# Selecting visits
target_age <- 73.7
df <- df1 %>%
  group_by(record_ID) %>%
  slice_min(order_by = abs(host_age - target_age), with_ties = FALSE) %>%
  ungroup()

#Fasting time
df$fasting_hrs_num <- as.numeric(
  na_if(df$fasting_hrs, "")
)

cor.test(
  df$fasting_hrs_num,
  df$host_age,
  use = "complete.obs",
  method = "pearson"
)

cor.test(
  df$fasting_hrs_num,
  df$cogresilience,
  use = "complete.obs",
  method = "pearson"
)

sum(complete.cases(
  df$fasting_hrs_num,
  df$host_age,
  df$cogresilience
))

### Normalization and age correction
#Selction of metabolites only
df_metab <- df %>%
  column_to_rownames(var = "filename") %>%  
  dplyr::select(3:22663) %>%                
  filter(rowSums(.) > 0)           

#Counting metabolites
n_present_metabolites <- sum(colSums(df_metab != 0) > 0)
n_present_metabolites

#RCLR conversion
df_rclr <- decostand(df_metab, method = "rclr")
df_rclr_matrix <- as.matrix(df_rclr)

# Identify participants with complete age + fasting time
keep <- complete.cases(df$host_age, df$fasting_hrs_num)

# Restrict metabolomics and outcome to those participants
X_fasting <- df_rclr_matrix[keep, , drop = FALSE]
Y_fasting <- df$cogresilience[keep]

# Covariates for the same participants
covariates <- df[keep, ] %>%
  dplyr::select(host_age, fasting_hrs_num)

# Residualize metabolites for age + fasting time
metab_resid_fasting <- apply(X_fasting, 2, function(x) {
  fit <- lm(x ~ host_age + fasting_hrs_num, data = covariates)
  resid(fit)
})

#Removing zero variance metabolites
var_resid <- apply(metab_resid_fasting, 2, var)
metab_resid_filtered <- metab_resid_fasting[, var_resid > 1e-6]

### Running model
X <- metab_resid_filtered
Y <- Y_fasting

var_resid <- apply(metab_resid_filtered, 2, var)
top_idx <- order(var_resid, decreasing = TRUE)[1:5000]
X_sub <- metab_resid_filtered[, top_idx]

spls_final <- spls(X_sub, Y, ncomp = 1, keepX = 3080) 

### Model performance (figure 2A)
# Get component scores
comp1_scores <- spls_final$variates$X[,1]

# Fit linear model: CR ~ sPLS Component 1
model <- lm(CR ~ Comp1, data = data.frame(Comp1 = comp1_scores, CR = df$cogresilience[keep]))
summary_model <- summary(model)

summary_model

# Extract statistics
beta <- round(summary_model$coefficients[2, 1], 2)   # slope
p_val <- signif(summary_model$coefficients[2, 4], 3) # p-value
r2 <- round(summary_model$r.squared, 2)             # R-squared

# Set position for annotation
x_pos <- max(comp1_scores) - 0.05 * diff(range(comp1_scores))
y_pos <- min(df$cogresilience[keep]) + 0.25 * diff(range(df$cogresilience[keep]))

# Create plot with annotation
sPLS_cogres <- ggplot(
  data.frame(Comp1 = comp1_scores, CR = df$cogresilience[keep]),
  aes(x = Comp1, y = CR, color = df$visit)
) +
  geom_point(alpha = 0.8, color = "#2C3E50", size = 3) +
  geom_smooth(method = "lm", color = "#0072B2", fill = "#0072B2",
              alpha = 0.2, se = TRUE, size = 1.2) +
  labs(
    x = "sPLS Component 1 Score",
    y = NULL
  ) +
  annotate(
    "label",
    x = x_pos, y = y_pos,
    label = paste0("R² = ", round(r2, 2),
                   "\nβ = ", round(beta, 2),
                   "\np = ", format.pval(p_val, 2)),
    hjust = 1, vjust = 1,
    size = 8,
    label.size = 0.3,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),       
    axis.ticks.length = unit(3, "pt"),               
    plot.margin = margin(15, 15, 15, 15)
  )


sPLS_cogres

#ggsave("sPLS_cogres.svg", plot = sPLS_cogres, dpi =300, width = 10, height = 7)


### Extracting top 3080, rank and weight
# Loadings for the first component
loadings <- spls_final$loadings$X[, 1]

# Keep only nonzero loadings
nonzero_idx <- loadings != 0

# Creating a df for extraction
df_loadings <- data.frame(
  Metabolite = names(loadings)[nonzero_idx],
  Weight = loadings[nonzero_idx]
)

#Adding direction to the df
df_loadings$Direction <- ifelse(df_loadings$Weight > 0, "Positive", "Negative")

# Sort df by absolute weight and add rank
df_loadings <- df_loadings %>%
  arrange(desc(abs(Weight))) %>%
  mutate(Rank = row_number())

#uncorrected ranking
metabs <- readRDS("ranked_metabolites.rds")

#saveRDS(df_loadings, file = "ranked_metabolites_storage_time_corrected.rds")

#Overlap and weight correlation
# Jaccard similarity
length(intersect(df_loadings$Metabolite, metabs$Metabolite)) /
  length(union(df_loadings$Metabolite, metabs$Metabolite))

# Keep metabolite IDs and weights, and give each signature distinct names
new_weights <- df_loadings %>%
  select(Metabolite, Weight) %>%
  rename(Weight_new = Weight)

old_weights <- metabs %>%
  select(Metabolite, Weight) %>%
  rename(Weight_old = Weight)

# Match metabolites present in both signatures
matched_weights <- inner_join(
  new_weights,
  old_weights,
  by = "Metabolite"
)

# How many metabolites overlap?
nrow(matched_weights)

# Correlate weights
cor.test(
  matched_weights$Weight_new,
  matched_weights$Weight_old,
  method = "pearson"
)

cor.test(
  matched_weights$Weight_new,
  matched_weights$Weight_old,
  method = "spearman"
)

#retainment of top 10%
n_top10_old <- ceiling(0.10 * nrow(metabs))

top10_old <- metabs %>%
  arrange(desc(abs(Weight))) %>%
  slice_head(n = n_top10_old)

top10_retained <- top10_old %>%
  inner_join(
    df_loadings %>%
      select(Metabolite, Weight_new = Weight),
    by = "Metabolite"
  )

nrow(top10_retained)
nrow(top10_retained) / nrow(top10_old)

#retainment of weight
cor.test(
  top10_retained$Weight,
  top10_retained$Weight_new,
  method = "pearson"
)

cor.test(
  top10_retained$Weight,
  top10_retained$Weight_new,
  method = "spearman"
)

r_all <- cor(
  matched_weights$Weight_old,
  matched_weights$Weight_new,
  method = "pearson"
)

r_all

#Plot
figS3a<-ggplot(matched_weights,
       aes(x = Weight_old, y = Weight_new)) +
  geom_point(alpha = 0.4, size = 1.8) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Pearson r = ", round(r_all, 2)),
    hjust = 1.1, vjust = 1.5,
    size = 5
  ) +
  labs(
    x = "Original loading weight",
    y = "Fasting time-adjusted loading weight"
  ) +
  theme_classic(base_size = 13)

figS3a

#ggsave("figS3a.svg", plot = figS3a, dpi =300, width = 10, height = 7)


#Cross-referencing down-stream highlighted metabolites
metabs_of_interest <- c(
  "82496",  # S1P
  "79204",  # DCA-GlcA
  "77192",  # Carnitine C12:1
  "78074",  # Carnitine C9:0
  "82799",  # Piperine
  "117986"  # Lutein
)

interest<-matched_weights %>%
  filter(Metabolite %in% metabs_of_interest) %>%
  select(Metabolite, Weight_new, Weight_old)


#Carnitines
carnitines <-read.csv("supl_table_2.csv")
carnitines <- carnitines %>%
  filter(!is.na(estimate))

carnitines<-carnitines %>% mutate(Metabolite = as.character(Metabolite))

#Removing non-significant carnitines from df_loadings
df_loadings_filtered <- df_loadings %>%
  semi_join(carnitines, by = "Metabolite")

#Overlap and weight correlation
# Jaccard similarity
length(intersect(df_loadings_filtered$Metabolite, carnitines$Metabolite)) /
  length(union(df_loadings_filtered$Metabolite, carnitines$Metabolite))

# Keep metabolite IDs and weights, and give each signature distinct names
new_weights_filtered <- df_loadings_filtered %>%
  select(Metabolite, Weight) %>%
  rename(Weight_new = Weight)

old_weights_carnitines <- carnitines %>%
  select(Metabolite, Weight) %>%
  rename(Weight_old = Weight)

# Match metabolites present in both signatures
matched_weights_carnitines <- inner_join(
  new_weights_filtered,
  old_weights_carnitines,
  by = "Metabolite")

# Correlate weights
cor.test(
  matched_weights_carnitines$Weight_new,
  matched_weights_carnitines$Weight_old,
  method = "pearson"
)

cor.test(
  matched_weights_carnitines$Weight_new,
  matched_weights_carnitines$Weight_old,
  method = "spearman"
)

r_all_carnitines <- cor(
  matched_weights_carnitines$Weight_old,
  matched_weights_carnitines$Weight_new,
  method = "pearson"
)

  
#Plot
ggplot(matched_weights_carnitines,
       aes(x = Weight_old, y = Weight_new)) +
  geom_point(alpha = 0.4, size = 1.8) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Pearson r = ", round(r_all_carnitines, 2)),
    hjust = 1.1, vjust = 1.5,
    size = 5
  ) +
  labs(
    x = "Original loading weight",
    y = "Storage-adjusted loading weight",
    title = "Comparison of metabolite loading weights"
  ) +
  theme_classic(base_size = 13)

