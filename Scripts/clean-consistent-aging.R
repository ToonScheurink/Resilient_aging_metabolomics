#Packages
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

#Loading data
df1 <- readRDS("clean_data_wellcome.rds")

metabs <- readRDS("ranked_metabolites.rds")

df1$host_age<-as.numeric(df1$host_age)
df1$host_body_mass_index<-as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl<-as.numeric(df1$physresilience_pctl)
df1$physresilience<-as.numeric(df1$physresilience)
df1$cogresilience<-as.numeric(df1$cogresilience)
df1$cogresilience_pctl<-as.numeric(df1$cogresilience_pctl)

#Selecting metabolites
df_metab <- df1 %>%
  column_to_rownames(var = "filename") %>%  
  dplyr::select(3:22663) %>%                
  filter(rowSums(.) > 0)       


#Selecting prevelance
prev <- colMeans(df_metab > 0, na.rm = TRUE)

hist(
  prev,
  breaks = 50,
  col = "grey",
  border = "white",
  main = "Metabolite prevalence (non-zero proportion)",
  xlab = "Proportion of samples with non-zero values"
)

plot(
  sort(prev),
  seq_along(prev) / length(prev),
  type = "l",
  xlab = "Prevalence threshold",
  ylab = "Proportion of metabolites retained",
  main = "Cumulative metabolite prevalence"
) +
  abline(v = c(0.5, 0.6, 0.7), col = c("blue", "red", "darkgreen"), lty = 2)

sum(prev >= 0.7)
sum(prev >= 0.6)
sum(prev >= 0.5)
sum(prev >= 0.4)
sum(prev >= 0.3)
sum(prev >= 0.2)
sum(prev >= 0.1)
sum(prev >= 0.01)

#Cut-off
keep_metabs <- names(prev[prev >= 0.7])
df_filt <- df_metab[, keep_metabs]

#Pseudocount for downstream stats
pseudocount <- apply(df_filt, 2, function(x) min(x[x > 0], na.rm = TRUE) / 2)
df_pc <- sweep(df_filt, 2, pseudocount, "+")

#Log transform
df_log <- log2(df_pc)

#Prepping metadata and calculating within-person age
meta <- df1 %>%
  dplyr::select(filename, record_ID, visit, host_age, sex)

##Changing variable names to fit LMM package
# Save original names
original_names <- colnames(df_log)

# Make valid R names
colnames(df_log) <- make.names(colnames(df_log))

# Create a mapping
name_map <- data.frame(original = original_names, valid = colnames(df_log))

#Combining with metadata
df_lmm <- cbind(meta[, c("record_ID", "visit", "host_age")], df_log)
df_lmm <- as.data.frame(df_lmm)

#centering age
df_lmm$age_c <- scale(df_lmm$host_age, center = TRUE, scale = FALSE)

#Running LMM Loop
results <- lapply(colnames(df_log), function(metab){
  formula <- as.formula(paste(metab, "~ age_c + (1 + age_c || record_ID)"))
  fit <- lmer(
    formula,
    data = df_lmm,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  tidy(fit) %>%
    filter(term == "age_c") %>%
    mutate(
      metabolite = metab,
      singular = isSingular(fit, tol = 1e-5)
    )
})s

results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p.value, method = "BH")


#Selecting metabolites
consistent_metabs <- results_df %>%
  filter(
    p_adj < 0.05,
    singular == TRUE
  )


consistent_metabs <- consistent_metabs %>%
  mutate(
    std_slope = estimate / sqrt(var(df_lmm$age_c))
  ) %>%
  arrange(desc(abs(std_slope)))


# Number of metabolites passing all criteria
n_final <- nrow(consistent_metabs)
cat("Number of QCed age-associated metabolites:", n_final, "\n")


# Calculate the proportion of positive slopes
positive_proportion <- mean(consistent_metabs$estimate > 0)

# Convert to percentage
positive_percentage <- positive_proportion * 100

positive_percentage

#Merging with cognitive resilience model
consistent_metabs_clean <- consistent_metabs %>%
  mutate(metabolite = sub("^X", "", metabolite))

consistent_metabs_clean <- consistent_metabs_clean %>%
  rename(Metabolite = metabolite)

df_merged <- consistent_metabs_clean %>%
  left_join(metabs, by = "Metabolite")

table<-df_merged %>%
  dplyr::select(Metabolite, Weight, Direction, Stability, std_slope, estimate, p_adj)


#write.csv(table, "supl_table_1.csv", row.names = FALSE)



