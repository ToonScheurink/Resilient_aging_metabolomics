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
library(inflection) 

#Loading data
df <- readRDS("clean_data_wellcome.rds")
df$host_age<-as.numeric(df$host_age)
df$host_body_mass_index<-as.numeric(df$host_body_mass_index)
df$physresilience_pctl<-as.numeric(df$physresilience_pctl)
df$physresilience<-as.numeric(df$physresilience)
df$cogresilience<-as.numeric(df$cogresilience)
df$cogresilience_pctl<-as.numeric(df$cogresilience_pctl)

###Harmonizing age
find_best_alignment <- function(df, target_grid = seq(50, 100, by = 0.1), verbose = TRUE) {
  
  best_sd <- Inf
  best_target <- NA
  best_df <- NULL
  
  for (t in target_grid) {
    chosen <- df %>%
      group_by(record_ID) %>%
      slice_min(order_by = abs(host_age - t), with_ties = FALSE) %>%
      ungroup()
    
    current_sd <- sd(chosen$host_age)
    
    if (current_sd < best_sd) {
      best_sd <- current_sd
      best_target <- t
      best_df <- chosen
      if (verbose) {
        message("New best target: ", round(t, 2), " with SD = ", round(best_sd, 3))
      }
    }
  }
  
  list(
    aligned_df = best_df,
    best_target = best_target,
    best_sd = best_sd
  )
}

result <- find_best_alignment(df)

#Setting target age and selecting visits in df
target_age <- 73.7

df <- df %>%
  group_by(record_ID) %>%
  slice_min(order_by = abs(host_age - target_age), with_ties = FALSE) %>%
  ungroup()

### Normalization and age correction
#Selction of metabolites only
df_metab <- df %>%
  column_to_rownames(var = "filename") %>%  
  select(3:22663) %>%                
  filter(rowSums(.) > 0)           

#RCLR conversion
df_rclr <- decostand(df_metab, method = "rclr")
df_rclr_matrix <- as.matrix(df_rclr)

#Add age as covariate
covariates <- df %>% select(host_age)
metab_resid <- apply(df_rclr_matrix, 2, function(x) {
  fit <- lm(x ~ ., data = covariates)
  resid(fit)
})

#Removing zero variance metabolites
var_resid <- apply(metab_resid, 2, var)
metab_resid_filtered <- metab_resid[, var_resid > 1e-6]

### Model tuning
X <- metab_resid_filtered
Y <- df$cogresilience

var_resid <- apply(metab_resid_filtered, 2, var)
top_idx <- order(var_resid, decreasing = TRUE)[1:5000]
X_sub <- metab_resid_filtered[, top_idx]

set.seed(42)
tune <- tune.spls(
  X = X_sub,
  Y = Y,
  ncomp = 5,
  test.keepX = c(50, 100, 200, 500),
  measure = "MSE",        # <--- specify explicitly
  validation = "Mfold",
  folds = 5,
  nrepeat = 5,
  progressBar = TRUE
)

tune$choice.ncomp

#Running model for later inflection
spls_final <- spls(X_sub, Y, ncomp = 1, keepX = 5000) 

# Extract absolute weights from sPLS model
loadings <- abs(spls_final$loadings$X[,1])

# Rank them
df_ranked <- data.frame(
  Metabolite = names(loadings),
  Weight = loadings
) %>%
  arrange(desc(Weight)) %>%
  mutate(Rank = row_number())


x <- 1:length(df_ranked$Weight)
y <- df_ranked$Weight

# Detect the elbow
elbows <- inflection::findiplist(x, y, index = 1)
elbow <- elbows[length(elbows)]

# Get cutoff value
cutoff <- y[elbow]

# Number of metabolites above cutoff
selected <- sum(y > cutoff)

cat("Elbow at metabolite rank:", elbow, "\n")
cat("Cutoff weight:", cutoff, "\n")
cat("Number of metabolites above cutoff:", selected, "\n")


