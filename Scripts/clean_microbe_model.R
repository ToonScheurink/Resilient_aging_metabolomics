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
library(miaViz)
library(ggrepel)
library(gridExtra)

### Loading data
df1 <- readRDS("clean_data_wellcome.rds")

df1$host_age<-as.numeric(df1$host_age)
df1$host_body_mass_index<-as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl<-as.numeric(df1$physresilience_pctl)
df1$physresilience<-as.numeric(df1$physresilience)
df1$cogresilience<-as.numeric(df1$cogresilience)
df1$cogresilience_pctl<-as.numeric(df1$cogresilience_pctl)

colnames(df1) <- as.character(colnames(df1))

microbe_scans <- fread("unique_microbial_scans.csv")[[1]]
microbe_scans <- as.character(microbe_scans)

#Checking missings
missing_scans <- setdiff(microbe_scans, colnames(df))
missing_scans

#Age harmonization
target_age <- 73.7

df <- df1 %>%
  group_by(record_ID) %>%
  slice_min(order_by = abs(host_age - target_age), with_ties = FALSE) %>%
  ungroup()

### Normalization and age correction
# Selection of metabolites only (based on microbe_scans list)
df_metab <- df %>%
  column_to_rownames(var = "filename") %>%
  select(all_of(microbe_scans)) %>%   
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

### Running model
X <- metab_resid_filtered
Y <- df$cogresilience

spls_final <- spls(X, Y, ncomp = 1, keepX = ncol(X))  


### Plotting model performance (not included in manuscript)
# Get component scores
comp1_scores <- spls_final$variates$X[,1]

# Fit linear model: CR ~ sPLS Component 1
model <- lm(CR ~ Comp1, data = data.frame(Comp1 = comp1_scores, CR = df$cogresilience))
summary_model <- summary(model)

# Extract statistics
beta <- round(summary_model$coefficients[2, 1], 2)   # slope
p_val <- signif(summary_model$coefficients[2, 4], 3) # p-value
r2 <- round(summary_model$r.squared, 2)             # R-squared

# Set position for annotation
x_pos <- max(comp1_scores) - 0.05 * diff(range(comp1_scores))
y_pos <- min(df$cogresilience) + 0.15 * diff(range(df$cogresilience))

# Create plot with annotation
sPLS_cogres<-ggplot(data.frame(Comp1 = comp1_scores, CR = df$cogresilience),
                    aes(x = Comp1, y = CR)) +
  geom_point(alpha = 0.8, color = "#2C3E50", size = 3) +
  geom_smooth(method = "lm", color = "#0072B2", fill = "#0072B2", alpha = 0.2, se = TRUE, size = 1.2) +
  labs(
    x = "sPLS Component 1 Score",
    y = "Cognitive Resilience"
  ) +
  annotate(
    "label",
    x = x_pos, y = y_pos,
    label = paste0("R² = ", round(r2, 2),
                   "\nβ = ", round(beta, 2),
                   "\np = ", format.pval(p_val, digits = 2, eps = 0.001)),
    hjust = 1, vjust = 1,
    size = 5,
    label.size = 0.3,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.margin = margin(15, 15, 15, 15)
  )

sPLS_cogres

### Internal validation
set.seed(500)

# Define cross-validation setup (e.g. 10-fold CV)
folds <- createFolds(Y, k = 10, list = TRUE)

cv_spls <- function(X, Y, folds, ncomp = 1, keepX) {
  results <- data.frame(R2 = NA, Spearman = NA, RMSE = NA)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(Y), test_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- Y[train_idx]
    X_test  <- X[test_idx, , drop = FALSE]
    y_test  <- Y[test_idx]
    
    # Fit sPLS on training set
    model <- spls(X_train, y_train, ncomp = ncomp, keepX = keepX)
    
    # Predict Y on test samples
    preds <- predict(model, X_test)
    y_pred <- preds$predict[, 1, 1]   # first Y, first component
    
    # Evaluate metrics in test data
    lm_fit <- lm(y_test ~ y_pred)
    R2  <- summary(lm_fit)$r.squared
    rho <- suppressWarnings(cor(y_test, y_pred, method = "spearman"))
    RMSE <- sqrt(mean((y_test - y_pred)^2))
    
    results[i, ] <- c(R2, rho, RMSE)
  }
  
  return(results)
}

cv_selected <- cv_spls(X, Y, folds, ncomp = 1, keepX = rep(742, 1))

summary(cv_selected)
mean(cv_selected$R2)
mean(cv_selected$Spearman)
mean(cv_selected$RMSE)

### Permutation checks
set.seed(500)

n_perm <- 100
perm_results <- numeric(n_perm)

for (p in 1:n_perm) {
  # Shuffle the outcome
  Y_perm <- sample(Y)
  
  # Run the same CV with permuted Y
  cv_perm <- cv_spls(X, Y_perm, folds, ncomp = 1, keepX = rep(742, 1))
  
  # Store mean R²
  perm_results[p] <- mean(cv_perm$R2)
  
  cat("Permutation", p, "done\n")
}

# Compare distributions
real_R2 <- mean(cv_selected$R2)

### Feature stability
set.seed(500)

folds <- createFolds(Y, k = 10, list = TRUE)

# Prepare matrix to record selections
var_selection <- matrix(0, nrow = ncol(X), ncol = length(folds))
rownames(var_selection) <- colnames(X)

keepX_value <- 50 

for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(Y), test_idx)
  
  X_train <- X[train_idx, , drop = FALSE]
  y_train <- Y[train_idx]
  
  model <- spls(X_train, y_train, ncomp = 1, keepX = rep(keepX_value, 1))
  
  selected_vars <- selectVar(model, comp = 1)$X$name
  
  # Mark selected variables
  var_selection[selected_vars, i] <- 1
}

# Compute selection frequencies
var_freq <- rowSums(var_selection) / length(folds)

##### External validation #####
###Prepping data
# Step 1: Get the visit closest to the target age
df_selected <- df1 %>%
  group_by(record_ID) %>%
  slice_min(order_by = abs(host_age - target_age), with_ties = FALSE) %>%
  ungroup()

# Step 2: Join that selection back to the full dataset
df_external <- df1 %>%
  dplyr::inner_join(df_selected %>% select(record_ID, selected_age = host_age),
                    by = "record_ID") %>%
  # Step 3: Exclude the originally selected visit
  filter(host_age != selected_age) %>%
  # Step 4: Pick the visit closest in age to the selected one
  group_by(record_ID) %>%
  slice_min(order_by = abs(host_age - selected_age), with_ties = FALSE) %>%
  ungroup()

### Normalization and age correction
#Selction of metabolites only
df_metab_external <- df %>%
  column_to_rownames(var = "filename") %>%
  select(all_of(microbe_scans)) %>%   
  filter(rowSums(.) > 0)

#RCLR conversion
df_rclr_external <- decostand(df_metab_external, method = "rclr")
df_rclr_matrix_external <- as.matrix(df_rclr_external)

#Add age as covariate(change covariates as needed)
covariates_external <- df_external %>% select(host_age)
metab_resid_external <- apply(df_rclr_matrix_external, 2, function(x) {
  fit_external <- lm(x ~ ., data = covariates_external)
  resid(fit_external)
})

#Removing zero variance metabolites
var_resid_external <- apply(metab_resid_external, 2, var)
metab_resid_filtered_external <- metab_resid_external[, var_resid_external > 1e-6]

#finding overlap
train_features <- colnames(X)   # these are the metabolites used in training

# Align columns between training and external data
common_features <- intersect(train_features, colnames(metab_resid_filtered_external))

# Subset the external matrix to those columns (in the same order!)
X_ext <- metab_resid_filtered_external[, common_features, drop = FALSE]

length(common_features)                  # how many overlap
length(setdiff(train_features, common_features))  # missing metabolites

### Perfomance model
preds_ext <- predict(spls_final, newdata = X_ext)

y_pred <- preds_ext$predict[, 1, 1]

df_pred <- data.frame(
  Observed = Y,
  Predicted = y_pred
)

# Fit linear model: CR ~ sPLS Component 1
model <- lm(Predicted ~ Observed, data = df_pred)
summary_model <- summary(model)
summary_model

# Extract statistics
beta <- round(summary_model$coefficients[2, 1], 2)   # slope
p_val <- signif(summary_model$coefficients[2, 4], 3) # p-value
r2 <- round(summary_model$r.squared, 2)             # R-squared

# Set position for annotation
x_pos <- max(df_pred$Predicted) - 2.4 * diff(range(df_pred$Predicted))
y_pos <- min(df_pred$Observed) + 0.65 * diff(range(df_pred$Observed))


### Plot (figure not included in manuscript)
external_validation<-ggplot(df_pred, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    x = "Observed cognitive resilience (external dataset)",
    y = "Predicted cognitive resilience"
  ) +
  annotate(
    "label",
    x = x_pos, y = y_pos,
    label = paste0("R² = ", round(r2, 2),
                   "\nβ = ", round(beta, 2),
                   "\np = ", format.pval(p_val, digits = 2, eps = 0.001)),
    hjust = 1, vjust = 1,
    size = 5,
    label.size = 0.3,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.margin = margin(15, 15, 15, 15)
  )

external_validation



