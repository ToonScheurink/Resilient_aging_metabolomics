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
library(kableExtra)

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

### Running model
X <- metab_resid_filtered
Y <- df$cogresilience

var_resid <- apply(metab_resid_filtered, 2, var)
top_idx <- order(var_resid, decreasing = TRUE)[1:5000]
X_sub <- metab_resid_filtered[, top_idx]

spls_final <- spls(X_sub, Y, ncomp = 1, keepX = 3063) 

### Model performance (figure 2A)
# Get component scores
comp1_scores <- spls_final$variates$X[,1]

# Fit linear model: CR ~ sPLS Component 1
model <- lm(CR ~ Comp1, data = data.frame(Comp1 = comp1_scores, CR = df$cogresilience))
summary_model <- summary(model)

summary_model

# Extract statistics
beta <- round(summary_model$coefficients[2, 1], 2)   # slope
p_val <- signif(summary_model$coefficients[2, 4], 3) # p-value
r2 <- round(summary_model$r.squared, 2)             # R-squared

# Set position for annotation
x_pos <- max(comp1_scores) - 0.05 * diff(range(comp1_scores))
y_pos <- min(df$cogresilience) + 0.25 * diff(range(df$cogresilience))

# Create plot with annotation
sPLS_cogres <- ggplot(
  data.frame(Comp1 = comp1_scores, CR = df$cogresilience),
  aes(x = Comp1, y = CR)
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

### Internal validation
set.seed(500)

# Define cross-validation setup 
folds <- createFolds(Y, k = 10, list = TRUE)

cv_spls <- function(X_sub, Y, folds, ncomp = 1, keepX) {
  results <- data.frame(R2 = NA, Spearman = NA, RMSE = NA)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(Y), test_idx)
    
    X_train <- X_sub[train_idx, , drop = FALSE]
    y_train <- Y[train_idx]
    X_test  <- X_sub[test_idx, , drop = FALSE]
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

cv_selected <- cv_spls(X_sub, Y, folds, ncomp = 1, keepX = rep(3063, 1))

summary(cv_selected)
mean(cv_selected$R2)
mean(cv_selected$Spearman)
mean(cv_selected$RMSE)


### Permutation testing
set.seed(500)

n_perm <- 100
perm_results <- numeric(n_perm)

for (p in 1:n_perm) {
  # Shuffle the outcome
  Y_perm <- sample(Y)
  
  # Run the same CV with permuted Y
  cv_perm <- cv_spls(X_sub, Y_perm, folds, ncomp = 1, keepX = rep(3063, 1))
  
  # Store mean R²
  perm_results[p] <- mean(cv_perm$R2)
  
  cat("Permutation", p, "done\n")
}

### Feature stability
set.seed(500)

folds <- createFolds(Y, k = 10, list = TRUE)

# Prepare matrix to record selections
var_selection <- matrix(0, nrow = ncol(X_sub), ncol = length(folds))
rownames(var_selection) <- colnames(X_sub)

keepX_value <- 50  

for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(Y), test_idx)
  
  X_train <- X_sub[train_idx, , drop = FALSE]
  y_train <- Y[train_idx]
  
  model <- spls(X_train, y_train, ncomp = 1, keepX = rep(keepX_value, 1))
  
  selected_vars <- selectVar(model, comp = 1)$X$name
  
  # Mark selected variables
  var_selection[selected_vars, i] <- 1
}

# Compute selection frequencies
var_freq <- rowSums(var_selection) / length(folds)

var_selection <- matrix(0, nrow = ncol(X_sub), ncol = length(folds))
rownames(var_selection) <- colnames(X_sub)

for (i in seq_along(folds)) {
  train_idx <- setdiff(seq_along(Y), folds[[i]])
  model <- spls(X_sub[train_idx, , drop = FALSE], Y[train_idx],
                ncomp = 1, keepX = rep(3063, 1))
  sel <- selectVar(model, comp = 1)$X$name
  if (!is.null(sel)) var_selection[sel, i] <- 1
}

# Compute frequency across folds
var_freq <- rowSums(var_selection) / length(folds)

# Create a dataframe
df_stability <- data.frame(
  Metabolite = names(var_freq),
  Stability = var_freq
)

### Extracting top 3063, rank and weight
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

#Adding stability
df_loadings_full <- df_loadings %>%
  dplyr::left_join(df_stability, by = "Metabolite") %>%
  arrange(desc(abs(Weight)))

#saveRDS(df_loadings_full, file = "ranked_metabolites.rds")

#### External validation ####
### External data prep
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

#Selction of metabolites only
df_metab_external <- df_external %>%
  column_to_rownames(var = "filename") %>%  
  select(3:22663) %>%                
  filter(rowSums(.) > 0)           

#RCLR conversion
df_rclr_external <- decostand(df_metab_external, method = "rclr")
df_rclr_matrix_external <- as.matrix(df_rclr_external)

#Add age as covariate
covariates_external <- df_external %>% select(host_age)
metab_resid_external <- apply(df_rclr_matrix_external, 2, function(x) {
  fit_external <- lm(x ~ ., data = covariates_external)
  resid(fit_external)
})

#Removing zero variance metabolites
var_resid_external <- apply(metab_resid_external, 2, var)
metab_resid_filtered_external <- metab_resid_external[, var_resid_external > 1e-6]

### Checking missing features
train_features <- colnames(X_sub)   # these are the top 5000 metabolites used in training

# Align columns between training and external data
common_features <- intersect(train_features, colnames(metab_resid_filtered_external))

# Subset the external matrix to those columns 
X_ext <- metab_resid_filtered_external[, common_features, drop = FALSE]

length(common_features)                  # how many overlap
length(setdiff(train_features, common_features))  # missing metabolites

# Get all training features
train_features <- colnames(X_sub)

# Create a copy of the external matrix
X_ext_aligned <- X_ext

# Identify missing metabolites
missing_features <- setdiff(train_features, colnames(X_ext_aligned))

# Convert to data frame
X_ext_aligned <- as.data.frame(X_ext)

# Add missing features as zeros
for (feat in missing_features) {
  X_ext_aligned[[feat]] <- 0
}

# Reorder columns to match training order
X_ext_aligned <- X_ext_aligned[, train_features]

# Convert back to matrix for prediction
X_ext_aligned <- as.matrix(X_ext_aligned)

# Reorder columns to match training set order exactly
X_ext_aligned <- X_ext_aligned[, train_features]

# Double-check
all(colnames(X_ext_aligned) == colnames(X_sub))  # should return TRUE

### External performance
preds_ext <- predict(spls_final, newdata = X_ext_aligned)

y_pred <- preds_ext$predict[, 1, 1]

df_pred <- data.frame(
  Observed = Y,
  Predicted = y_pred
)

# Fit linear model: CR ~ sPLS Component 1
model <- lm(Predicted ~ Observed, data = df_pred)
summary_model <- summary(model)

cor.test(df_pred$Observed, df_pred$Predicted, method = "pearson")  # correlation
summary(lm(Predicted ~ Observed, data = df_pred))    # R²
sqrt(mean((df_pred$Observed - df_pred$Predicted)^2))

### Permutation testing external model
# --- Step 1: Extract observed values ---
y_true <- df_external$cogresilience
y_pred <- preds_ext$predict[, 1, 1]

# Observed correlation
obs_cor <- cor(y_true, y_pred, method = "pearson")

# --- Step 2: Permutation test ---
set.seed(123)          # for reproducibility
n_perm <- 100         # number of permutations
perm_cor <- numeric(n_perm)

for (i in 1:n_perm) {
  y_perm <- sample(y_true)  # shuffle outcomes
  # predict with the same model (model is fixed, no retraining)
  perm_cor[i] <- cor(y_perm, y_pred, method = "pearson")
}

# --- Step 3: Empirical p-value ---
p_value <- mean(abs(perm_cor) >= abs(obs_cor))

# --- Step 4: Summary and visualization ---
cat("Observed correlation:", round(obs_cor, 3), "\n")
cat("Empirical p-value from permutation:", p_value, "\n")

# Convert permutation correlations to a data frame
perm_cor_df <- data.frame(Perm_Cor = perm_cor)

# Observed correlation
obs_cor_val <- obs_cor

# Create ggplot
perm_cor_model <- ggplot(perm_cor_df, aes(x = Perm_Cor)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#B0BEC5", color = "white") +
  geom_density(fill = "#90A4AE", alpha = 0.3) +  # smooth density overlay
  geom_vline(xintercept = obs_cor_val, color = "#E53935", linetype = "dashed", size = 1.5) +
  annotate("text",
           x = obs_cor_val,
           y = max(density(perm_cor)$y)*0.4,
           label = paste0("Observed r = ", round(obs_cor_val, 3)),
           color = "#E53935",
           angle = 90,
           vjust = -0.5,
           hjust = 0,
           size = 5) +
  labs(
    x = "Correlation with predicted values",
    y = "Density"  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30")
  )

# Display
perm_cor_model

### External stability of features
# Running the same model parameters on the external dataset
train_selected <- selectVar(spls_final, comp = 1)$X
spls_ext <- spls(X_ext_aligned, y_true, ncomp = 1, keepX = 3063)
ext_selected <- selectVar(spls_ext, comp = 1)$X

# Extract corresponding weights
train_features <- train_selected$name
train_weights  <- train_selected$value$value.var

ext_features <- ext_selected$name
ext_weights  <- ext_selected$value$value.var

common_features <- intersect(train_features, ext_features)

train_weights_common <- train_weights[match(common_features, train_features)]
ext_weights_common   <- ext_weights[match(common_features, ext_features)]

# Feature overlap (Jaccard index)
jaccard <- length(common_features) / length(union(train_features, ext_features))

# Correlation of loadings
loading_cor <- cor(train_weights_common, ext_weights_common)

cat("Feature overlap (Jaccard):", round(jaccard, 3), "\n")
cat("Loading correlation:", round(loading_cor, 3), "\n")

### Stability of high rank features
quant <- 0.95  # change to change selected percentile (eg 0.8 = top 20%)

thr_train <- quantile(abs(train_weights), probs = quant)

train_high <- train_features[abs(train_weights) >= thr_train]

# Presence in external set
present_in_external <- train_high %in% ext_features

# Count
num_present <- sum(present_in_external)
num_total   <- length(train_high)
prop_present <- num_present / num_total

cat("High-weight training features present in external set:", 
    num_present, "of", num_total, 
    sprintf("(%.3f)", prop_present), "\n")

### Plotting feature stability (Figure Suplementary X)
# Ensure weights are numeric vectors
train_weights_num <- train_selected$value$value.var
ext_weights_num   <- ext_selected$value$value.var

# Combine
weights_df <- data.frame(
  feature = train_features,
  train_weight = train_weights
)

weights_df$ext_weight <- NA

# Find indices of external features that exist in training
ext_idx <- match(ext_features, weights_df$feature)

# Keep only non-NA indices
valid_idx <- !is.na(ext_idx)
ext_idx <- ext_idx[valid_idx]
ext_weights_valid <- ext_weights[valid_idx]

# Assign the weights
weights_df$ext_weight[ext_idx] <- ext_weights_valid

# Prepare data
weights_df <- weights_df %>%
  mutate(
    missing_in_external = ifelse(is.na(ext_weight), "Missing in external", "Present in external"),
    abs_train = abs(train_weight)
  )

# Plot
feature_stability_ext<-ggplot(weights_df, aes(x = abs_train, fill = missing_in_external)) +
  geom_histogram(
    aes(y = ..density..),
    bins = 30,                  # small bins for smoother appearance
    position = "stack",         # stacked histogram
    color = "black", size = 0.3
  ) +
  # Density curve for missing group only
  geom_density(
    data = subset(weights_df, missing_in_external == "Missing in external"),
    aes(y = ..density..),
    color = "black",
    fill = NA,
    size = 1.3
  ) +
  scale_fill_manual(
    values = c(
      "Present in external" = "#1b9e77",
      "Missing in external" = "#d95f02"
    )
  ) +
  labs(
    x = "Absolute training loadings",
    y = "Density",
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top"
  )

feature_stability_ext

