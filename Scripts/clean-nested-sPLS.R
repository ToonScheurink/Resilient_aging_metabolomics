############################################################
# NESTED CROSS-VALIDATION FOR sPLS
#
# Reviewer-requested analysis:
# - Outer 10-fold CV estimates predictive performance
# - All preprocessing / feature selection is performed
#   using training data only
# - Inner 5-fold CV tunes keepX
# - ncomp = 1, matching the final model
# - Performance is calculated from pooled outer
#   out-of-fold predictions
############################################################


### Packages ------------------------------------------------

library(data.table)
library(dplyr)
library(tibble)
library(caret)
library(mixOmics)
library(tidyverse)
library(vegan)


### 1. Load and prepare data --------------------------------

df1 <- readRDS("clean_data_wellcome.rds")

df1$host_age <- as.numeric(df1$host_age)
df1$host_body_mass_index <- as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl <- as.numeric(df1$physresilience_pctl)
df1$physresilience <- as.numeric(df1$physresilience)
df1$cogresilience <- as.numeric(df1$cogresilience)
df1$cogresilience_pctl <- as.numeric(df1$cogresilience_pctl)


### 2. Select the primary visit -----------------------------

target_age <- 73.7

df <- df1 %>%
  group_by(record_ID) %>%
  slice_min(
    order_by = abs(host_age - target_age),
    with_ties = FALSE
  ) %>%
  ungroup()


### 3. Extract metabolomics data ----------------------------

# Keep exactly the same metabolite columns as your
# original analysis.

df_metab <- df %>%
  column_to_rownames(var = "filename") %>%
  dplyr::select(3:22663)

# Remove samples with no detected metabolites.
# This follows your original preprocessing.
keep_samples <- rowSums(df_metab) > 0

df_metab <- df_metab[keep_samples, , drop = FALSE]

# Keep outcome and age aligned with retained samples
Y <- df$cogresilience[keep_samples]
age <- df$host_age[keep_samples]

# Sanity check
stopifnot(nrow(df_metab) == length(Y))
stopifnot(nrow(df_metab) == length(age))


### 4. RCLR transformation -------------------------------

# RCLR is sample-wise, so it does not estimate parameters
# from the outcome and therefore can be performed separately
# within each training/test set.

df_rclr <- decostand(
  df_metab,
  method = "rclr"
)

X_rclr <- as.matrix(df_rclr)

cat(
  "Number of participants:",
  nrow(X_rclr),
  "\n"
)

cat(
  "Number of metabolites:",
  ncol(X_rclr),
  "\n"
)


############################################################
# 5. Helper functions
############################################################

# ----------------------------------------------------------
# Age correction
# ----------------------------------------------------------
#
# This reproduces lm(x ~ host_age) but calculates the
# coefficients matrix-wise, which is much faster when there
# are ~20,000 metabolites.
#
# The important point is that coefficients are estimated
# ONLY in the training data and then applied to test data.
# ----------------------------------------------------------

fit_age_correction <- function(X, age) {
  
  age <- as.numeric(age)
  
  if (anyNA(age)) {
    stop("Missing values in age.")
  }
  
  age_mean <- mean(age)
  
  age_centered <- age - age_mean
  
  denom <- sum(age_centered^2)
  
  if (denom == 0) {
    stop("No variation in age in training data.")
  }
  
  X_mean <- colMeans(X, na.rm = FALSE)
  
  beta <- as.numeric(
    crossprod(age_centered, X) / denom
  )
  
  intercept <- X_mean - beta * age_mean
  
  list(
    intercept = intercept,
    beta = beta
  )
}


# ----------------------------------------------------------
# Apply age correction
# ----------------------------------------------------------

apply_age_correction <- function(X, age, age_model) {
  
  X_corrected <- sweep(
    X,
    2,
    age_model$intercept +
      age_model$beta * age,
    FUN = "-"
  )
  
  X_corrected
}


# ----------------------------------------------------------
# Complete preprocessing for a TRAINING dataset
# ----------------------------------------------------------
#
# Steps:
# 1. age correction
# 2. remove zero-variance features
# 3. select top 5000 by residual variance
#
# Returns the processed matrix AND everything needed to
# apply the same processing to an independent test dataset.
# ----------------------------------------------------------

fit_preprocessing <- function(X, age, n_features = 5000) {
  
  # Age correction
  age_model <- fit_age_correction(
    X = X,
    age = age
  )
  
  X_age_corrected <- apply_age_correction(
    X = X,
    age = age,
    age_model = age_model
  )
  
  # Remove zero / near-zero variance features
  feature_variance <- apply(
    X_age_corrected,
    2,
    var
  )
  
  keep_variance <- is.finite(feature_variance) &
    feature_variance > 1e-6
  
  if (sum(keep_variance) == 0) {
    stop("No features remain after variance filtering.")
  }
  
  X_filtered <- X_age_corrected[
    ,
    keep_variance,
    drop = FALSE
  ]
  
  feature_variance_filtered <- feature_variance[
    keep_variance
  ]
  
  # Number of features cannot exceed available features
  n_select <- min(
    n_features,
    ncol(X_filtered)
  )
  
  top_idx <- order(
    feature_variance_filtered,
    decreasing = TRUE
  )[seq_len(n_select)]
  
  selected_features <- colnames(X_filtered)[top_idx]
  
  X_selected <- X_filtered[
    ,
    selected_features,
    drop = FALSE
  ]
  
  list(
    X = X_selected,
    age_model = age_model,
    variance_features = colnames(X_filtered),
    selected_features = selected_features
  )
}


# ----------------------------------------------------------
# Apply training preprocessing to TEST data
# ----------------------------------------------------------
#
# Critically:
# - age correction uses training coefficients
# - variance selection comes from training
# - top 5000 selection comes from training
# ----------------------------------------------------------

apply_preprocessing <- function(
    X,
    age,
    preprocessing
) {
  
  # Apply training-derived age correction
  X_age_corrected <- apply_age_correction(
    X = X,
    age = age,
    age_model = preprocessing$age_model
  )
  
  # Keep only features that survived training variance filter
  available_features <- intersect(
    preprocessing$variance_features,
    colnames(X_age_corrected)
  )
  
  X_filtered <- X_age_corrected[
    ,
    available_features,
    drop = FALSE
  ]
  
  # Add missing training-selected features as zero
  missing_features <- setdiff(
    preprocessing$selected_features,
    colnames(X_filtered)
  )
  
  if (length(missing_features) > 0) {
    
    missing_matrix <- matrix(
      0,
      nrow = nrow(X_filtered),
      ncol = length(missing_features)
    )
    
    colnames(missing_matrix) <- missing_features
    
    X_filtered <- cbind(
      X_filtered,
      missing_matrix
    )
  }
  
  # Reorder exactly as training data
  X_selected <- X_filtered[
    ,
    preprocessing$selected_features,
    drop = FALSE
  ]
  
  X_selected
}


# ----------------------------------------------------------
# Evaluate pooled predictions
# ----------------------------------------------------------

evaluate_predictions <- function(
    observed,
    predicted
) {
  
  complete <- complete.cases(
    observed,
    predicted
  )
  
  observed <- observed[complete]
  predicted <- predicted[complete]
  
  # R² from pooled predictions
  lm_fit <- lm(
    predicted ~ observed
  )
  
  r2 <- summary(lm_fit)$r.squared
  
  # RMSE
  rmse <- sqrt(
    mean(
      (observed - predicted)^2
    )
  )
  
  # Spearman correlation
  spearman <- suppressWarnings(
    cor(
      observed,
      predicted,
      method = "spearman"
    )
  )
  
  # Pearson correlation
  pearson <- suppressWarnings(
    cor(
      observed,
      predicted,
      method = "pearson"
    )
  )
  
  data.frame(
    R2 = r2,
    RMSE = rmse,
    Spearman = spearman,
    Pearson = pearson,
    N = length(observed)
  )
}


############################################################
# 6. Settings
############################################################

set.seed(500)

n_outer <- 10
n_inner <- 5

n_features <- 5000

# IMPORTANT:
# Include 3080 so we can see where it performs relative
# to other values.
keepX_grid <- c(
  50,
  100,
  200,
  500,
  1000,
  1500,
  2000,
  2500,
  3000,
  3080,
  3500,
  4000
)

ncomp <- 1


############################################################
# 7. Create OUTER folds
############################################################

outer_folds <- createFolds(
  Y,
  k = n_outer,
  list = TRUE,
  returnTrain = FALSE
)

# Storage for final outer predictions
oof_predictions <- rep(
  NA_real_,
  length(Y)
)

# Store which keepX was selected in each outer fold
outer_selected_keepX <- rep(
  NA_integer_,
  n_outer
)

# Store inner-CV performance
inner_results_all <- list()

# Store information about selected features
outer_feature_info <- vector(
  "list",
  n_outer
)


############################################################
# 8. OUTER CV
############################################################

for (outer_i in seq_along(outer_folds)) {
  
  cat("\n")
  cat("========================================\n")
  cat("OUTER FOLD", outer_i, "of", n_outer, "\n")
  cat("========================================\n")
  
  
  # --------------------------------------------------------
  # Outer train/test indices
  # --------------------------------------------------------
  
  outer_test_idx <- outer_folds[[outer_i]]
  
  outer_train_idx <- setdiff(
    seq_along(Y),
    outer_test_idx
  )
  
  
  # --------------------------------------------------------
  # Raw data
  # --------------------------------------------------------
  
  X_outer_train <- X_rclr[
    outer_train_idx,
    ,
    drop = FALSE
  ]
  
  X_outer_test <- X_rclr[
    outer_test_idx,
    ,
    drop = FALSE
  ]
  
  y_outer_train <- Y[
    outer_train_idx
  ]
  
  y_outer_test <- Y[
    outer_test_idx
  ]
  
  age_outer_train <- age[
    outer_train_idx
  ]
  
  age_outer_test <- age[
    outer_test_idx
  ]
  
  
  # --------------------------------------------------------
  # OUTER TRAINING PREPROCESSING
  # --------------------------------------------------------
  #
  # This preprocessing is fitted ONLY on the outer training
  # set.
  #
  # The resulting feature set is used for the final model
  # within this outer fold.
  # --------------------------------------------------------
  
  outer_preprocessing <- fit_preprocessing(
    X = X_outer_train,
    age = age_outer_train,
    n_features = n_features
  )
  
  X_outer_train_processed <- outer_preprocessing$X
  
  X_outer_test_processed <- apply_preprocessing(
    X = X_outer_test,
    age = age_outer_test,
    preprocessing = outer_preprocessing
  )
  
  
  cat(
    "Features after outer preprocessing:",
    ncol(X_outer_train_processed),
    "\n"
  )
  
  
  ##########################################################
  # 9. INNER CV
  ##########################################################
  
  set.seed(500 + outer_i)
  
  inner_folds <- createFolds(
    y_outer_train,
    k = n_inner,
    list = TRUE,
    returnTrain = FALSE
  )
  
  
  # Storage for inner-CV predictions
  #
  # Rows = observations in outer training set
  # Columns = keepX values
  #
  # Each observation gets a prediction from a model that
  # did NOT use that observation during its inner training.
  inner_predictions <- matrix(
    NA_real_,
    nrow = length(y_outer_train),
    ncol = length(keepX_grid)
  )
  
  colnames(inner_predictions) <- paste0(
    "keepX_",
    keepX_grid
  )
  
  
  ##########################################################
  # INNER FOLDS
  ##########################################################
  
  for (inner_i in seq_along(inner_folds)) {
    
    cat(
      "  Inner fold",
      inner_i,
      "of",
      n_inner,
      "\n"
    )
    
    
    inner_test_idx <- inner_folds[[inner_i]]
    
    inner_train_idx <- setdiff(
      seq_along(y_outer_train),
      inner_test_idx
    )
    
    
    # Raw inner training/test data
    X_inner_train <- X_outer_train[
      inner_train_idx,
      ,
      drop = FALSE
    ]
    
    X_inner_test <- X_outer_train[
      inner_test_idx,
      ,
      drop = FALSE
    ]
    
    y_inner_train <- y_outer_train[
      inner_train_idx
    ]
    
    y_inner_test <- y_outer_train[
      inner_test_idx
    ]
    
    age_inner_train <- age_outer_train[
      inner_train_idx
    ]
    
    age_inner_test <- age_outer_train[
      inner_test_idx
    ]
    
    
    # ------------------------------------------------------
    # INNER preprocessing
    # ------------------------------------------------------
    #
    # This is deliberately repeated inside the inner CV.
    #
    # Thus even the top-5000 selection does not see the
    # inner test observations.
    # ------------------------------------------------------
    
    inner_preprocessing <- fit_preprocessing(
      X = X_inner_train,
      age = age_inner_train,
      n_features = n_features
    )
    
    X_inner_train_processed <-
      inner_preprocessing$X
    
    X_inner_test_processed <-
      apply_preprocessing(
        X = X_inner_test,
        age = age_inner_test,
        preprocessing = inner_preprocessing
      )
    
    
    # ------------------------------------------------------
    # Test every keepX
    # ------------------------------------------------------
    
    for (k_i in seq_along(keepX_grid)) {
      
      keepX_value <- keepX_grid[k_i]
      
      # Cannot select more variables than are available
      keepX_actual <- min(
        keepX_value,
        ncol(X_inner_train_processed)
      )
      
      
      # Fit one-component sPLS
      fit <- tryCatch(
        
        spls(
          X_inner_train_processed,
          y_inner_train,
          ncomp = ncomp,
          keepX = keepX_actual
        ),
        
        error = function(e) {
          NULL
        }
      )
      
      
      # If model fails, leave predictions as NA
      if (is.null(fit)) {
        
        cat(
          "    keepX =",
          keepX_value,
          "FAILED\n"
        )
        
        next
      }
      
      
      # Predict inner test observations
      pred_obj <- predict(
        fit,
        X_inner_test_processed
      )
      
      y_pred <- pred_obj$predict[
        ,
        1,
        1
      ]
      
      inner_predictions[
        inner_test_idx,
        k_i
      ] <- y_pred
    }
  }
  
  
  ##########################################################
  # 10. Select best keepX from INNER CV
  ##########################################################
  
  inner_performance <- data.frame(
    keepX = keepX_grid,
    MSE = NA_real_,
    RMSE = NA_real_,
    R2 = NA_real_,
    Spearman = NA_real_
  )
  
  
  for (k_i in seq_along(keepX_grid)) {
    
    pred <- inner_predictions[
      ,
      k_i
    ]
    
    complete <- complete.cases(
      y_outer_train,
      pred
    )
    
    if (sum(complete) < 2) {
      next
    }
    
    y_obs <- y_outer_train[
      complete
    ]
    
    y_pred <- pred[
      complete
    ]
    
    inner_performance$MSE[k_i] <- mean(
      (y_obs - y_pred)^2
    )
    
    inner_performance$RMSE[k_i] <- sqrt(
      inner_performance$MSE[k_i]
    )
    
    inner_performance$R2[k_i] <- summary(
      lm(y_pred ~ y_obs)
    )$r.squared
    
    inner_performance$Spearman[k_i] <-
      suppressWarnings(
        cor(
          y_obs,
          y_pred,
          method = "spearman"
        )
      )
  }
  
  
  # --------------------------------------------------------
  # Choose keepX based on lowest inner-CV MSE
  # --------------------------------------------------------
  
  valid <- which(
    is.finite(inner_performance$MSE)
  )
  
  if (length(valid) == 0) {
    stop(
      paste(
        "No valid inner-CV models in outer fold",
        outer_i
      )
    )
  }
  
  best_idx <- valid[
    which.min(
      inner_performance$MSE[valid]
    )
  ]
  
  best_keepX <- inner_performance$keepX[
    best_idx
  ]
  
  outer_selected_keepX[
    outer_i
  ] <- best_keepX
  
  
  cat(
    "\nSelected keepX for outer fold",
    outer_i,
    ":",
    best_keepX,
    "\n"
  )
  
  
  print(
    inner_performance
  )
  
  
  # Save inner results
  inner_results_all[[outer_i]] <-
    inner_performance
  
  
  ##########################################################
  # 11. Fit final model within this OUTER fold
  ##########################################################
  #
  # Now use ALL outer-training observations.
  #
  # Preprocessing was already fitted using all outer-training
  # observations above.
  ##########################################################
  
  final_outer_model <- tryCatch(
    
    spls(
      X_outer_train_processed,
      y_outer_train,
      ncomp = ncomp,
      keepX = min(
        best_keepX,
        ncol(X_outer_train_processed)
      )
    ),
    
    error = function(e) {
      stop(
        paste(
          "Final outer model failed in fold",
          outer_i,
          ":",
          e$message
        )
      )
    }
  )
  
  
  ##########################################################
  # 12. Predict OUTER TEST SET
  ##########################################################
  
  pred_outer <- predict(
    final_outer_model,
    X_outer_test_processed
  )
  
  y_pred_outer <- pred_outer$predict[
    ,
    1,
    1
  ]
  
  
  # Store predictions in their original observation positions
  oof_predictions[
    outer_test_idx
  ] <- y_pred_outer
  
  
  ##########################################################
  # 13. Store feature information
  ##########################################################
  
  selected_outer_features <-
    selectVar(
      final_outer_model,
      comp = 1
    )$X$name
  
  outer_feature_info[[outer_i]] <-
    selected_outer_features
  
  
  ##########################################################
  # 14. Print fold-specific performance
  ##########################################################
  
  fold_performance <- evaluate_predictions(
    observed = y_outer_test,
    predicted = y_pred_outer
  )
  
  cat(
    "\nOuter fold",
    outer_i,
    "performance:\n"
  )
  
  print(
    fold_performance
  )
}


############################################################
# 15. COMBINED OUT-OF-FOLD PERFORMANCE
############################################################

cat("\n")
cat("========================================\n")
cat("FINAL NESTED CV PERFORMANCE\n")
cat("========================================\n")


# Check that every observation has exactly one prediction
cat(
  "Number of observations:",
  length(Y),
  "\n"
)

cat(
  "Number with OOF predictions:",
  sum(!is.na(oof_predictions)),
  "\n"
)

stopifnot(
  all(!is.na(oof_predictions))
)


# This is the key analysis requested by the reviewer:
#
# all outer-fold predictions are combined first,
# and performance is calculated ONCE.

nested_cv_performance <- evaluate_predictions(
  observed = Y,
  predicted = oof_predictions
)

print(
  nested_cv_performance
)


############################################################
# 16. Save OOF predictions
############################################################

oof_results <- data.frame(
  record_ID = df$record_ID,
  host_age = df$host_age,
  observed = Y,
  predicted = oof_predictions
)

oof_results <- oof_results %>%
  mutate(
    residual = observed - predicted
  )

head(oof_results)


############################################################
# 17. Distribution of selected keepX values
############################################################

keepX_summary <- data.frame(
  outer_fold = seq_len(n_outer),
  selected_keepX = outer_selected_keepX
)

print(
  keepX_summary
)


print(
  table(outer_selected_keepX)
)


############################################################
# 18. Combined inner-CV tuning results
############################################################

inner_results_df <- bind_rows(
  lapply(
    seq_along(inner_results_all),
    function(i) {
      
      inner_results_all[[i]] %>%
        mutate(
          outer_fold = i
        )
    }
  )
)

print(
  inner_results_df
)


############################################################
# 19. Summarize inner-CV performance across outer folds
############################################################

inner_summary <- inner_results_df %>%
  group_by(keepX) %>%
  summarise(
    mean_MSE = mean(MSE, na.rm = TRUE),
    sd_MSE = sd(MSE, na.rm = TRUE),
    
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_RMSE = sd(RMSE, na.rm = TRUE),
    
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    
    mean_Spearman = mean(
      Spearman,
      na.rm = TRUE
    ),
    
    sd_Spearman = sd(
      Spearman,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(
  inner_summary
)


############################################################
# 20. Which keepX was selected most often?
############################################################

keepX_frequency <- as.data.frame(
  table(
    selected_keepX = outer_selected_keepX
  )
)

keepX_frequency$Percentage <-
  100 *
  keepX_frequency$Freq /
  sum(keepX_frequency$Freq)

print(
  keepX_frequency
)


############################################################
# 21. Feature selection stability across OUTER folds
############################################################

# Count how often each metabolite was selected in the
# final model of each outer fold.

all_outer_features <- unique(
  unlist(
    outer_feature_info
  )
)

outer_feature_selection <- matrix(
  0,
  nrow = length(all_outer_features),
  ncol = n_outer
)

rownames(
  outer_feature_selection
) <- all_outer_features

for (i in seq_len(n_outer)) {
  
  selected <- outer_feature_info[[i]]
  
  if (length(selected) > 0) {
    
    outer_feature_selection[
      selected,
      i
    ] <- 1
  }
}


outer_feature_frequency <-
  rowSums(
    outer_feature_selection
  ) / n_outer


nested_feature_stability <- data.frame(
  Metabolite = names(
    outer_feature_frequency
  ),
  Stability = as.numeric(
    outer_feature_frequency
  )
) %>%
  arrange(
    desc(Stability)
  )

head(
  nested_feature_stability,
  20
)


############################################################
# 22. Compare nested CV with original 3080 model
############################################################

cat("\n")
cat("Nested CV results:\n")

cat(
  "R2:",
  round(
    nested_cv_performance$R2,
    3
  ),
  "\n"
)

cat(
  "RMSE:",
  round(
    nested_cv_performance$RMSE,
    3
  ),
  "\n"
)

cat(
  "Spearman:",
  round(
    nested_cv_performance$Spearman,
    3
  ),
  "\n"
)

cat(
  "Pearson:",
  round(
    nested_cv_performance$Pearson,
    3
  ),
  "\n"
)



############################################################
# 24. Save results if desired
############################################################

#saveRDS(list(performance = nested_cv_performance,oof_predictions = oof_results, selected_keepX = keepX_summary,inner_results = inner_results_df,inner_summary = inner_summary,keepX_frequency = keepX_frequency, feature_stability = nested_feature_stability),file = "nested_cv_spls_results.rds")

# ============================================================
# 25. RANDOM-FEATURE NESTED CROSS-VALIDATION
# ============================================================

# The original nested CV is our reference model.
#
# Original model:
#   - age correction fitted within training data
#   - variance filtering
#   - top 5000 metabolites by variance
#   - inner CV tunes keepX
#   - outer CV evaluates prediction
#
# Random-feature model:
#   - age correction fitted within training data
#   - variance filtering
#   - RANDOMLY select 5000 metabolites
#   - inner CV tunes keepX
#   - outer CV evaluates prediction
#
# Everything else is kept identical.

n_features <- 5000
ncomp <- 1

# ------------------------------------------------------------
# 25.1 Function to randomly select candidate metabolites
# ------------------------------------------------------------

fit_random_preprocessing <- function(X, age,
                                     n_features = 5000,
                                     seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Age correction
  age_model <- fit_age_correction(X, age)
  
  X_age_corrected <- apply_age_correction(
    X,
    age,
    age_model
  )
  
  # Remove zero / near-zero variance features
  feature_variance <- apply(
    X_age_corrected,
    2,
    var
  )
  
  keep_variance <- is.finite(feature_variance) &
    feature_variance > 1e-6
  
  if (sum(keep_variance) == 0) {
    stop("No features remaining after variance filtering.")
  }
  
  X_filtered <- X_age_corrected[
    ,
    keep_variance,
    drop = FALSE
  ]
  
  # Randomly select 5000 features
  n_select <- min(
    n_features,
    ncol(X_filtered)
  )
  
  selected_features <- sample(
    colnames(X_filtered),
    size = n_select,
    replace = FALSE
  )
  
  X_selected <- X_filtered[
    ,
    selected_features,
    drop = FALSE
  ]
  
  list(
    X = X_selected,
    age_model = age_model,
    variance_features = colnames(X_filtered),
    selected_features = selected_features
  )
}


# ------------------------------------------------------------
# 25.2 Apply random preprocessing to new data
# ------------------------------------------------------------

apply_random_preprocessing <- function(X, age, preprocessing) {
  
  X_age_corrected <- apply_age_correction(
    X,
    age,
    preprocessing$age_model
  )
  
  # Keep exactly the features selected in the training set
  missing_features <- setdiff(
    preprocessing$selected_features,
    colnames(X_age_corrected)
  )
  
  if (length(missing_features) > 0) {
    X_age_corrected <- cbind(
      X_age_corrected,
      matrix(
        0,
        nrow = nrow(X_age_corrected),
        ncol = length(missing_features),
        dimnames = list(
          rownames(X_age_corrected),
          missing_features
        )
      )
    )
  }
  
  X_selected <- X_age_corrected[
    ,
    preprocessing$selected_features,
    drop = FALSE
  ]
  
  X_selected
}


# ------------------------------------------------------------
# 25.3 Outer folds
# ------------------------------------------------------------

# IMPORTANT:
# Use the EXACT SAME outer folds as the original analysis.
# This makes the comparison between original and random
# features directly paired.

set.seed(500)

outer_folds_random <- createFolds(
  Y,
  k = 10,
  list = TRUE,
  returnTrain = FALSE
)

# Verify that these are identical to the original folds
identical(
  outer_folds,
  outer_folds_random
)

# ============================================================
# 26. RUN RANDOM-FEATURE NESTED CV
# ============================================================

random_oof_predictions <- rep(
  NA_real_,
  length(Y)
)

random_outer_results <- vector(
  "list",
  length(outer_folds_random)
)

random_selected_keepX <- numeric(
  length(outer_folds_random)
)

random_selected_features <- vector(
  "list",
  length(outer_folds_random)
)


for (outer_i in seq_along(outer_folds_random)) {
  
  cat(
    "\n========================================\n",
    "Random-feature outer fold:",
    outer_i,
    "\n========================================\n"
  )
  
  test_idx <- outer_folds_random[[outer_i]]
  
  train_idx <- setdiff(
    seq_along(Y),
    test_idx
  )
  
  X_outer_train <- X_rclr[
    train_idx,
    ,
    drop = FALSE
  ]
  
  X_outer_test <- X_rclr[
    test_idx,
    ,
    drop = FALSE
  ]
  
  y_outer_train <- Y[train_idx]
  y_outer_test <- Y[test_idx]
  
  age_outer_train <- age[train_idx]
  age_outer_test <- age[test_idx]
  
  
  # ----------------------------------------------------------
  # 26.1 Inner CV
  # ----------------------------------------------------------
  
  set.seed(500 + outer_i)
  
  inner_folds <- createFolds(
    y_outer_train,
    k = 5,
    list = TRUE,
    returnTrain = FALSE
  )
  
  inner_results <- matrix(
    NA_real_,
    nrow = length(keepX_grid),
    ncol = length(inner_folds)
  )
  
  rownames(inner_results) <- keepX_grid
  
  
  for (k in seq_along(keepX_grid)) {
    
    keepX_value <- keepX_grid[k]
    
    for (inner_i in seq_along(inner_folds)) {
      
      inner_test_idx <- inner_folds[[inner_i]]
      
      inner_train_idx <- setdiff(
        seq_along(y_outer_train),
        inner_test_idx
      )
      
      X_inner_train <- X_outer_train[
        inner_train_idx,
        ,
        drop = FALSE
      ]
      
      X_inner_test <- X_outer_train[
        inner_test_idx,
        ,
        drop = FALSE
      ]
      
      y_inner_train <- y_outer_train[
        inner_train_idx
      ]
      
      y_inner_test <- y_outer_train[
        inner_test_idx
      ]
      
      age_inner_train <- age_outer_train[
        inner_train_idx
      ]
      
      age_inner_test <- age_outer_train[
        inner_test_idx
      ]
      
      
      # ------------------------------------------------------
      # Random preprocessing is fitted ONLY on inner training
      # ------------------------------------------------------
      
      random_seed <- 100000 +
        outer_i * 1000 +
        inner_i
      
      inner_preprocessing <- fit_random_preprocessing(
        X = X_inner_train,
        age = age_inner_train,
        n_features = n_features,
        seed = random_seed
      )
      
      X_inner_train_processed <-
        inner_preprocessing$X
      
      X_inner_test_processed <-
        apply_random_preprocessing(
          X = X_inner_test,
          age = age_inner_test,
          preprocessing = inner_preprocessing
        )
      
      
      # ------------------------------------------------------
      # sPLS
      # ------------------------------------------------------
      
      set.seed(
        200000 +
          outer_i * 1000 +
          inner_i * 100 +
          k
      )
      
      model <- spls(
        X_inner_train_processed,
        y_inner_train,
        ncomp = ncomp,
        keepX = keepX_value
      )
      
      
      # ------------------------------------------------------
      # Prediction
      # ------------------------------------------------------
      
      prediction <- predict(
        model,
        X_inner_test_processed
      )
      
      y_pred <- prediction$predict[
        ,
        1,
        1
      ]
      
      # MSE used for tuning
      inner_results[k, inner_i] <-
        mean(
          (y_inner_test - y_pred)^2,
          na.rm = TRUE
        )
    }
  }
  
  
  # ----------------------------------------------------------
  # 26.2 Select keepX from inner CV
  # ----------------------------------------------------------
  
  mean_inner_mse <- rowMeans(
    inner_results,
    na.rm = TRUE
  )
  
  best_k <- which.min(
    mean_inner_mse
  )
  
  selected_keepX <- keepX_grid[best_k]
  
  random_selected_keepX[outer_i] <-
    selected_keepX
  
  cat(
    "Selected keepX:",
    selected_keepX,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 26.3 Fit random preprocessing on ALL outer training data
  # ----------------------------------------------------------
  
  outer_random_seed <- 300000 + outer_i
  
  outer_preprocessing <- fit_random_preprocessing(
    X = X_outer_train,
    age = age_outer_train,
    n_features = n_features,
    seed = outer_random_seed
  )
  
  X_outer_train_processed <-
    outer_preprocessing$X
  
  X_outer_test_processed <-
    apply_random_preprocessing(
      X = X_outer_test,
      age = age_outer_test,
      preprocessing = outer_preprocessing
    )
  
  random_selected_features[[outer_i]] <-
    outer_preprocessing$selected_features
  
  
  # ----------------------------------------------------------
  # 26.4 Fit final model on outer training set
  # ----------------------------------------------------------
  
  set.seed(
    400000 + outer_i
  )
  
  final_model <- spls(
    X_outer_train_processed,
    y_outer_train,
    ncomp = ncomp,
    keepX = selected_keepX
  )
  
  
  # ----------------------------------------------------------
  # 26.5 Predict outer test set
  # ----------------------------------------------------------
  
  prediction <- predict(
    final_model,
    X_outer_test_processed
  )
  
  y_pred <- prediction$predict[
    ,
    1,
    1
  ]
  
  random_oof_predictions[
    test_idx
  ] <- y_pred
  
  
  # ----------------------------------------------------------
  # 26.6 Fold-specific performance
  # ----------------------------------------------------------
  
  fold_performance <- evaluate_predictions(
    observed = y_outer_test,
    predicted = y_pred
  )
  
  random_outer_results[[outer_i]] <-
    fold_performance
  
  cat(
    "Fold R²:",
    round(fold_performance$R2, 4),
    "\n"
  )
}

# ============================================================
# 27. RANDOM-FEATURE OOF PERFORMANCE
# ============================================================

random_nested_performance <- evaluate_predictions(
  observed = Y,
  predicted = random_oof_predictions
)

random_nested_performance

