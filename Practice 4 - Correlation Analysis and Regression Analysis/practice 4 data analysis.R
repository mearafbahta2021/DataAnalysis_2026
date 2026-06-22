# ============================================
# COMPLETE ANALYSIS FOR data_for_analysis
# WITHOUT wPerm PACKAGE
# ============================================

# Load required libraries
# Install packages (run once)
install.packages("tidyverse")
install.packages("pROC")
install.packages("corrplot")

# Load packages (run every session)
library(tidyverse)
library(pROC)
library(corrplot)

# Load data
data <- read.csv("data_for_analysis.csv")

# Convert outcome to factor for logistic regression
data$outcome <- as.factor(data$outcome)

# Display structure
summary(data)
str(data)

# ============================================
# 1. CORRELATION ANALYSIS WITH PERMUTATION TEST
# (Manual implementation without wPerm)
# ============================================

# Function for permutation test of Spearman correlation
permutation_correlation <- function(x, y, R = 1000) {
  # Remove NAs
  complete_idx <- complete.cases(x, y)
  x_clean <- x[complete_idx]
  y_clean <- y[complete_idx]
  
  if (length(x_clean) < 5) return(NULL)
  
  # Observed correlation
  obs_cor <- cor(x_clean, y_clean, method = "spearman")
  
  # Permutation test
  perm_cors <- numeric(R)
  for (i in 1:R) {
    y_perm <- sample(y_clean, length(y_clean))
    perm_cors[i] <- cor(x_clean, y_perm, method = "spearman")
  }
  
  # Two-tailed p-value
  p_value <- mean(abs(perm_cors) >= abs(obs_cor))
  
  return(list(observed = obs_cor, p_value = p_value))
}

# Identify numeric variables (exclude outcome and categorical/factor variables)
numeric_vars <- names(data)[sapply(data, is.numeric)]
numeric_vars <- numeric_vars[!numeric_vars %in% c("record_id")]  # Remove ID
cat("Numeric variables for correlation:", paste(numeric_vars, collapse=", "), "\n")
cat("Total numeric variables:", length(numeric_vars), "\n")

# Initialize results table
correlation_results <- data.frame(
  var1 = character(),
  var2 = character(),
  correlation = numeric(),
  p_value = numeric(),
  significance = character(),
  stringsAsFactors = FALSE
)

# Perform pairwise correlation with permutation test
set.seed(123)  # For reproducibility

cat("\nRunning correlation analysis...\n")
total_pairs <- length(numeric_vars) * (length(numeric_vars) - 1) / 2
pair_count <- 0

for (i in 1:(length(numeric_vars)-1)) {
  for (j in (i+1):length(numeric_vars)) {
    pair_count <- pair_count + 1
    if (pair_count %% 50 == 0) {
      cat("Progress:", pair_count, "/", total_pairs, "\n")
    }
    
    var_x <- numeric_vars[i]
    var_y <- numeric_vars[j]
    
    # Perform permutation correlation
    perm_result <- permutation_correlation(data[[var_x]], data[[var_y]], R = 500)
    
    if (!is.null(perm_result)) {
      # Store results
      correlation_results <- rbind(correlation_results, data.frame(
        var1 = var_x,
        var2 = var_y,
        correlation = perm_result$observed,
        p_value = perm_result$p_value,
        significance = ifelse(perm_result$p_value < 0.001, "***",
                              ifelse(perm_result$p_value < 0.01, "**",
                                     ifelse(perm_result$p_value < 0.05, "*", "ns"))),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Display correlation results
cat("\n========== CORRELATION ANALYSIS RESULTS ==========\n")
print(head(correlation_results, 20))

# Save full results
write.csv(correlation_results, "correlation_results_full.csv", row.names = FALSE)

# Filter significant correlations (p < 0.05)
significant_cors <- correlation_results[correlation_results$p_value < 0.05, ]
significant_cors <- significant_cors[order(abs(significant_cors$correlation), decreasing = TRUE), ]

cat("\n========== SIGNIFICANT CORRELATIONS (p < 0.05) ==========\n")
cat("Total significant correlations:", nrow(significant_cors), "\n")
print(head(significant_cors, 30))

# Create correlation matrix for visualization
cor_matrix <- cor(data[, numeric_vars], method = "spearman", use = "pairwise.complete.obs")

# Visualize correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.5,
         title = "Spearman Correlation Matrix", mar = c(0,0,1,0))

# ============================================
# 2. REGRESSION ANALYSIS BETWEEN OTHER VARIABLES
# Select best model by BIC
# ============================================

# Prepare data for regression (using lipids1 as dependent variable)
reg_data <- data[, c("lipids1", "lipids2", "lipids3", "lipids4", "lipids5", 
                     "hormone1", "hormone2", "hormone3", "hormone4")]
reg_data <- na.omit(reg_data)

# Define candidate models
models <- list()

# Model 1: Simple linear (lipids1 ~ lipids2)
models[["Linear"]] <- lm(lipids1 ~ lipids2, data = reg_data)

# Model 2: Quadratic
models[["Quadratic"]] <- lm(lipids1 ~ poly(lipids2, 2), data = reg_data)

# Model 3: Cubic
models[["Cubic"]] <- lm(lipids1 ~ poly(lipids2, 3), data = reg_data)

# Model 4: Multiple linear (all lipids)
models[["Multiple_Lipids"]] <- lm(lipids1 ~ lipids2 + lipids3 + lipids4 + lipids5, data = reg_data)

# Model 5: Multiple linear (all hormones)
models[["Multiple_Hormones"]] <- lm(lipids1 ~ hormone1 + hormone2 + hormone3 + hormone4, data = reg_data)

# Model 6: Full model (all predictors)
models[["Full"]] <- lm(lipids1 ~ ., data = reg_data)

# Model 7: Stepwise selected (using AIC)
full_model <- lm(lipids1 ~ ., data = reg_data)
models[["Stepwise"]] <- step(full_model, direction = "both", trace = 0)

# Model 8: Log transformation of response
models[["Log_Response"]] <- lm(log(lipids1) ~ lipids2, data = reg_data)

# Model 9: Log transformation of predictor
models[["Log_Predictor"]] <- lm(lipids1 ~ log(lipids2 + 0.1), data = reg_data)

# Compare models by BIC
bic_results <- data.frame(
  Model = names(models),
  BIC = sapply(models, BIC),
  AIC = sapply(models, AIC),
  R_squared = sapply(models, function(m) summary(m)$r.squared),
  Adj_R_squared = sapply(models, function(m) summary(m)$adj.r.squared)
)

# Sort by BIC (lower is better)
bic_results <- bic_results[order(bic_results$BIC), ]

cat("\n========== REGRESSION MODEL COMPARISON (by BIC) ==========\n")
print(bic_results)

# Best model
best_model_name <- bic_results$Model[1]
best_model <- models[[best_model_name]]

cat("\n========== BEST MODEL ==========\n")
cat("Best model:", best_model_name, "\n")
cat("BIC:", round(bic_results$BIC[1], 2), "\n")
cat("R-squared:", round(bic_results$R_squared[1], 4), "\n")
cat("Adjusted R-squared:", round(bic_results$Adj_R_squared[1], 4), "\n")
print(summary(best_model))

# ============================================
# 3. LOGISTIC REGRESSION USING HORMONE VARIABLES (CORRECTED)
# Predict binary outcome, compare by AIC/BIC, compute odds ratios
# ============================================

# Identify hormone variables (CORRECTED - no hormone9)
hormone_vars <- c("hormone1", "hormone2", "hormone3", "hormone4", 
                  "hormone5", "hormone6", "hormone7", "hormone8", "hormone10_generated")

# Prepare data for logistic regression
logit_data <- data[, c("outcome", hormone_vars)]
logit_data <- na.omit(logit_data)
logit_data$outcome <- as.numeric(as.character(logit_data$outcome))

# Check class balance
cat("\n========== OUTCOME DISTRIBUTION ==========\n")
print(table(logit_data$outcome))
cat("Proportion of cases:", prop.table(table(logit_data$outcome))[2], "\n")

# Define logistic regression models
logit_models <- list()

# Model 1: Single hormone (hormone1)
logit_models[["Hormone1_only"]] <- glm(outcome ~ hormone1, data = logit_data, family = binomial)

# Model 2: Single hormone (hormone2)
logit_models[["Hormone2_only"]] <- glm(outcome ~ hormone2, data = logit_data, family = binomial)

# Model 3: Single hormone (hormone3)
logit_models[["Hormone3_only"]] <- glm(outcome ~ hormone3, data = logit_data, family = binomial)

# Model 4: Single hormone (hormone4)
logit_models[["Hormone4_only"]] <- glm(outcome ~ hormone4, data = logit_data, family = binomial)

# Model 5: Single hormone (hormone5)
logit_models[["Hormone5_only"]] <- glm(outcome ~ hormone5, data = logit_data, family = binomial)

# Model 6: Single hormone (hormone6)
logit_models[["Hormone6_only"]] <- glm(outcome ~ hormone6, data = logit_data, family = binomial)

# Model 7: Single hormone (hormone7)
logit_models[["Hormone7_only"]] <- glm(outcome ~ hormone7, data = logit_data, family = binomial)

# Model 8: Single hormone (hormone8)
logit_models[["Hormone8_only"]] <- glm(outcome ~ hormone8, data = logit_data, family = binomial)

# Model 9: Single hormone (hormone10_generated)
logit_models[["Hormone10_only"]] <- glm(outcome ~ hormone10_generated, data = logit_data, family = binomial)

# Model 10: Two hormones (hormone1 + hormone2)
logit_models[["Two_Hormones_H1H2"]] <- glm(outcome ~ hormone1 + hormone2, data = logit_data, family = binomial)

# Model 11: Two hormones (hormone3 + hormone4)
logit_models[["Two_Hormones_H3H4"]] <- glm(outcome ~ hormone3 + hormone4, data = logit_data, family = binomial)

# Model 12: Three hormones (hormone1, hormone2, hormone3)
logit_models[["Three_Hormones"]] <- glm(outcome ~ hormone1 + hormone2 + hormone3, data = logit_data, family = binomial)

# Model 13: All hormones
logit_models[["All_Hormones"]] <- glm(outcome ~ ., data = logit_data, family = binomial)

# Model 14: Stepwise selected
step_model <- step(logit_models[["All_Hormones"]], direction = "both", trace = 0)
logit_models[["Stepwise"]] <- step_model

# Model 15: Interaction model (hormone1 * hormone2)
logit_models[["Interaction_H1H2"]] <- glm(outcome ~ hormone1 * hormone2, data = logit_data, family = binomial)

# Model 16: Hormone1 + hormone3
logit_models[["H1_H3"]] <- glm(outcome ~ hormone1 + hormone3, data = logit_data, family = binomial)

# Model 17: Hormone2 + hormone4
logit_models[["H2_H4"]] <- glm(outcome ~ hormone2 + hormone4, data = logit_data, family = binomial)

# Model 18: Hormone5 + hormone6
logit_models[["H5_H6"]] <- glm(outcome ~ hormone5 + hormone6, data = logit_data, family = binomial)

# Model 19: Hormone7 + hormone8
logit_models[["H7_H8"]] <- glm(outcome ~ hormone7 + hormone8, data = logit_data, family = binomial)

# Compare model performance
logit_comparison <- data.frame(
  Model = names(logit_models),
  AIC = sapply(logit_models, AIC),
  BIC = sapply(logit_models, BIC),
  LogLik = sapply(logit_models, logLik),
  DF = sapply(logit_models, function(m) length(coef(m)))
)

# Calculate McFadden's pseudo R-squared
null_model <- glm(outcome ~ 1, data = logit_data, family = binomial)
for (i in 1:length(logit_models)) {
  logit_comparison$McFadden_R2[i] <- 1 - (as.numeric(logLik(logit_models[[i]])) / as.numeric(logLik(null_model)))
}

# Sort by BIC (lower is better)
logit_comparison <- logit_comparison[order(logit_comparison$BIC), ]

cat("\n========== LOGISTIC REGRESSION MODEL COMPARISON ==========\n")
print(logit_comparison)

# Select best model (lowest BIC)
best_logit_name <- logit_comparison$Model[1]
best_logit_model <- logit_models[[best_logit_name]]

cat("\n========== BEST LOGISTIC MODEL ==========\n")
cat("Best model by BIC:", best_logit_name, "\n")
cat("AIC:", round(logit_comparison$AIC[1], 2), "\n")
cat("BIC:", round(logit_comparison$BIC[1], 2), "\n")
cat("McFadden R²:", round(logit_comparison$McFadden_R2[1], 4), "\n")
print(summary(best_logit_model))

# ============================================
# 4. ODDS RATIOS FOR BEST MODEL
# ============================================

# Calculate odds ratios and confidence intervals
odds_ratios <- exp(coef(best_logit_model))
ci_odds <- exp(confint(best_logit_model))

odds_ratio_table <- data.frame(
  Predictor = names(odds_ratios),
  Odds_Ratio = odds_ratios,
  CI_lower = ci_odds[, 1],
  CI_upper = ci_odds[, 2],
  Coefficient = coef(best_logit_model),
  Std_Error = summary(best_logit_model)$coefficients[, 2],
  Z_value = summary(best_logit_model)$coefficients[, 3],
  P_value = summary(best_logit_model)$coefficients[, 4]
)

# Add significance stars
odds_ratio_table$Significance <- ifelse(odds_ratio_table$P_value < 0.001, "***",
                                        ifelse(odds_ratio_table$P_value < 0.01, "**",
                                               ifelse(odds_ratio_table$P_value < 0.05, "*", "ns")))

cat("\n========== ODDS RATIOS FOR BEST MODEL ==========\n")
print(odds_ratio_table)
# ============================================
# 5. MODEL PERFORMANCE METRICS FOR BEST MODEL (ROBUST VERSION)
# ============================================

# Predict probabilities
predicted_prob <- predict(best_logit_model, type = "response")
predicted_class <- ifelse(predicted_prob > 0.5, 1, 0)

# Confusion matrix with all possible levels
conf_matrix <- table(Actual = factor(logit_data$outcome, levels = c(0,1)), 
                     Predicted = factor(predicted_class, levels = c(0,1)))
cat("\n========== CONFUSION MATRIX ==========\n")
print(conf_matrix)

# Initialize metrics with NA values
accuracy <- NA
sensitivity <- NA
specificity <- NA
precision <- NA
negative_pred_value <- NA

# Calculate metrics only if both classes are present in confusion matrix
if (nrow(conf_matrix) == 2 & ncol(conf_matrix) == 2) {
  accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  sensitivity <- conf_matrix[2,2] / sum(conf_matrix[2,])  # True positive rate
  specificity <- conf_matrix[1,1] / sum(conf_matrix[1,])  # True negative rate
  precision <- conf_matrix[2,2] / sum(conf_matrix[,2])    # Positive predictive value
  negative_pred_value <- conf_matrix[1,1] / sum(conf_matrix[1,])  # NPV
  
  cat("\nConfusion matrix has both classes.\n")
} else {
  cat("\nWARNING: Confusion matrix has only one class.\n")
  cat("The model predicted all cases as the same outcome.\n")
  cat("Confusion matrix dimensions:", dim(conf_matrix), "\n")
}

# ROC and AUC (robust)
if (length(unique(logit_data$outcome)) > 1 & length(unique(predicted_prob)) > 1) {
  roc_curve <- roc(logit_data$outcome, predicted_prob, quiet = TRUE)
  auc_value <- auc(roc_curve)
  
  # Plot ROC curve
  plot(roc_curve, main = paste("ROC Curve - AUC =", round(auc_value, 3)),
       col = "blue", lwd = 2)
  abline(a=0, b=1, lty=2, col="gray")
} else {
  auc_value <- NA
  cat("\nWARNING: Cannot compute ROC curve - need both classes and varied predictions.\n")
}

# Create performance metrics data frame
performance_metrics <- data.frame(
  Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", 
             "Precision (PPV)", "Negative Predictive Value", "AUC"),
  Value = c(accuracy, sensitivity, specificity, precision, negative_pred_value, as.numeric(auc_value))
)

cat("\n========== PERFORMANCE METRICS ==========\n")
print(performance_metrics)

# ============================================
# 6. ADDITIONAL METRICS FOR ALL MODELS (ROBUST)
# ============================================

# Calculate additional metrics for all models
logit_comparison$Accuracy <- NA
logit_comparison$AUC <- NA

for (i in 1:length(logit_models)) {
  # Predict
  pred_prob <- predict(logit_models[[i]], type = "response")
  pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  
  # Confusion matrix with both levels
  cm <- table(Actual = factor(logit_data$outcome, levels = c(0,1)), 
              Predicted = factor(pred_class, levels = c(0,1)))
  
  # Calculate accuracy only if both classes are present
  if (nrow(cm) == 2 & ncol(cm) == 2) {
    logit_comparison$Accuracy[i] <- sum(diag(cm)) / sum(cm)
    
    # Calculate AUC
    if (length(unique(logit_data$outcome)) > 1 & length(unique(pred_prob)) > 1) {
      roc_obj <- tryCatch(roc(logit_data$outcome, pred_prob, quiet = TRUE), 
                          error = function(e) NULL)
      if (!is.null(roc_obj)) {
        logit_comparison$AUC[i] <- as.numeric(auc(roc_obj))
      }
    }
  }
}

# Reorder by BIC and display
logit_comparison <- logit_comparison[order(logit_comparison$BIC), ]

cat("\n========== COMPLETE MODEL COMPARISON WITH PERFORMANCE ==========\n")
print(logit_comparison)

# ============================================
# 7. ENHANCED SUMMARY WITH INTERPRETATION
# ============================================

cat("\n")
cat("==================== LOGISTIC REGRESSION SUMMARY ====================\n")
cat("\n1. DATA OVERVIEW:\n")
cat("   Total observations:", nrow(logit_data), "\n")
cat("   Outcome distribution:\n")
print(table(logit_data$outcome))
cat("   Proportion of cases (outcome=1):", 
    round(prop.table(table(logit_data$outcome))[2], 3), "\n")

cat("\n2. BEST MODEL SELECTION (by BIC):\n")
cat("   Model:", best_logit_name, "\n")
cat("   BIC:", round(min(logit_comparison$BIC, na.rm = TRUE), 2), "\n")
cat("   AIC:", round(logit_comparison$AIC[1], 2), "\n")
cat("   McFadden R²:", round(logit_comparison$McFadden_R2[1], 4), "\n")

if (!is.na(logit_comparison$Accuracy[1])) {
  cat("   Accuracy:", round(logit_comparison$Accuracy[1], 3), "\n")
}
if (!is.na(logit_comparison$AUC[1])) {
  cat("   AUC:", round(logit_comparison$AUC[1], 3), "\n")
}

cat("\n3. MODEL COEFFICIENTS (BEST MODEL):\n")
coef_summary <- summary(best_logit_model)$coefficients
print(coef_summary)

cat("\n4. ODDS RATIOS WITH 95% CI:\n")
# Recalculate odds ratios safely
odds_ratios <- exp(coef(best_logit_model))
ci_odds <- tryCatch(exp(confint(best_logit_model)), error = function(e) {
  # If confint fails, use normal approximation
  se <- summary(best_logit_model)$coefficients[, 2]
  coefs <- coef(best_logit_model)
  lower <- exp(coefs - 1.96 * se)
  upper <- exp(coefs + 1.96 * se)
  return(cbind(lower, upper))
})

odds_ratio_table <- data.frame(
  Predictor = names(odds_ratios),
  Odds_Ratio = odds_ratios,
  CI_lower = ci_odds[, 1],
  CI_upper = ci_odds[, 2],
  P_value = summary(best_logit_model)$coefficients[, 4]
)

# Add significance stars
odds_ratio_table$Significance <- ifelse(odds_ratio_table$P_value < 0.001, "***",
                                        ifelse(odds_ratio_table$P_value < 0.01, "**",
                                               ifelse(odds_ratio_table$P_value < 0.05, "*", "ns")))

print(odds_ratio_table)

cat("\n5. SIGNIFICANT PREDICTORS (p < 0.05):\n")
significant_ors <- odds_ratio_table[odds_ratio_table$P_value < 0.05 & odds_ratio_table$Predictor != "(Intercept)", ]
if(nrow(significant_ors) > 0) {
  for(i in 1:nrow(significant_ors)) {
    cat(sprintf("   %s: OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f %s\n",
                significant_ors$Predictor[i],
                significant_ors$Odds_Ratio[i],
                significant_ors$CI_lower[i],
                significant_ors$CI_upper[i],
                significant_ors$P_value[i],
                significant_ors$Significance[i]))
  }
} else {
  cat("   No statistically significant predictors found (p < 0.05).\n")
}

cat("\n6. TOP 5 MODELS BY BIC:\n")
top5 <- head(logit_comparison[, c("Model", "BIC", "AIC", "McFadden_R2")], 5)
print(top5)

# Add interpretation for McFadden R²
cat("\n7. INTERPRETATION:\n")
cat("   McFadden's R² interpretation:\n")
cat("   - 0.2-0.4: Excellent fit (similar to R² of 0.7-0.9 in OLS)\n")
cat("   - 0.1-0.2: Good fit\n")
cat("   - 0.05-0.1: Acceptable fit\n")
cat("   - <0.05: Poor fit\n")

best_r2 <- logit_comparison$McFadden_R2[1]
if(best_r2 > 0.2) {
  cat(sprintf("   -> Current best model McFadden R² = %.4f (EXCELLENT fit)\n", best_r2))
} else if(best_r2 > 0.1) {
  cat(sprintf("   -> Current best model McFadden R² = %.4f (GOOD fit)\n", best_r2))
} else if(best_r2 > 0.05) {
  cat(sprintf("   -> Current best model McFadden R² = %.4f (ACCEPTABLE fit)\n", best_r2))
} else {
  cat(sprintf("   -> Current best model McFadden R² = %.4f (POOR fit)\n", best_r2))
}

cat("\n==================== END OF SUMMARY ====================\n")

# Save results
write.csv(logit_comparison, "logistic_model_comparison.csv", row.names = FALSE)
write.csv(odds_ratio_table, "odds_ratios.csv", row.names = FALSE)
write.csv(performance_metrics, "performance_metrics.csv", row.names = FALSE)

cat("\nResults saved to CSV files:\n")
cat("  - logistic_model_comparison.csv\n")
cat("  - odds_ratios.csv\n")
cat("  - performance_metrics.csv\n")