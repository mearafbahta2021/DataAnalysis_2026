#--------------------start-------------------------------
# Get current working directory
getwd()

#----------------read dataset--------------------------
list.files()

# For the distribution.csv file
example_df <- read.csv("distribution.csv", header = TRUE, dec = ',', sep = ";")

# For factor_data.csv
factor_df <- read.csv("factor_data.csv")

# For imputed_data.csv
imputed_df <- read.csv("imputed_data.csv")
# Display structure with variable types
str(example_df)
str(factor_df)
str(imputed_df)

#---------------merge two files-------------------------
data_for_analysis <- merge(
  factor_df, 
  imputed_df, 
  by = "record_id",        # column for merge
  all = FALSE       # FALSE = INNER JOIN (only coincidences), TRUE = FULL JOIN
)
str(data_for_analysis)

# save data_for_analysis in CSV
write.csv(data_for_analysis, "data_for_analysis.csv", row.names = FALSE)  

#------------------Probability Distributions----------------------- 
install.packages("MASS", dependencies = T)
library(MASS)

#----------------example-----------------------------------------
summary(example_df)
example_df$value <- as.numeric(example_df$value)
summary(example_df)

#building histograms for example
# normal distribution
val <- example_df[example_df$distribution == "norm", ]$value

mean(val)
sd(val)
hist(val)

fit <- fitdistr(val, densfun = "normal")
fit

#lognormal distribution
val <- example_df[example_df$distribution == "lognorm", ]$value

mean(val)
sd(val)
hist(val)

fit <- fitdistr(val, densfun = "lognormal")
fit

unname(fit$estimate[1])
unname(fit$estimate[2])

m_log <- exp(unname(fit$estimate[1])) * sqrt(exp(unname(fit$estimate[2])^2))
m_log
sd_log <- sqrt(exp(2 * unname(fit$estimate[1])) * (exp(unname(fit$estimate[2])^2) - 1) * sqrt(exp(unname(fit$estimate[2])^2)))
sd_log

#exponential distribution
val <- example_df[example_df$distribution == "exp", ]$value

mean(val)
sd(val)
hist(val)

fit <- fitdistr(val, densfun = "exponential")
fit

unname(fit$estimate[1])

m_exp <- 1 / unname(fit$estimate[1])
m_exp

#Poisson distribution
val <- example_df[example_df$distribution == "pois", ]$value

mean(val)
sd(val)
hist(val)

fit <- fitdistr(val, densfun = "Poisson")
fit

unname(fit$estimate[1])

sd_pois <- sqrt(unname(fit$estimate[1]))
sd_pois

#Selecting a Distribution Model
val <- example_df[example_df$distribution == "lognorm", ]$value

fit_1 <- fitdistr(val, densfun = "normal")
fit_2 <- fitdistr(val, densfun = "lognormal")
fit_3 <- fitdistr(val, densfun = "exponential")

#Bayesian Information Criterion calculation
BIC(fit_3)

#calculation of the Bayesian information criterion for all models
BIC_value <- c(BIC(fit_1), BIC(fit_2), BIC(fit_3))

#forming a vector with the name of the models
distribution <- c("normal", "lognormal", "exponential")

#combining the results into a final table
rez <- data.frame(BIC_value = BIC_value, distribution = distribution)

#sort table in ascending order of Bayesian Information Criterion value
rez <- rez[order(rez$BIC_value, decreasing = FALSE), ]
rez

#calculation of absolute values of the confidence interval for the mean of a lognormal distribution
error_min <- unname(fit_2$estimate[1]) - unname(fit_2$sd[1])
error_max <- unname(fit_2$estimate[1]) + unname(fit_2$sd[1])

error_min
error_max

m <- exp(unname(fit_2$estimate[1])) * sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_min <- exp(error_min) * sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_max <- exp(error_max) * sqrt(exp(unname(fit_2$estimate[2])^2))

value_error_min
m
value_error_max

#--------------data for analysis--------------------------
#building histograms
value_d1 <- data_for_analysis$lipids1
hist(value_d1)
value_d2 <- data_for_analysis$lipids2
hist(value_d2)
value_d3 <- data_for_analysis$lipids3
hist(value_d3)
value_d4 <- data_for_analysis$lipids4
hist(value_d4)

# d1 distribution estimate
fit_d1_1 <- fitdistr(value_d1, densfun = "normal")
fit_d1_2 <- fitdistr(value_d1, densfun = "lognormal")
fit_d1_3 <- fitdistr(value_d1, densfun = "exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d1
BIC_value_d1 <- c(BIC(fit_d1_1), BIC(fit_d1_2), BIC(fit_d1_3))
distribution <- c("normal", "lognormal", "exponential")
result_d1 <- data.frame(BIC_value_d1 = BIC_value_d1, distribution = distribution)
result_d1
min(result_d1$BIC_value_d1)
distribution_d1 <- result_d1[result_d1$BIC_value_d1 == min(result_d1$BIC_value_d1), ]$distribution
distribution_d1

# Finding parameters for d1
fit_d1_1$estimate[1:2]

# d2 distribution estimate
fit_d2_1 <- fitdistr(value_d2, densfun = "normal")
fit_d2_2 <- fitdistr(value_d2, densfun = "lognormal")
fit_d2_3 <- fitdistr(value_d2, densfun = "exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d2
BIC_value_d2 <- c(BIC(fit_d2_1), BIC(fit_d2_2), BIC(fit_d2_3))
distribution <- c("normal", "lognormal", "exponential")
result_d2 <- data.frame(BIC_value_d2 = BIC_value_d2, distribution = distribution)
result_d2
min(result_d2$BIC_value_d2)
distribution_d2 <- result_d2[result_d2$BIC_value_d2 == min(result_d2$BIC_value_d2), ]$distribution
distribution_d2

# Finding parameters for d2
fit_d2_1$estimate[1:2]

# ================ FIX MISSING DATA FOR LIPIDS5 (EXTRA POINTS) ================
# Check for missing values in lipids5
if("lipids5" %in% names(data_for_analysis)) {
  missing_count <- sum(is.na(data_for_analysis$lipids5))
  cat("\nMissing values in lipids5 before fixing:", missing_count, "\n")
  
  if(missing_count > 0) {
    # Fix missing data by imputing with the mean
    data_for_analysis$lipids5[is.na(data_for_analysis$lipids5)] <- mean(data_for_analysis$lipids5, na.rm = TRUE)
    cat("Missing values in lipids5 after fixing:", sum(is.na(data_for_analysis$lipids5)), "\n")
  }
}

# ================ DISTRIBUTION ESTIMATION BY GROUP ================
# Identify continuous variables (lipid measurements)
lipid_vars <- names(data_for_analysis)[grep("lipids", names(data_for_analysis))]
print(paste("Lipid variables:", paste(lipid_vars, collapse = ", ")))

# CORRECTED FUNCTION - with consistent column structure
estimate_distribution <- function(values, var_name, group_name = NULL) {
  # Remove NA values
  values_clean <- values[!is.na(values)]
  
  if(length(values_clean) < 2) {
    return(data.frame(
      Variable = var_name,
      Group = ifelse(is.null(group_name), "All", group_name),
      Distribution = "Insufficient data",
      Parameter1_Name = NA_character_,
      Parameter1 = NA_real_,
      Parameter2_Name = NA_character_,
      Parameter2 = NA_real_,
      Mean = mean(values_clean, na.rm = TRUE),
      SD = sd(values_clean, na.rm = TRUE),
      N = length(values_clean),
      BIC = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # Try fitting different distributions
  tryCatch({
    fit_normal <- fitdistr(values_clean, densfun = "normal")
    fit_lognormal <- fitdistr(values_clean, densfun = "lognormal")
    fit_exponential <- fitdistr(values_clean, densfun = "exponential")
    
    # Calculate BIC
    bic_values <- c(BIC(fit_normal), BIC(fit_lognormal), BIC(fit_exponential))
    distributions <- c("normal", "lognormal", "exponential")
    
    # Select best distribution (lowest BIC)
    best_idx <- which.min(bic_values)
    best_dist <- distributions[best_idx]
    
    # Extract parameters based on best distribution
    if(best_dist == "normal") {
      return(data.frame(
        Variable = var_name,
        Group = ifelse(is.null(group_name), "All", group_name),
        Distribution = "normal",
        Parameter1_Name = "mean",
        Parameter1 = unname(fit_normal$estimate[1]),
        Parameter2_Name = "sd",
        Parameter2 = unname(fit_normal$estimate[2]),
        Mean = mean(values_clean, na.rm = TRUE),
        SD = sd(values_clean, na.rm = TRUE),
        N = length(values_clean),
        BIC = bic_values[1],
        stringsAsFactors = FALSE
      ))
    } else if(best_dist == "lognormal") {
      return(data.frame(
        Variable = var_name,
        Group = ifelse(is.null(group_name), "All", group_name),
        Distribution = "lognormal",
        Parameter1_Name = "meanlog",
        Parameter1 = unname(fit_lognormal$estimate[1]),
        Parameter2_Name = "sdlog",
        Parameter2 = unname(fit_lognormal$estimate[2]),
        Mean = mean(values_clean, na.rm = TRUE),
        SD = sd(values_clean, na.rm = TRUE),
        N = length(values_clean),
        BIC = bic_values[2],
        stringsAsFactors = FALSE
      ))
    } else { # exponential
      return(data.frame(
        Variable = var_name,
        Group = ifelse(is.null(group_name), "All", group_name),
        Distribution = "exponential",
        Parameter1_Name = "rate",
        Parameter1 = unname(fit_exponential$estimate[1]),
        Parameter2_Name = NA_character_,
        Parameter2 = NA_real_,
        Mean = mean(values_clean, na.rm = TRUE),
        SD = sd(values_clean, na.rm = TRUE),
        N = length(values_clean),
        BIC = bic_values[3],
        stringsAsFactors = FALSE
      ))
    }
  }, error = function(e) {
    return(data.frame(
      Variable = var_name,
      Group = ifelse(is.null(group_name), "All", group_name),
      Distribution = paste("Error"),
      Parameter1_Name = NA_character_,
      Parameter1 = NA_real_,
      Parameter2_Name = NA_character_,
      Parameter2 = NA_real_,
      Mean = mean(values_clean, na.rm = TRUE),
      SD = sd(values_clean, na.rm = TRUE),
      N = length(values_clean),
      BIC = NA_real_,
      stringsAsFactors = FALSE
    ))
  })
}

# Check if we have an outcome/grouping variable
group_var <- NULL
if("outcome" %in% names(data_for_analysis)) {
  group_var <- "outcome"
} else if("group" %in% names(data_for_analysis)) {
  group_var <- "group"
} else if("treatment" %in% names(data_for_analysis)) {
  group_var <- "treatment"
}

# Create descriptive statistics table by group
results_list <- list()

if(!is.null(group_var) && length(unique(data_for_analysis[[group_var]])) > 1) {
  # Estimate distribution for each lipid variable by group
  for(lipid in lipid_vars) {
    for(grp in unique(data_for_analysis[[group_var]])) {
      group_values <- data_for_analysis[[lipid]][data_for_analysis[[group_var]] == grp]
      results_list[[length(results_list) + 1]] <- estimate_distribution(
        group_values, lipid, as.character(grp)
      )
    }
  }
  
  # Also calculate overall statistics (without grouping)
  for(lipid in lipid_vars) {
    results_list[[length(results_list) + 1]] <- estimate_distribution(
      data_for_analysis[[lipid]], lipid, "Overall"
    )
  }
} else {
  # If no grouping variable, just estimate overall distributions
  for(lipid in lipid_vars) {
    results_list[[length(results_list) + 1]] <- estimate_distribution(
      data_for_analysis[[lipid]], lipid, NULL
    )
  }
}

# Combine all results
final_results <- do.call(rbind, results_list)

# Display the results table
print("========== DESCRIPTIVE STATISTICS AND DISTRIBUTION PARAMETERS ==========")
print(final_results)

# Save results to CSV
write.csv(final_results, "distribution_analysis_results.csv", row.names = FALSE)

# ================ FIXED VISUALIZATION WITH MARGIN CONTROL ================
# Reset graphics parameters to default first
dev.off()  # Close any open graphics devices

# Set smaller margins for cloud environment
par(mar = c(4, 4, 2, 1))  # bottom, left, top, right margins
par(mfrow = c(1, 1))  # Reset to single plot

# Create histograms for each lipid variable by group with error handling
if(!is.null(group_var) && length(unique(data_for_analysis[[group_var]])) > 1) {
  # Check number of groups and adjust layout
  n_groups <- length(unique(data_for_analysis[[group_var]]))
  n_plots <- length(lipid_vars)
  
  # Calculate appropriate layout
  if(n_plots <= 2) {
    par(mfrow = c(1, n_plots))
  } else if(n_plots <= 4) {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(ceiling(n_plots/2), 2))
  }
  
  for(lipid in lipid_vars) {
    tryCatch({
      # Get min and max for consistent x-axis
      xlim_range <- range(data_for_analysis[[lipid]], na.rm = TRUE)
      
      # Create histogram for each group
      groups <- unique(data_for_analysis[[group_var]])
      colors <- rainbow(length(groups))
      
      # First plot to set up
      first_group <- groups[1]
      hist(data_for_analysis[[lipid]][data_for_analysis[[group_var]] == first_group],
           main = paste(lipid, "by", group_var),
           xlab = lipid, col = colors[1], 
           xlim = xlim_range, 
           probability = TRUE,
           border = "white")
      
      # Add other groups
      for(i in 2:length(groups)) {
        hist(data_for_analysis[[lipid]][data_for_analysis[[group_var]] == groups[i]],
             add = TRUE, col = colors[i], probability = TRUE,
             border = "white")
      }
      
      # Add legend with smaller size
      legend("topright", legend = groups, fill = colors, cex = 0.8)
    }, error = function(e) {
      cat("Could not plot", lipid, ":", e$message, "\n")
    })
  }
} else {
  # Simple histograms if no grouping
  n_plots <- length(lipid_vars)
  if(n_plots <= 2) {
    par(mfrow = c(1, n_plots))
  } else if(n_plots <= 4) {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(ceiling(n_plots/2), 2))
  }
  
  for(lipid in lipid_vars) {
    tryCatch({
      hist(data_for_analysis[[lipid]], main = lipid, xlab = lipid, 
           col = "lightblue", probability = TRUE)
      
      # Add density curve
      lines(density(data_for_analysis[[lipid]], na.rm = TRUE), col = "red", lwd = 2)
    }, error = function(e) {
      cat("Could not plot", lipid, ":", e$message, "\n")
    })
  }
}

# Reset to single plot
par(mfrow = c(1, 1))

# ================ SUMMARY STATISTICS TABLE ================
# Create a comprehensive summary table
summary_table <- data.frame(
  Variable = character(),
  Group = character(),
  N = integer(),
  Mean = numeric(),
  Median = numeric(),
  SD = numeric(),
  Min = numeric(),
  Max = numeric(),
  Skewness = numeric(),
  Best_Distribution = character(),
  stringsAsFactors = FALSE
)

# Function to calculate skewness
skewness <- function(x) {
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  if(n < 3) return(NA)
  (sum((x_clean - mean(x_clean))^3) / n) / (sum((x_clean - mean(x_clean))^2) / n)^(3/2)
}

if(!is.null(group_var) && length(unique(data_for_analysis[[group_var]])) > 1) {
  for(lipid in lipid_vars) {
    for(grp in unique(data_for_analysis[[group_var]])) {
      values <- data_for_analysis[[lipid]][data_for_analysis[[group_var]] == grp]
      values_clean <- values[!is.na(values)]
      
      # Get best distribution from previous results
      best_dist <- final_results$Distribution[final_results$Variable == lipid & 
                                                final_results$Group == as.character(grp)]
      if(length(best_dist) == 0) best_dist <- "Unknown"
      
      summary_table <- rbind(summary_table, data.frame(
        Variable = lipid,
        Group = as.character(grp),
        N = length(values_clean),
        Mean = round(mean(values_clean), 3),
        Median = round(median(values_clean), 3),
        SD = round(sd(values_clean), 3),
        Min = round(min(values_clean), 3),
        Max = round(max(values_clean), 3),
        Skewness = round(skewness(values_clean), 3),
        Best_Distribution = best_dist[1],
        stringsAsFactors = FALSE
      ))
    }
  }
} else {
  for(lipid in lipid_vars) {
    values <- data_for_analysis[[lipid]]
    values_clean <- values[!is.na(values)]
    
    best_dist <- final_results$Distribution[final_results$Variable == lipid]
    if(length(best_dist) == 0) best_dist <- "Unknown"
    
    summary_table <- rbind(summary_table, data.frame(
      Variable = lipid,
      Group = "Overall",
      N = length(values_clean),
      Mean = round(mean(values_clean), 3),
      Median = round(median(values_clean), 3),
      SD = round(sd(values_clean), 3),
      Min = round(min(values_clean), 3),
      Max = round(max(values_clean), 3),
      Skewness = round(skewness(values_clean), 3),
      Best_Distribution = best_dist[1],
      stringsAsFactors = FALSE
    ))
  }
}

# Print summary table
print("========== COMPREHENSIVE SUMMARY TABLE ==========")
print(summary_table)

# Save summary table
write.csv(summary_table, "descriptive_statistics_by_group.csv", row.names = FALSE)

# Print parameter estimates for each distribution
print("========== DETAILED PARAMETER ESTIMATES ==========")
print(final_results[, c("Variable", "Group", "Distribution", 
                        "Parameter1_Name", "Parameter1", 
                        "Parameter2_Name", "Parameter2", "BIC")])

# ================ FIXED QQ PLOTS ================
# Create QQ plots with margin control
if(length(lipid_vars) <= 4) {
  # Reset graphics
  dev.off()
  par(mar = c(4, 4, 2, 1))
  par(mfrow = c(2, 2))
  
  for(lipid in lipid_vars) {
    tryCatch({
      qqnorm(data_for_analysis[[lipid]], main = paste("Q-Q Plot:", lipid))
      qqline(data_for_analysis[[lipid]], col = "red")
    }, error = function(e) {
      cat("Could not create QQ plot for", lipid, ":", e$message, "\n")
    })
  }
  
  par(mfrow = c(1, 1))
} else {
  cat("Too many variables for QQ plots, skipping...\n")
}

# ================ ALTERNATIVE: SAVE PLOTS TO FILES ================
# If you still have margin issues, save plots directly to files:
cat("\n========== SAVING PLOTS TO FILES ==========\n")

# Save histograms to PDF
pdf("histograms.pdf", width = 8, height = 6)
par(mar = c(4, 4, 2, 1))

if(!is.null(group_var) && length(unique(data_for_analysis[[group_var]])) > 1) {
  for(lipid in lipid_vars) {
    xlim_range <- range(data_for_analysis[[lipid]], na.rm = TRUE)
    groups <- unique(data_for_analysis[[group_var]])
    colors <- rainbow(length(groups))
    
    hist(data_for_analysis[[lipid]][data_for_analysis[[group_var]] == first_group],
         breaks = "Sturges",
         main = paste(lipid, "by", group_var),
         xlab = lipid, col = colors[1], 
         xlim = xlim_range, 
         probability = TRUE,
         border = "white")
    
    # Add other groups
    for(i in 2:length(groups)) {
      group_data <- data_for_analysis[[lipid]][data_for_analysis[[group_var]] == groups[i]]
      group_data <- group_data[!is.na(group_data)]
      
      # Check if data has variation before plotting
      if(length(unique(group_data)) > 1) {
        hist(group_data,
             add = TRUE, col = colors[i], probability = TRUE,
             border = "white", breaks = "FD")
      }
    }
    legend("topright", legend = groups, fill = colors, cex = 0.8)
  }
} else {
  for(lipid in lipid_vars) {
    hist(data_for_analysis[[lipid]], main = lipid, xlab = lipid, 
         col = "lightblue", probability = TRUE)
    lines(density(data_for_analysis[[lipid]], na.rm = TRUE), col = "red", lwd = 2)
  }
}
dev.off()
cat("Histograms saved to histograms.pdf\n")

# Save QQ plots to PDF
pdf("qqplots.pdf", width = 8, height = 6)
par(mar = c(4, 4, 2, 1))
for(lipid in lipid_vars) {
  qqnorm(data_for_analysis[[lipid]], main = paste("Q-Q Plot:", lipid))
  qqline(data_for_analysis[[lipid]], col = "red")
}
dev.off()
cat("QQ plots saved to qqplots.pdf\n")

# ================ FINAL OUTPUT ================
cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Files generated:\n")
cat("1. data_for_analysis.csv - Merged dataset\n")
cat("2. distribution_analysis_results.csv - Distribution parameters by variable\n")
cat("3. descriptive_statistics_by_group.csv - Summary statistics table\n")
cat("4. histograms.pdf - Histogram plots\n")
cat("5. qqplots.pdf - Q-Q plots\n")
