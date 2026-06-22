#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
data_for_analysis<-read.csv("data_for_analysis.csv")

#-----------descriptive statistics------------------
summary(data_for_analysis)
#-----------for publication tables-----------------
#---------------Creating a custom table--------------
# Homework: Creating a custom table with descriptive statistics results
#install.packages("gtsummary")
#install.packages(c("cardx", "cards"))
#library(cardx)
#library(gtsummary)

#tbl_summary(data_for_analysis)  # Automatic table
#tbl_summary(data_for_analysis, by = outcome)  # By groups

#--------------Statistical Tests---------------------
value_outcome1<-data_for_analysis[data_for_analysis$outcome=="1",]$lipids1
hist(value_outcome1, col = "lightblue")

qqnorm(value_outcome1, main = "Q-Q Plot")
qqline(value_outcome1, col = "red", lwd = 2)

# Shapiro-Wilk test (for n < 5000)
shapiro.test(value_outcome1)


value_outcome0<-data_for_analysis[data_for_analysis$outcome=="0",]$lipids1
hist(value_outcome0, col = "lightgreen")

qqnorm(value_outcome0, main = "Q-Q Plot")
qqline(value_outcome0, col = "red", lwd = 2)

# Shapiro-Wilk test (for n < 5000)
shapiro.test(value_outcome0)

#-------Levene's Test for Homogeneity of Variance--------------
install.packages("car")
library(car)
str(data_for_analysis)
data_for_analysis$outcome<- as.factor(data_for_analysis$outcome)
car::leveneTest(lipids1 ~ outcome, data = data_for_analysis)
#---------------Application of the Brunner-Munzel test----------
install.packages("lawstat")
library(lawstat)
group1 <- data_for_analysis$lipids1[data_for_analysis$outcome == "0"]
group2 <- data_for_analysis$lipids1[data_for_analysis$outcome == "1"]

brunner.munzel.test(group1, group2)
#-------------comparison of results with other tests--------------
t.test(group1, group2)
wilcox.test(group1, group2)

#----------------------------EDA----------------------------------
install.packages("DataExplorer")
library(DataExplorer)
create_report(data_for_analysis)  # Generates HTML report with graphs and statistics
create_report(
  data = data_for_analysis,
  output_file = "EDA_Report.html",  
  output_dir = getwd(),                
  report_title = "EDA Report"          
)

# Load required libraries
library(car)
library(lawstat)
library(dplyr)
library(knitr)

# Convert outcome to factor
data_for_analysis$outcome <- as.factor(data_for_analysis$outcome)

# List of hormone variables (all hormone columns)
hormone_vars <- c("hormone1", "hormone2", "hormone3", "hormone4", "hormone5", 
                  "hormone6", "hormone7", "hormone8", "hormone10_generated")

# Create descriptive statistics table
descriptives_table <- function(data, vars, group_var) {
  results <- data.frame()
  
  for(var in vars) {
    # Group 0 (outcome = 0)
    group0 <- data[data[[group_var]] == "0", var]
    group0_clean <- group0[!is.na(group0)]
    
    # Group 1 (outcome = 1)
    group1 <- data[data[[group_var]] == "1", var]
    group1_clean <- group1[!is.na(group1)]
    
    # Shapiro-Wilk tests
    shapiro0 <- if(length(group0_clean) >= 3 & length(group0_clean) <= 5000) 
      shapiro.test(group0_clean)$p.value else NA
    shapiro1 <- if(length(group1_clean) >= 3 & length(group1_clean) <= 5000) 
      shapiro.test(group1_clean)$p.value else NA
    
    # Determine distribution type (normal if p > 0.05)
    dist0 <- ifelse(!is.na(shapiro0) & shapiro0 > 0.05, "Normal", "Non-normal")
    dist1 <- ifelse(!is.na(shapiro1) & shapiro1 > 0.05, "Normal", "Non-normal")
    
    # Calculate statistics
    row <- data.frame(
      Hormone = var,
      # Group 0 (outcome=0)
      n0 = length(group0_clean),
      Mean0 = round(mean(group0_clean, na.rm = TRUE), 3),
      SD0 = round(sd(group0_clean, na.rm = TRUE), 3),
      Median0 = round(median(group0_clean, na.rm = TRUE), 3),
      IQR0 = round(IQR(group0_clean, na.rm = TRUE), 3),
      Distribution0 = dist0,
      Shapiro_p0 = round(shapiro0, 4),
      # Group 1 (outcome=1)
      n1 = length(group1_clean),
      Mean1 = round(mean(group1_clean, na.rm = TRUE), 3),
      SD1 = round(sd(group1_clean, na.rm = TRUE), 3),
      Median1 = round(median(group1_clean, na.rm = TRUE), 3),
      IQR1 = round(IQR(group1_clean, na.rm = TRUE), 3),
      Distribution1 = dist1,
      Shapiro_p1 = round(shapiro1, 4)
    )
    results <- rbind(results, row)
  }
  return(results)
}

# Generate table
hormone_desc_table <- descriptives_table(data_for_analysis, hormone_vars, "outcome")
print(hormone_desc_table)

#2
# Perform Levene's test and Shapiro-Wilk test for all hormones
hormone_test_results <- data.frame()

for(var in hormone_vars) {
  # Levene's test for homogeneity of variance
  levene_result <- leveneTest(as.formula(paste(var, "~ outcome")), data = data_for_analysis)
  
  # Shapiro-Wilk test for each group
  group0 <- data_for_analysis[data_for_analysis$outcome == "0", var]
  group1 <- data_for_analysis[data_for_analysis$outcome == "1", var]
  
  shapiro0 <- if(sum(!is.na(group0)) >= 3) shapiro.test(group0)$p.value else NA
  shapiro1 <- if(sum(!is.na(group1)) >= 3) shapiro.test(group1)$p.value else NA
  
  hormone_test_results <- rbind(hormone_test_results, data.frame(
    Hormone = var,
    Levene_F = round(levene_result$`F value`[1], 4),
    Levene_p = round(levene_result$`Pr(>F)`[1], 4),
    Shapiro_p0 = round(shapiro0, 4),
    Shapiro_p1 = round(shapiro1, 4)
  ))
}

print(hormone_test_results)


# Set up plotting parameters
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Function to create histograms and Q-Q plots for each hormone
plot_diagnostics <- function(data, hormone_name) {
  group0 <- data[data$outcome == "0", hormone_name]
  group1 <- data[data$outcome == "1", hormone_name]
  
  # Histogram for outcome=0
  hist(group0, main = paste(hormone_name, "- Outcome=0"), 
       xlab = hormone_name, col = "lightblue", border = "white", 
       breaks = 20, freq = FALSE)
  lines(density(group0, na.rm = TRUE), col = "red", lwd = 2)
  
  # Q-Q plot for outcome=0
  qqnorm(group0, main = paste(hormone_name, "- Q-Q Plot (Outcome=0)"))
  qqline(group0, col = "red", lwd = 2)
  
  # Histogram for outcome=1
  hist(group1, main = paste(hormone_name, "- Outcome=1"), 
       xlab = hormone_name, col = "lightgreen", border = "white", 
       breaks = 20, freq = FALSE)
  lines(density(group1, na.rm = TRUE), col = "red", lwd = 2)
  
  # Q-Q plot for outcome=1
  qqnorm(group1, main = paste(hormone_name, "- Q-Q Plot (Outcome=1)"))
  qqline(group1, col = "red", lwd = 2)
}

# Generate plots for all hormones
for(hormone in hormone_vars) {
  plot_diagnostics(data_for_analysis, hormone)
}



# Function to perform all three tests for a given variable
compare_tests <- function(data, var_name) {
  group0 <- data[data$outcome == "0", var_name]
  group1 <- data[data$outcome == "1", var_name]
  
  # Remove NAs
  group0 <- group0[!is.na(group0)]
  group1 <- group1[!is.na(group1)]
  
  # Brunner-Munzel test (robust to non-normality and heteroscedasticity)
  bm_test <- tryCatch({
    brunner.munzel.test(group0, group1)
  }, error = function(e) list(p.value = NA))
  
  # t-test (parametric, assumes normality)
  t_test <- t.test(group0, group1)
  
  # Wilcoxon rank-sum test (non-parametric)
  wilcox_test <- wilcox.test(group0, group1)
  
  return(data.frame(
    Hormone = var_name,
    Brunner_Munzel_p = round(bm_test$p.value, 5),
    t_test_p = round(t_test$p.value, 5),
    Wilcoxon_p = round(wilcox_test$p.value, 5)
  ))
}

# Perform tests for all hormones
all_test_results <- data.frame()
for(hormone in hormone_vars) {
  all_test_results <- rbind(all_test_results, compare_tests(data_for_analysis, hormone))
}

print(all_test_results)



# Load required libraries
install.packages("corrplot")
library(corrplot)
library(ggplot2)

# Load library
library(corrplot)

# Calculate correlation matrices (Spearman for non-normal data)
data_outcome0 <- data_for_analysis[data_for_analysis$outcome == "0", hormone_vars]
data_outcome1 <- data_for_analysis[data_for_analysis$outcome == "1", hormone_vars]
cor_matrix0 <- cor(data_outcome0, method = "spearman", use = "pairwise.complete.obs")
cor_matrix1 <- cor(data_outcome1, method = "spearman", use = "pairwise.complete.obs")

# Create high-quality PDF
pdf("Correlation_Heatmaps.pdf", width = 14, height = 7)

par(mfrow = c(1, 2), mar = c(1, 1, 4, 1))

# Outcome = 0
corrplot(cor_matrix0, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.9,
         title = "Outcome = 0 (Spearman)",
         mar = c(0, 0, 5, 0),
         diag = FALSE,
         col = colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", 
                                  "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", 
                                  "#F46D43", "#D73027", "#A50026"))(100),
         addCoef.col = "black",
         number.cex = 0.7,
         cl.cex = 0.8)

# Outcome = 1
corrplot(cor_matrix1, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.9,
         title = "Outcome = 1 (Spearman)",
         mar = c(0, 0, 5, 0),
         diag = FALSE,
         col = colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", 
                                  "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", 
                                  "#F46D43", "#D73027", "#A50026"))(100),
         addCoef.col = "black",
         number.cex = 0.7,
         cl.cex = 0.8)

dev.off()

# Print correlation matrices
print(round(cor_matrix0, 3))
print(round(cor_matrix1, 3))

# Combine all results into a comprehensive table
final_results <- merge(hormone_test_results, all_test_results, by = "Hormone")
final_results <- merge(hormone_desc_table[, c("Hormone", "Distribution0", "Distribution1")], 
                       final_results, by = "Hormone")

print(final_results)