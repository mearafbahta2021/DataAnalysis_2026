# Practice 6: Principal Component Analysis of Morphometric Data
# Dataset: data_morphometry.txt (plant morphometric measurements, 6 collection sites)

# Install packages (run once)
install.packages("Hmisc")
install.packages("knitr")
install.packages("vegan")
install.packages("factoextra")
install.packages("plotly")

# Load required libraries
library(Hmisc)       # for rcorr (Spearman correlation with p-values)
library(knitr)        # for kable (nice tables)
library(vegan)        # for decostand (range standardization)
library(factoextra)   # for fviz_pca_biplot (PCA visualization)
library(plotly)       # for interactive 3D biplot
###############################################################
# Working directory
###############################################################

getwd()



###############################################################
data <- read.table(
  "data_morphometry.txt",
  header = TRUE,
  sep = "\t",
  fileEncoding = "Windows-1251",
  check.names = FALSE
)

###############################################################
# Explore the data
###############################################################

head(data)
str(data)
summary(data)

###############################################################
# Separate grouping factor and numeric variables
###############################################################

data_factor <- data[,1,drop=FALSE]

data_morph <- data[,2:ncol(data)]

###############################################################
# Standardization (Range Normalization)
###############################################################

data_std <- decostand(data_morph,
                      method="range",
                      MARGIN=2)

summary(data_std)

###############################################################
# Spearman Correlation Analysis
###############################################################

data_matrix <- as.matrix(data_std)

rcorr_result <- rcorr(data_matrix,
                      type="spearman")

DD <- rcorr_result$r

DP <- rcorr_result$P

DD[DP>0.05] <- 0

diag(DD) <- 1

kable(DD,
      digits=3,
      caption="Significant Spearman Correlation Coefficients (p < 0.05)")

###############################################################
# Principal Component Analysis
###############################################################

fit <- prcomp(data_std)

###############################################################
# PCA Summary
###############################################################

summary(fit)

###############################################################
# Eigenvalues
###############################################################

eig.val <- get_eigenvalue(fit)

kable(eig.val,
      digits=3,
      caption="Eigenvalues")

###############################################################
# PCA Loadings
###############################################################

loadings <- as.data.frame(fit$rotation)

kable(loadings,
      digits=3,
      caption="Principal Component Loadings")

###############################################################
# PCA Scores
###############################################################

scores <- as.data.frame(fit$x)

scores$Group <- data_factor[,1]

head(scores)

###############################################################
# Scree Plot
###############################################################

fviz_eig(fit,
         addlabels=TRUE)

###############################################################
# Variable Contribution Plot
###############################################################

fviz_pca_var(fit,
             col.var="contrib",
             repel=TRUE)

###############################################################
# PCA Biplot
###############################################################

fviz_pca_biplot(fit,
                habillage=data_factor[,1],
                repel=TRUE)

###############################################################
# PCA Biplot with 95% Confidence Ellipses
###############################################################

fviz_pca_biplot(fit,
                habillage=data_factor[,1],
                addEllipses=TRUE,
                ellipse.level=0.95,
                repel=TRUE)

###############################################################
# Individual Plot
###############################################################

fviz_pca_ind(fit,
             geom="point",
             habillage=data_factor[,1],
             addEllipses=TRUE,
             ellipse.level=0.95)

###############################################################
# Variable Plot
###############################################################

fviz_pca_var(fit,
             repel=TRUE)

###############################################################
# Interactive 3D PCA
###############################################################

plot_ly(scores,
        x=~PC1,
        y=~PC2,
        z=~PC3,
        color=~Group,
        colors="Set1",
        type="scatter3d",
        mode="markers") %>%
  layout(title="Interactive 3D PCA")

###############################################################
# Save PCA Biplot
###############################################################

png("PCA_Biplot.png",
    width=1200,
    height=900)

fviz_pca_biplot(fit,
                habillage=data_factor[,1],
                repel=TRUE)

dev.off()

###############################################################
# Save PCA Biplot with Ellipses
###############################################################

png("PCA_Biplot_Ellipses.png",
    width=1200,
    height=900)

fviz_pca_biplot(fit,
                habillage=data_factor[,1],
                addEllipses=TRUE,
                ellipse.level=0.95,
                repel=TRUE)

dev.off()

###############################################################
# Save Variable Contribution Plot
###############################################################

png("PCA_Variables.png",
    width=1200,
    height=900)

fviz_pca_var(fit,
             col.var="contrib",
             repel=TRUE)

dev.off()

###############################################################
# Save Scree Plot
###############################################################

png("Scree_Plot.png",
    width=1200,
    height=900)

fviz_eig(fit,
         addlabels=TRUE)

dev.off()

###############################################################
# Save Correlation Matrix
###############################################################

write.csv(DD,
          "Significant_Spearman_Correlations.csv")

###############################################################
# Save PCA Scores
###############################################################

write.csv(scores,
          "PCA_Scores.csv",
          row.names=FALSE)

###############################################################
# Save PCA Loadings
###############################################################

write.csv(loadings,
          "PCA_Loadings.csv")

