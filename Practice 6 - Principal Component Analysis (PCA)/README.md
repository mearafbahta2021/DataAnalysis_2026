# Principal Component Analysis of Morphometric Data

## Assignment
Practice 6 — Data Analysis course. PCA and Spearman correlation analysis 
of plant morphometric measurements grouped by collection site.

## Data
- File: data_morphometry.txt
- Format: TXT, tab-separated, Windows-1251 encoding, header = TRUE
- 100 individuals, 6 collection sites (grouping factor: Точка_сбора)
- 11 numeric morphometric traits (shoot height, leaf length/width x2, 
  tepal length/width x2, stamen height, pistil height)

## R version and packages
- R version: [run sessionInfo() in your R console and paste it here]
- Packages: Hmisc, knitr, vegan, factoextra, plotly

## Procedures
1. Range standardisation of traits (decostand, method = "range")
2. Spearman rank correlation matrix with significance testing (rcorr)
3. PCA on standardised data (prcomp)
4. 2D biplots (with/without 95% confidence ellipses) by group (factoextra)
5. Interactive 3D biplot of PC1-PC3 with loading vectors (plotly)

## Files
- DA- Practice 6.R — analysis code
- data_morphometry.txt — source data
- Significant_Spearman_Correlations.csv — correlation table output
- PCA_Scores.csv / PCA_Loadings.csv — PCA numeric output
- PCA_Biplot.png / PCA_Biplot_Ellipses.png / PCA_Variables.png / Scree_Plot.png — static figures
- PCA_3D_Interactive_Biplot.html — interactive 3D PCA biplot
- PCA_Morphometry_Report.docx — full written report
