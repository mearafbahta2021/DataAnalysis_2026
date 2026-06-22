install.packages("vegan")
library(vegan)
#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset----------------------------
data<-read.table("data.txt",header=TRUE,sep="\t",check.names = FALSE)
summary (data)
data
rownames(data)<-data[,1]
data<-data[,-1]
data
# Non-metric multidimensional scaling (NMDS) using Euclidean distances
# Perform NMDS ordination based on Euclidean distances.
# Euclidean distance is suitable for quantitative data but may be less appropriate
# for ecological community data due to sensitivity to double zeros.
# The metaMDS function automatically runs multiple random starts and scales the solution.
ord <- metaMDS(data, distance = "euclidean")

#Visualisation of NMDS results
plot(ord, type = "n")                     # Create empty plot (axes only)
points(ord, disp = "sites", pch = 21, 
       cex = 2.5, lwd = 2.5, col = "red") # Add site points (samples)
text(ord, display = "site", cex = 0.7, 
     col = "red", pos = 3)                         # Label sites with their names


# Cluster analysis with Euclidean distances
# Compute Euclidean distance matrix between sites
d<-vegdist(data,method="euclidean")
# Perform hierarchical agglomerative clustering using average linkage (UPGMA)
fit<-hclust(d, method="average")

# Visualise the dendrogram with labels aligned at the baseline
plot(fit, hang =-1)
# Simple dendrogram plot (default parameters)
plot(fit)


#/////////////////////////////////////////////////////////////
#Non-metric multidimensional scaling (NMDS) with Bray-Curtis dissimilarities
# Bray-Curtis is a standard ecological distance measure that ignores double zeros
# and is robust for abundance data.
# The metaMDS function automatically runs multiple random starts and scales the final solution.
ord <- metaMDS(data, distance = "bray")

#Visualisation of NMDS results
plot(ord, type = "n")
points(ord, disp="sites", pch=21, cex=2.5, lwd=2.5, col = "red")
text(ord, display = "site", cex=0.7, col="red", pos = 3)

# Cluster analysis with Bray-Curtis distances
# Compute Bray-Curtis dissimilarity matrix (standard for ecological community data)
d<-vegdist(data,method="bray")
# Perform hierarchical agglomerative clustering using average linkage (UPGMA)
fit<-hclust(d, method="average")

# Visualise dendrogram with labels aligned at the same horizontal level
plot(fit, hang =-1)
# Simple dendrogram plot (default R style)
plot(fit)


#/////////////////////////////////////////////////////////////
# Non-metric multidimensional scaling (NMDS) with Jaccard distance
ord <- metaMDS(data, distance = "jaccard")


plot(ord, type = "n")
points(ord, disp="sites", pch=21, cex=2.5, lwd=2.5, col = "red")
text(ord, display = "site", cex=0.7, col="red", pos = 3)

# Cluster analysis with Jaccard distance
# Compute Jaccard dissimilarity matrix
d<-vegdist(data,method="jaccard")
# Hierarchical clustering using average linkage (UPGMA)
fit<-hclust(d, method="average")
# Visualise dendrogram with labels aligned
plot(fit, hang =-1)
# Default dendrogram plot
plot(fit)


#/////////////////////////////////////////////////////////////
# TASK 5: Detailed analysis with Bray-Curtis distance
#/////////////////////////////////////////////////////////////

# -----------------------------------------------------------
# STEP 1: Compute NMDS with Bray-Curtis dissimilarity
# -----------------------------------------------------------
set.seed(42)  # for reproducibility
ord <- metaMDS(data, distance = "bray", k = 2, trymax = 100)

# Print stress value
cat("NMDS stress:", ord$stress, "\n")

# Stressplot (Shepard diagram)
stressplot(ord, main = "Shepard / Stress plot (Bray-Curtis NMDS)")

# -----------------------------------------------------------
# STEP 2: Fit significant species vectors with envfit (p <= 0.05)
# -----------------------------------------------------------
fit_sp <- envfit(ord, data, permutations = 999, na.rm = TRUE)
print(fit_sp)

# -----------------------------------------------------------
# STEP 3: UPGMA hierarchical clustering with Bray-Curtis
# -----------------------------------------------------------
d_bray <- vegdist(data, method = "bray")
fit_clust <- hclust(d_bray, method = "average")  # UPGMA

plot(fit_clust, hang = -1,
     main = "UPGMA Dendrogram (Bray-Curtis)",
     xlab = "Sites", sub = "", ylab = "Dissimilarity")

# -----------------------------------------------------------
# STEP 4: Cut dendrogram into 2-3 clusters
# -----------------------------------------------------------
k <- 2
clusters <- cutree(fit_clust, k = k)

cat("\nCluster membership (k =", k, "):\n")
print(clusters)

# Highlight clusters on dendrogram
plot(fit_clust, hang = -1,
     main = paste("UPGMA Dendrogram -", k, "clusters (Bray-Curtis)"),
     xlab = "Sites", sub = "", ylab = "Dissimilarity")
rect.hclust(fit_clust, k = k, border = c("steelblue", "tomato")[1:k])

# -----------------------------------------------------------
# STEP 5: NMDS plot with colours, ellipses, arrows, labels
# -----------------------------------------------------------
cluster_cols <- c("steelblue", "tomato", "forestgreen")
site_cols    <- cluster_cols[clusters]

plot(ord, type = "n",
     main = "NMDS (Bray-Curtis) - clusters, ellipses & species vectors")

# 95% confidence ellipses per cluster
ordiellipse(ord, groups = clusters, kind = "se", conf = 0.95,
            col = cluster_cols, lwd = 2, label = FALSE)

# Site points coloured by cluster
points(ord, display = "sites",
       pch = 21, bg = site_cols, col = "black", cex = 2.5, lwd = 1.5)

# Non-overlapping site labels
orditorp(ord, display = "sites",
         col = site_cols, cex = 0.8, air = 0.9, pch = NA)

# Significant species arrows (p <= 0.05)
plot(fit_sp, p.max = 0.05, col = "darkgreen", cex = 0.7, arrow.mul = 0.8, add = TRUE)

# Legend
legend("topright",
       legend = paste("Cluster", 1:k),
       pt.bg  = cluster_cols[1:k],
       pch = 21, pt.cex = 1.5, bty = "n")

# -----------------------------------------------------------
# STEP 6: PERMANOVA - test whether clusters differ significantly
# -----------------------------------------------------------
perm_result <- adonis2(d_bray ~ clusters, permutations = 999)

cat("\n--- PERMANOVA results ---\n")
print(perm_result)

R2   <- perm_result$R2[1]
pval <- perm_result$`Pr(>F)`[1]

cat(sprintf("\nConclusion: PERMANOVA R2 = %.3f, p = %.3f\n", R2, pval))

if (pval <= 0.05) {
  cat("The clusters differ SIGNIFICANTLY in species composition (p <= 0.05).\n")
} else {
  cat("No significant difference between clusters (p > 0.05).\n")
  cat("This is expected with only 6 sites - insufficient power for PERMANOVA.\n")
}


# Save all Task 5 plots as PNG files

# Stressplot
png("plot1_stressplot.png", width=800, height=600)
stressplot(ord, main = "Shepard / Stress plot (Bray-Curtis NMDS)")
dev.off()

# UPGMA Dendrogram plain
png("plot2_dendrogram.png", width=800, height=600)
plot(fit_clust, hang = -1,
     main = "UPGMA Dendrogram (Bray-Curtis)",
     xlab = "Sites", sub = "", ylab = "Dissimilarity")
dev.off()

# UPGMA Dendrogram with cluster boxes
png("plot3_dendrogram_clusters.png", width=800, height=600)
plot(fit_clust, hang = -1,
     main = "UPGMA Dendrogram - 2 clusters (Bray-Curtis)",
     xlab = "Sites", sub = "", ylab = "Dissimilarity")
rect.hclust(fit_clust, k = k, border = c("steelblue", "tomato")[1:k])
dev.off()

# Combined NMDS plot
png("plot4_nmds_clusters.png", width=800, height=600)
plot(ord, type = "n",
     main = "NMDS (Bray-Curtis) - clusters, ellipses & species vectors")
ordiellipse(ord, groups = clusters, kind = "se", conf = 0.95,
            col = cluster_cols, lwd = 2, label = FALSE)
points(ord, display = "sites",
       pch = 21, bg = site_cols, col = "black", cex = 2.5, lwd = 1.5)
orditorp(ord, display = "sites",
         col = site_cols, cex = 0.8, air = 0.9, pch = NA)
plot(fit_sp, p.max = 0.05, col = "darkgreen", cex = 0.7, arrow.mul = 0.8, add = TRUE)
legend("topright", legend = paste("Cluster", 1:k),
       pt.bg = cluster_cols[1:k], pch = 21, pt.cex = 1.5, bty = "n")
dev.off()