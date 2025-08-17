###########################################################################
# Title       : Hierarchical Clustering and Circular Dendrogram 
# Author      : Ehtisham Khokhar
# University  : New Mexico State University
# Email       : ehtishamshakeel@gmail.com 
# Date        : 2025
# Purpose     : Perform hierarchical clustering of genotype trait data, 
#               visualize clusters with circular dendrograms, calculate 
#               Within-Cluster Sum of Squares (WCSS), R² values, and cluster means.
###########################################################################

# Set working directory to the path where data and output files will be stored
setwd("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Mechanical Harvest Paper/Phenotype manuscript/single location/PCA")

# Import dataset (adjusted means)
data <- read.csv("mydata.csv")

# Load necessary libraries
library(circlize)      # For circular dendrogram visualization
library(dendextend)    # For customizing dendrogram aesthetics
library(dplyr)         # For data manipulation (e.g., summarizing by clusters)

# Check dataset structure and summary statistics
head(data)             # View first few rows
str(data)              # Display structure (variable types, column names)
summary(data)          # Summary statistics for each variable
options(max.print = 99999)
print(data)            # Print complete dataset

# Standardize data (zero mean, unit variance) for clustering
data_scaled <- scale(data)

# Calculate distance matrix and perform hierarchical clustering
dd <- dist(scale(data[, -1]), method = "euclidean")  # Exclude identifier column
hc <- hclust(dd, method = "ward.D2")                 # Ward's method for clustering
dend <- as.dendrogram(hc)                            # Convert to dendrogram for customization

# Define function to calculate Within-Cluster Sum of Squares (WCSS) for a specified number of clusters
calc_wcss <- function(hc, k, data) {
  clusters <- cutree(hc, k = k)                # Assign observations to clusters
  total_wcss <- 0
  for (cluster_id in unique(clusters)) {       # Loop through each cluster
    cluster_data <- data[clusters == cluster_id, -1]  # Exclude identifier column
    cluster_mean <- colMeans(cluster_data, na.rm = TRUE)
    total_wcss <- total_wcss + sum((cluster_data - cluster_mean)^2, na.rm = TRUE)
  }
  return(total_wcss)
}

# Compute WCSS for cluster numbers 1 to 10 for the elbow method
wcss_values <- sapply(1:10, function(k) calc_wcss(hc, k, data))

# Plot WCSS values to visualize the "elbow" for optimal clusters
plot(1:10, wcss_values, type = "b", pch = 19, col = "red", lwd = 3, frame = FALSE,
     xlab = "Number of Clusters",
     ylab = "Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Optimal Clusters")

# Add grid lines to the plot
grid()

# Define number of clusters and colors
num_clades <- 8
cluster_colors <- c("firebrick4", "forestgreen", "black", "blue", "orange", "purple", "darkgoldenrod", "deeppink")

# Color branches and labels based on clusters
dend <- dend %>%
  color_branches(k = num_clades, col = cluster_colors) %>%
  color_labels(k = num_clades, col = cluster_colors) %>%
  set("branches_lwd", 2)

# Assign genotype names to dendrogram labels
labels(dend) <- data$Treatment[order.dendrogram(dend)]

# Plot circular dendrogram
circlize_dendrogram(dend, dend_track_height = 0.7)

# Assign clusters and add cluster information to data frame
clusters <- cutree(hc, k = num_clades)
data$Cluster <- clusters
print(data)       # Print updated data frame with clusters

# Export clustered data to CSV
write.csv(data, file = "data_cluster.csv", row.names = FALSE)

# Calculate and summarize means for each trait within clusters
cluster_means8 <- data %>%
  group_by(Cluster) %>%
  summarise(across(c("PHT", "PWDT", "HTFPB", "DTFN", "NBB", "DSRG", "DSRR",
                     "DSFG", "DSFR", "PER", "ARA", "WMH", "MAXW", "HMW",
                     "MAXH", "CURH", "TLY", "RED", "GRN", "PODWT"), mean, na.rm = TRUE))

# Print cluster mean values and export to CSV
print(cluster_means8)
write.csv(cluster_means8, file = "cluster_means8.csv", row.names = FALSE)

# Define function to calculate R² for specified cluster count
calc_r_squared <- function(hc, k) {
  clusters <- cutree(hc, k = k)
  
  # Total sum of squares (TSS) for the entire dataset
  overall_mean <- colMeans(data[, -1], na.rm = TRUE)
  tss <- sum(rowSums((data[, -1] - overall_mean) ^ 2, na.rm = TRUE))
  
  # Calculate WCSS for clusters and derive R²
  total_wcss <- 0
  for (cluster_id in unique(clusters)) {
    cluster_data <- data[clusters == cluster_id, -1]
    cluster_mean <- colMeans(cluster_data, na.rm = TRUE)
    total_wcss <- total_wcss + sum(rowSums((cluster_data - cluster_mean) ^ 2, na.rm = TRUE))
  }
  
  r_squared <- 1 - (total_wcss / tss)
  return(r_squared)
}

# Calculate R² values for 1 to 10 clusters
r_squared_values <- sapply(1:10, function(k) calc_r_squared(hc, k))

# Plot R² values for visualizing the quality of clustering solutions
plot(1:10, r_squared_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters",
     ylab = "R² Value",
     main = "R² Values for Different Number of Clusters")

# Save R² values to CSV for record-keeping
r_squared_df <- data.frame(Num_Clusters = 1:10, R_Squared = r_squared_values)
write.csv(r_squared_df, file = "r_squared_values.csv", row.names = FALSE)
