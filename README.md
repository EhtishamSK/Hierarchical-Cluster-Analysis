# Hierarchical Clustering and Circular Dendrogram for Genotype Traits

## Overview
This R script performs **hierarchical clustering** on genotype trait data and visualizes the results using **circular dendrograms**. It is designed for plant breeders and researchers who want to explore patterns of variation, identify clusters of genotypes, and summarize trait values within clusters.

The workflow includes:

- Standardization of data for clustering.
- Calculation of **distance matrices** and **hierarchical clustering** using Ward's method.
- Visualization of clusters with **circular dendrograms**.
- Calculation of **Within-Cluster Sum of Squares (WCSS)** to evaluate cluster compactness.
- Calculation of **R² values** to assess clustering quality.
- Assignment of genotypes to clusters and summarization of **cluster means** for all traits.
- Exporting clustered data and summaries to **CSV files**.

---

## Features

- Handles multiple genotype traits simultaneously.
- Generates **circular dendrograms** colored by cluster.
- Provides **cluster assignments** for each genotype.
- Calculates **WCSS** for elbow method to determine optimal cluster number.
- Calculates **R² values** for clustering quality assessment.
- Summarizes **cluster means** for all traits.
- Exports results to **CSV** for further analysis.

---

## Input Requirements

- A CSV file with columns including:
  - `Treatment` or genotype identifiers
  - Numeric columns for traits
- Data can include adjusted means from previous analyses.

---

## Output

The script produces:

- **Circular dendrogram** visualizations with clusters colored.
- **Cluster assignments** for each genotype.
- **Cluster means** for each trait.
- **WCSS values** for cluster evaluation (elbow method).
- **R² values** for assessing clustering quality.
- Exported CSV files for clustered data and cluster summaries.

---

## Example Plot

The hierarchical clustering and circular dendrogram output is saved as `HCA.png` in the repository.  
![Hierarchical Clustering Dendrogram](HCA.png)

---

## Author

**Ehtisham Khokhar**  
New Mexico State University  
Email: ehtisham.nmsu@example.com

---

## Short Description (for GitHub)

R script for hierarchical clustering of genotype traits with circular dendrograms, cluster summaries, and R² evaluation.
