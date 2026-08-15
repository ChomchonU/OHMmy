
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OHMmy 🧬

<!-- badges: start -->

[![R-CMD-check](https://github.com/ChomchonU/OHMmy/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChomchonU/OHMmy/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

**OHMmy** is a comprehensive R toolkit engineered for high-resolution
single-cell RNA sequencing (scRNA-seq) analysis. Built natively for
**Seurat v5**, it bridges the gap between raw count matrices and
publication-ready figures.

Designed to handle complex, multi-modal datasets—such as characterizing
highly heterogeneous T cell and NK cell subsets or processing noisy
samples from viral infection models—`OHMmy` provides streamlined,
reproducible wrappers for ambient RNA decontamination, data processing,
functional enrichment, and advanced cell-state visualization.

## 🎯 Core Modules

- **🧼 Ambient RNA Decontamination:** End-to-end integration with
  `SoupX` to accurately model and aggressively remove background mRNA
  contamination prior to clustering.
- **⚙️ Streamlined Seurat v5 Processing:** Automated wrappers for both
  standard Log-normalization and `SCTransform` workflows. Fully
  compatible with Seurat v5 layer structures, minimizing technical noise
  across samples.
- **📊 Advanced Visualizations:**
  - Custom blend plots and multi-gene density plotting (powered by
    `Nebulosa`).
  - Split dotplots and gene-marker heatmaps integrated with hierarchical
    dendrograms for complex sub-clustering.
  - Automated statistical comparisons of cell abundances and metadata
    distributions across experimental groups (`rstatix` & `ggpubr`
    integration).
- **🔬 Functional & Enrichment Profiling:** Built-in pipelines for
  global Differential Expression (DE) screening, Over-Representation
  Analysis (ORA), and Gene Set Enrichment Analysis (GSEA).

## 💻 Installation

You can install the development version of `OHMmy` directly from GitHub
using `pak` or `remotes`.

*(Note: `OHMmy` relies on a modern scRNA-seq stack including
Bioconductor packages. We recommend using `renv` if deploying on
containerized HPC environments to manage dependencies).*

``` r

# Using pak (Recommended)
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("ChomchonU/OHMmy")

# OR using remotes
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("ChomchonU/OHMmy")
```

## 🚀 Quick Start Workflow

`OHMmy` functions are designed to chain together intuitively. Here is an
example of a standard preprocessing and visualization workflow:

``` r

library(Seurat)
library(OHMmy)

# ---------------------------------------------------------
# 1. Initialize & Clean Ambient RNA (SoupX Wrappers)
# ---------------------------------------------------------
soup_channels <- create_soup_channels(raw_counts_list, filtered_counts_list)
soup_channels <- estimate_contamination(soup_channels, initial_clusters)
clean_seurat <- adjustCounts(soup_channels[[1]])

# ---------------------------------------------------------
# 2. Process & Cluster (Seurat v5 Wrappers)
# ---------------------------------------------------------
# Run automated SCTransform workflow
seurat_obj <- ProcessSeuratSCT(clean_seurat)

# Compute PCA, UMAP, and find clusters
seurat_obj <- ClusterAndUMAP(seurat_obj, resolution = 0.5)

# ---------------------------------------------------------
# 3. Analyze & Visualize
# ---------------------------------------------------------
# Statistically compare cell type abundance across infection groups
plot_cell_abundance(seurat_obj, group.by = "Condition", test = "wilcox")

# Generate advanced feature blend plots for subset markers (e.g., CD8A / XCL1)
plot_blend_nebulosa(seurat_obj, features = c("CD8A", "XCL1"))
```

## 📦 Package Structure

- `R/`: Core functional scripts categorized by module (`Processing`,
  `Technical_noise_screening`, `Marker_visualizations`,
  `GSEA_ORA_DEseq2`, etc.)

- `man/`: Roxygen2 documentation files for all exported functions.

- `tests/`: Integrated `testthat` suites ensuring pipeline
  reproducibility and matrix integrity.

## 📄 License

This project is licensed under the MIT License - see the
[LICENSE](https://www.google.com/search?q=LICENSE) file for details.
