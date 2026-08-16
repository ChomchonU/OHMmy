
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
# 1. Ambient RNA Decontamination (SoupX)
# ---------------------------------------------------------
# Automatically process raw/filtered matrices and transfer cleaned data
soup_results <- process_soupx_samples(filtered_dir = "data/filtered/", raw_dir = "data/raw/")
seurat_obj <- addSoupXMetaToSeurat(seurat_raw, soup_results$final_seurat)

# ---------------------------------------------------------
# 2. Process, Integrate, & Cluster (Seurat v5)
# ---------------------------------------------------------
# Regress technical noise and integrate batches (e.g., CCA)
seurat_obj <- ProcessSeuratLOG(
  seurat_obj, 
  batch_col = "batch", 
  vars_to_regress = c("pct_counts_mt", "S.Score"),
  integration_method = "CCAIntegration"
)

# Compute UMAP and cluster across multiple resolutions
seurat_obj <- ClusterAndUMAP(seurat_obj, reduction = "integrated.cca.log")$seurat

# ---------------------------------------------------------
# 3. Analyze & Visualize
# ---------------------------------------------------------
# Statistically compare cell type abundance across disease states
plot_cell_abundance(
  seurat_obj, 
  sample_col = "donor_id", 
  condition_col = "disease_state", 
  celltype_col = "seurat_clusters"
)

# Generate advanced feature density blend plots for co-expressing subsets
plot_blend_nebulosa(
  seurat_obj, 
  cluster_col = "seurat_clusters", 
  gene_pairs = list(c("CD8A", "GZMB"))
)
```

## 📚 Documentation & Tutorials

This is just a fraction of what `OHMmy` can do! To see the complete,
end-to-end tutorials covering all functions—including complex
hierarchical dot plots, PCA diagnostics, and functional enrichment—visit
our official package website:

🔗 [**OHMmy Official Documentation
Website**](https://www.google.com/search?q=https://chomchonu.github.io/OHMmy/)

Available Guides:

- [Comprehensive End-to-End
  Workflow](https://www.google.com/search?q=https://chomchonu.github.io/OHMmy/articles/OHMmy-comprehensive.html)

- [Data Processing &
  Clustering](https://www.google.com/search?q=https://chomchonu.github.io/OHMmy/articles/OHMmy-processing-clustering.html)

- [Advanced Visualizations & Marker
  Analysis](https://www.google.com/search?q=https://chomchonu.github.io/OHMmy/articles/OHMmy-advanced-visualization.html)

- [Statistical Abundance &
  Correlations](https://www.google.com/search?q=https://chomchonu.github.io/OHMmy/articles/OHMmy-statistical-analysis.html)

- [Differential Expression & Pathway
  Enrichment](https://www.google.com/search?q=https://chomchonu.github.io/OHMmy/articles/OHMmy-differential-expression.html)

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
