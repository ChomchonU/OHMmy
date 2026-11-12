# Diagnostic and Variable Evaluation

## Overview

Before trusting a clustering, check what the principal components
actually capture. Two questions matter:

1.  **How many PCs carry biological signal?**
2.  **Are any PCs driven by technical artefacts** such as mitochondrial
    or ribosomal transcripts, nuclear lncRNAs, sequencing depth or cell
    cycle?

OHMmy provides five diagnostic functions. Each saves its figures to disk
so the same checks can run unattended on an HPC node.

| Function | Question it answers |
|----|----|
| [`generate_dimheatmaps()`](https://chomchonu.github.io/OHMmy/reference/generate_dimheatmaps.md) | Which genes and cells drive each PC? |
| [`plot_vizdimloadings()`](https://chomchonu.github.io/OHMmy/reference/plot_vizdimloadings.md) | What are the gene loadings of each PC? |
| [`plot_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_technical_contribution.md) | What fraction of each PC’s top loadings is technical? |
| [`plot_stacked_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_stacked_technical_contribution.md) | Does that fraction change with loading depth? |
| [`plot_pc_metadata_correlation()`](https://chomchonu.github.io/OHMmy/reference/plot_pc_metadata_correlation.md) | Do PC scores correlate with QC covariates? |

``` r

library(Seurat)
library(OHMmy)
source("ohmmy-vignette-helpers.R")

pbmc <- load_pbmc3k_ohmmy()   # built in "Data Processing, Integration, and Clustering"
out_dir <- file.path(tempdir(), "OHMmy-diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
```

The diagnostics use the **pre-integration PCA** (`pca.log`), which still
carries gene loadings.

## 1. Variance and loadings

### Elbow plot

``` r

ElbowPlot(pbmc, reduction = "pca.log", ndims = 40)
```

![](OHMmy-diagnostics-evaluation_files/figure-html/elbow-1.png)

The standard deviation drops steeply over the first 7 PCs and is
essentially flat after PC 10. The 20 PCs used downstream in the
processing article are a conservative choice for pbmc3k.

### Dimension heatmaps

[`generate_dimheatmaps()`](https://chomchonu.github.io/OHMmy/reference/generate_dimheatmaps.md)
splits the PCs into windows (`pc_windows`) and saves one balanced
`DimHeatmap` per window. Each window is wrapped in error handling, so
one failed plot does not stop a batch run.

``` r

generate_dimheatmaps(
  seurat_obj  = pbmc,
  sample_name = "pbmc3k",
  reduction   = "pca.log",
  pc_windows  = list(1:9, 10:18),
  output_dir  = out_dir,
  nfeatures   = 20,
  cells       = 300,
  width_in    = 10,
  height_in   = 12,
  res         = 100
)
```

![](OHMmy-diagnostics-evaluation_files/figure-html/saved-pbmc3k_DimHeatmap_pca.log_PC1-9_20261112_064340.jpg)

PC 1 separates myeloid genes (*LYZ*, *FCN1*, *CST3*) from lymphoid genes
(*IL7R*, *LTB*). Note that *MALAT1* is among its top lymphoid loadings,
which the technical-contribution checks below quantify. PC 2 captures
the cytotoxic NK/CD8 programme against B-cell genes, PCs 3 and 4
platelets, and PC 9 a mixture of interferon-stimulated and S-phase
genes.

### Loading plots

[`plot_vizdimloadings()`](https://chomchonu.github.io/OHMmy/reference/plot_vizdimloadings.md)
stitches one `VizDimLoadings` panel per PC into a grid and returns the
names of any windows that failed. An empty vector means everything was
plotted.

``` r

failed <- plot_vizdimloadings(
  seurat_obj  = pbmc,
  sample_name = "pbmc3k",
  reduction   = "pca.log",
  pc_windows  = list(1:6),
  output_dir  = out_dir,
  nfeatures   = 15,
  width_in    = 14,
  height_in   = 9,
  res         = 100,
  ncol        = 3
)
```

``` r

failed
#> character(0)
```

![](OHMmy-diagnostics-evaluation_files/figure-html/saved-pbmc3k_VizDimLoadings_PC1-6_20261112_064341.jpg)

## 2. Quantifying technical gene contributions

[`plot_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_technical_contribution.md)
takes the top `n_top_genes / 2` positive and negative loadings of each
PC and computes the **loading-weighted percentage** that comes from
genes matching `technical_keywords`. A dashed line marks the `cutoff`
above which a PC deserves a closer look.

``` r

tech_keys <- c("^MT-", "^RPL", "^RPS", "^IG[HKL]", "MALAT1", "NEAT1", "XIST")
```

``` r

plot_technical_contribution(
  seurat_obj         = pbmc,
  sample_name        = "pbmc3k",
  reduction          = "pca.log",
  technical_keywords = tech_keys,
  max_pcs            = 30,
  n_top_genes        = 200,
  cutoff             = 15,
  output_dir         = out_dir
)
```

![](OHMmy-diagnostics-evaluation_files/figure-html/saved-tech_contrib_pca_log_posneg_split_pbmc3k_20261112_064343.jpg)

In pbmc3k every PC stays well below the 15 % cutoff. The largest value
is a few percent on the negative side of PC 1. This is expected, because
PCA loadings only exist for the variable features, and mitochondrial and
ribosomal genes rarely pass HVG selection after QC filtering. In noisier
data, such as stressed tissue, low-quality droplets or nuclei, a PC far
above the line is a strong hint to regress the corresponding covariate.

The percentage can change sharply with the number of genes considered.
[`plot_stacked_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_stacked_technical_contribution.md)
repeats the calculation over a sequence of depths (`gene_depths`) and
facets the result by PC. Bars above the cutoff are highlighted, so a PC
that is “clean” at the top but technical deeper down is easy to spot.

``` r

plot_stacked_technical_contribution(
  seurat_obj         = pbmc,
  sample_name        = "pbmc3k",
  reduction          = "pca.log",
  technical_keywords = tech_keys,
  gene_depths        = seq(0, 300, by = 50),
  max_pcs            = 20,
  cutoff             = 15,
  output_dir         = out_dir,
  plot_width         = 16,
  plot_height        = 10,
  dpi                = 100
)
```

![](OHMmy-diagnostics-evaluation_files/figure-html/saved-stacked_technical_contrib_posneg_pbmc3k_20261112_064345.jpg)

## 3. Correlating PCs with QC covariates

[`plot_pc_metadata_correlation()`](https://chomchonu.github.io/OHMmy/reference/plot_pc_metadata_correlation.md)
computes Spearman correlations between each PC’s cell embeddings and the
metadata columns in `vars_to_test`, then draws a clustered `pheatmap`
with the coefficients printed in each cell. Expression of single nuclear
transcripts, such as *MALAT1*, can be tested by first copying it into
the metadata.

``` r

pbmc$MALAT1_expr <- FetchData(pbmc, vars = "MALAT1", layer = "data")[, 1]
```

``` r

plot_pc_metadata_correlation(
  seurat_obj     = pbmc,
  sample_name    = "pbmc3k",
  vars_to_test   = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt", "percent.ribo",
                     "S.Score", "G2M.Score", "MALAT1_expr"),
  reduction      = "pca.log",
  n_pcs          = 20,
  output_dir     = out_dir,
  plot_width_in  = 13,
  plot_height_in = 5,
  res_dpi        = 100
)
```

![](OHMmy-diagnostics-evaluation_files/figure-html/saved-corr_pc_pbmc3k_20261112_064347.jpg)

## Interpreting the diagnostics

- A PC whose top loadings are dominated by `MT-`, `RPL`/`RPS` or lncRNA
  genes, and whose scores correlate strongly (\|rho\| \> 0.5) with
  `pct_counts_mt`, `percent.ribo` or `nCount_RNA`, is technical. Regress
  that covariate in `ProcessSeuratLOG(vars_to_regress = ...)`, switch to
  [`ProcessSeuratSCT()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratSCT.md),
  or leave the PC out of `dims`.
- Ribosomal genes often load on lymphocyte-versus-myeloid PCs in PBMCs
  simply because naive T cells are ribosome-rich. Before removing such a
  PC, check whether it also separates known cell types, for example with
  the dimension heatmaps above.
- Correlation with `S.Score` or `G2M.Score` indicates a cycling
  compartment. Regress the cell-cycle difference only when proliferation
  is not of interest.

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] OHMmy_1.0.2.9000   Seurat_6.0.0       SeuratObject_5.4.0 sp_2.2-4          
#> 
#> loaded via a namespace (and not attached):
#>   [1] pbapply_1.7-5       gridExtra_2.3.1     rlang_1.3.0        
#>   [4] magrittr_2.0.5      RcppAnnoy_0.0.23    otel_0.3.0         
#>   [7] matrixStats_1.5.0   ggridges_0.5.7      compiler_4.6.2     
#>  [10] png_0.1-9           systemfonts_1.3.2   vctrs_0.7.3        
#>  [13] reshape2_1.4.5      stringr_1.6.0       crayon_1.5.3       
#>  [16] pkgconfig_2.0.3     fastmap_1.2.0       backports_1.5.1    
#>  [19] labeling_0.4.3      ggraph_2.2.2        rmarkdown_2.32     
#>  [22] SoupX_1.6.2         ragg_1.5.2          purrr_1.2.2        
#>  [25] xfun_0.62           cachem_1.1.0        jsonlite_2.0.0     
#>  [28] progress_1.2.3      tweenr_2.0.3        jpeg_0.1-11        
#>  [31] prettyunits_1.2.0   broom_1.0.13        irlba_2.4.1        
#>  [34] parallel_4.6.2      cluster_2.1.8.3     R6_2.6.1           
#>  [37] ica_1.0-3           bslib_0.12.0        stringi_1.8.9      
#>  [40] RColorBrewer_1.1-3  parallelly_1.50.0   car_3.1-5          
#>  [43] jquerylib_0.1.4     scattermore_1.2     Rcpp_1.1.2         
#>  [46] knitr_1.52          future.apply_1.21.0 sctransform_0.4.3  
#>  [49] Matrix_1.7-6        igraph_2.3.4        tidyselect_1.2.1   
#>  [52] viridis_0.6.5       abind_1.4-8         yaml_2.3.12        
#>  [55] codetools_0.2-20    listenv_1.1.0       lattice_0.23-1     
#>  [58] tibble_3.3.1        plyr_1.8.9          withr_3.0.3        
#>  [61] S7_0.2.2            evaluate_1.0.5      Rtsne_0.17         
#>  [64] future_1.82.0       fastDummies_1.7.6   desc_1.4.3         
#>  [67] polyclip_1.10-7     pillar_1.11.1       ggpubr_1.0.0       
#>  [70] carData_3.0-6       leidenbase_0.1.37   KernSmooth_2.23-27 
#>  [73] generics_0.1.4      RcppHNSW_0.7.0      hms_1.1.4          
#>  [76] ggplot2_4.0.3       scales_1.4.0        globals_0.20.0     
#>  [79] glue_1.8.1          pheatmap_1.0.13     tools_4.6.2        
#>  [82] RSpectra_0.16-2     ggsignif_0.6.4      RANN_2.6.3         
#>  [85] graphlayouts_1.2.5  fs_2.1.0            dotCall64_1.2      
#>  [88] tidygraph_1.3.1     cowplot_1.2.0       grid_4.6.2         
#>  [91] tidyr_1.3.2         patchwork_1.3.2     ggforce_0.5.0      
#>  [94] Formula_1.2-6       cli_3.6.6           textshaping_1.0.5  
#>  [97] spam_2.11-4         viridisLite_0.4.3   dplyr_1.2.1        
#> [100] uwot_0.2.5          gtable_0.3.6        pls_2.9-0          
#> [103] rstatix_1.1.0       sass_0.4.10         digest_0.6.39      
#> [106] progressr_1.1.0     ggrepel_0.9.8       htmlwidgets_1.6.4  
#> [109] farver_2.1.2        memoise_2.0.1       htmltools_0.5.9    
#> [112] pkgdown_2.2.1       lifecycle_1.0.5     httr_1.4.9         
#> [115] MASS_7.3-66
```
