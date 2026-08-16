# Extract and Bin Average Gene Expression into Terciles

Leverages the internal data extraction of Seurat's `DotPlot` function to
calculate the unscaled average expression and percent expressed for a
specified list of genes across cell groups. It then evaluates each gene
independently, calculates its expression terciles (33rd and 67th
percentiles), and categorizes each cluster's expression into "Low",
"Int" (Intermediate), or "High" bins. This is highly useful for
converting continuous transcriptomic data into discrete categories for
simplified metadata assignment or categorical plotting.

## Usage

``` r
extract_binned_expression(seurat_obj, gene_list, group_col = "seurat_clusters")
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- gene_list:

  Character vector. A list of specific gene names to extract and bin
  (e.g., `c("CD3E", "CD8A")`).

- group_col:

  Character. The metadata column name defining the cell groups or
  clusters to aggregate the expression by. Default is "seurat_clusters".

## Value

A data frame (`tibble`) containing the following columns:

- `Cluster`: The identity class / group.

- `Gene`: The feature name.

- `AvgExpression`: The raw average expression within the group.

- `PctExpress`: The percentage of cells in the group expressing the
  gene.

- `Expression_Level`: An ordered factor (`"Low", "Int", "High"`)
  representing the binned category.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define a few key functional markers
target_genes <- c("CD8A", "GZMB", "PRF1", "IFNG")

# Extract and bin the expression across fine-resolution clusters
binned_data <- extract_binned_expression(
  seurat_obj = my_seurat,
  gene_list = target_genes,
  group_col = "T_Cell_Subsets"
)

# View the resulting data frame
head(binned_data)

# Filter to find which clusters have "High" expression of GZMB
subset(binned_data, Gene == "GZMB" & Expression_Level == "High")
} # }
```
