# Clean and Standardize Seurat Dimensional Reductions

Automatically detects, renames, and rebuilds dimensional reductions
within a Seurat object to ensure they comply with Seurat's strict
internal naming conventions. It converts reduction names containing dots
or underscores (e.g., "pca.log") into standard camelCase (e.g.,
"pcaLog"). Crucially, it extracts the raw matrices and completely
rebuilds the `DimReduc` object to guarantee that the internal key and
the column names of both the embeddings and loadings correctly match the
new name.

## Usage

``` r
CleanSeuratReductions(seurat_obj)
```

## Arguments

- seurat_obj:

  A Seurat object containing one or more dimensional reductions.

## Value

An updated `Seurat` object where all dimensional reductions have been
cleanly renamed and rebuilt. The old misnamed reductions are safely
removed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Inspect current reduction names (e.g., "pca.log", "umap_integration")
Reductions(my_seurat)

# Clean and rebuild all reductions
my_seurat <- CleanSeuratReductions(my_seurat)

# Inspect the new names (e.g., "pcaLog", "umapIntegration")
Reductions(my_seurat)
} # }
```
