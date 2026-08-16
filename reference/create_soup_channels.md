# Create SoupChannel Objects for Ambient RNA Decontamination

Batch initializes `SoupChannel` objects for multiple samples using lists
of raw (unfiltered) and filtered count matrices. This is a required
preparation step for the `SoupX` pipeline, which uses the empty droplets
in the raw matrix ("Table of Drops" or `tod`) to estimate the ambient
RNA profile contaminating the true cells in the filtered matrix ("Table
of Cells" or `toc`).

## Usage

``` r
create_soup_channels(cts_raw_list, cts_filtered_list, sample_names)
```

## Arguments

- cts_raw_list:

  A named list of raw sparse count matrices (containing all droplets),
  typically from the `raw` element of
  [`load_counts()`](https://chomchonu.github.io/OHMmy/reference/load_counts.md).

- cts_filtered_list:

  A named list of filtered sparse count matrices (containing only called
  cells), typically from the `filtered` element of
  [`load_counts()`](https://chomchonu.github.io/OHMmy/reference/load_counts.md).

- sample_names:

  Character vector. A list of sample identifiers used as keys to match
  the matrices in both lists and name the output list elements.

## Value

A named list of `SoupChannel` objects, ready for downstream
contamination fraction estimation and expression adjustment.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have loaded matrices using the previous pipeline function
counts_list <- load_counts(
  sample_names = my_samples,
  filtered_dir = "data/filtered_outputs",
  raw_dir = "data/raw_outputs"
)

# Initialize the SoupChannels
soup_channels <- create_soup_channels(
  cts_raw_list = counts_list$raw,
  cts_filtered_list = counts_list$filtered,
  sample_names = my_samples
)

# The result is a list ready for autoEstCont() and adjustCounts()
print(soup_channels[["Patient1"]])
} # }
```
