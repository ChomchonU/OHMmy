test_that("plot_cell_abundance works correctly", {
  library(Seurat)
  library(ggplot2)

  dummy <- pbmc_small
  # Inject mock metadata
  set.seed(42)
  dummy$Sample <- sample(paste0("Patient_", 1:6), ncol(dummy), replace = TRUE)
  dummy$Condition <- ifelse(dummy$Sample %in% c("Patient_1", "Patient_2", "Patient_3"), "Control", "Treated")
  dummy$CellType <- sample(c("T_cell", "B_cell", "Monocyte"), ncol(dummy), replace = TRUE)

  # Suppress warnings about small sample sizes during testing
  suppressWarnings({
    p <- plot_cell_abundance(
      seurat_obj = dummy,
      sample_col = "Sample",
      condition_col = "Condition",
      celltype_col = "CellType",
      global_test = "kruskal.test",
      output_dir = tempdir() # Save to temporary directory so we don't clutter your project
    )
  })

  expect_s3_class(p, "ggplot")
})

test_that("plot_metadata_stats executes both continuous and categorical tests", {
  library(Seurat)

  dummy <- pbmc_small
  set.seed(42)
  dummy$Sample <- sample(paste0("Patient_", 1:10), ncol(dummy), replace = TRUE)
  dummy$Condition <- ifelse(dummy$Sample %in% paste0("Patient_", 1:5), "WT", "KO")
  dummy$Age <- sample(20:60, ncol(dummy), replace = TRUE) # Continuous
  dummy$Sex <- sample(c("M", "F"), ncol(dummy), replace = TRUE) # Categorical

  suppressWarnings({
    res <- plot_metadata_stats(
      seurat_obj = dummy,
      sample_col = "Sample",
      condition_col = "Condition",
      metadata_vars = c("Age", "Sex"),
      output_dir = tempdir()
    )
  })

  # Check structure
  expect_type(res, "list")
  expect_named(res, c("plots", "stats"))

  # Check continuous routing
  expect_s3_class(res$plots[["Age"]], "ggplot")
  expect_s3_class(res$stats[["Age"]]$Global, "data.frame")

  # Check categorical routing
  expect_s3_class(res$plots[["Sex"]], "ggplot")
})
