test_that("get_sample_names correctly parses 10x output folders", {
  # Create a fake file structure in the temporary directory
  fake_dir <- file.path(tempdir(), "cellranger_out")
  dir.create(fake_dir, showWarnings = FALSE)

  dir.create(file.path(fake_dir, "Patient1_filtered_feature_bc_matrix"))
  dir.create(file.path(fake_dir, "Patient2_filtered_feature_bc_matrix"))
  dir.create(file.path(fake_dir, "Control_feature_bc_matrix")) # Testing alternate naming

  extracted_names <- get_sample_names(fake_dir)

  # Assertions
  expect_type(extracted_names, "character")
  expect_length(extracted_names, 3)
  expect_true(all(c("Patient1", "Patient2", "Control") %in% extracted_names))

  # Cleanup
  unlink(fake_dir, recursive = TRUE)
})
