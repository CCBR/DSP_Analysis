library(testthat)

test_that("Q3 normalized counts match reference", {
  
  skip_if_not(file.exists(reference.files$q3.normalized.counts),
              "Reference not found -- run generate_references.R first")
  skip_if_not(file.exists(current.files$q3.normalized.counts),
              paste("Current output file not found:", current.files$q3.normalized.counts))
  
  current <- read_csv_for_test(current.files$q3.normalized.counts)
  ref  <- read_csv_for_test(reference.files$q3.normalized.counts)
  
  # Adjust "gene" if this file's ID column has a different name
  # (e.g. "TargetName", "RTS_ID").
  expect_table_equal(current, ref, id.col = "gene")
})

test_that("StandR logcounts match reference", {
  
  skip_if_not(file.exists(reference.files$standr.logcounts),
              "Reference not found -- run generate_references.R first")
  skip_if_not(file.exists(current.files$standr.logcounts),
              paste("Current output file not found:", current.files$standr.logcounts))
  
  current <- read_csv_for_test(current.files$standr.logcounts)
  ref  <- read_csv_for_test(reference.files$standr.logcounts)
  
  expect_table_equal(current, ref, id.col = "gene")
})
