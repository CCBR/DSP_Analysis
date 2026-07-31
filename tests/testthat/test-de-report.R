library(testthat)

test_that("DKD Glomerulus vs Tubule (default method) DE results match reference", {
  
  skip_if_not(file.exists(ref.files$de.glom_v_tubule.default),
              "Reference not found -- run generate_references.R first")
  skip_if_not(file.exists(current.files$de.glom_v_tubule.default),
              paste("Current output file not found:", current.files$de.glom_v_tubule.default))
  
  current <- read_csv_for_test(current.files$de.glom_v_tubule.default)
  ref  <- read_csv_for_test(ref.files$de.glom_v_tubule.default)
  
  # ADJUST "gene" below to match the actual gene/target ID column name in
  # your DE results CSV.
  expect_table_equal(current, ref, id.col = "gene")
})

test_that("DKD Glomerulus vs Tubule (standr method) DE results match reference", {
  
  skip_if_not(file.exists(ref.files$de.glom_v_tubule.standr),
              "Reference not found -- run generate_references.R first")
  skip_if_not(file.exists(current.files$de.glom_v_tubule.standr),
              paste("Current output file not found:", current.files$de.glom_v_tubule.standr))
  
  current <- read_csv_for_test(current.files$de.glom_v_tubule.standr)
  ref  <- read_csv_for_test(ref.files$de.glom_v_tubule.standr)
  
  expect_table_equal(current, ref, id.col = "gene")
})