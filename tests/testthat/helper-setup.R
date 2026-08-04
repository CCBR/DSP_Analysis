# Shared setup, loaded automatically by testthat before any test-*.R file runs.
#
# WORKFLOW:
#   1. Run the pipeline (QC + DE reports) to produce output CSVs.
#   2. Run generate_references.R once to copy those CSVs into
#      tests/testthat/references/ as the trusted baseline. Commit these to git.
#   3. Run tests any time after that with testthat::test_dir("tests/testthat").
#      Each test re-reads the current output CSV and compares it to the
#      reference copy.
#   4. If a change is expected to alter results, then re-run
#      generate_references.R and commit
#      that update alongside the code/param change that caused it.

library(testthat)

ref.dir <- testthat::test_path("references")

# Paths to the current pipeline output files being tested
project.root <- rprojroot::find_root(rprojroot::has_file("run_test.sh"))

invisible(current.files <- list(
  q3.normalized.counts = file.path(project.root, "test_datasets/Human_Kidney/qc/human_kidney_test_q3normalized_counts.csv"),
  standr.logcounts      = file.path(project.root, "test_datasets/Human_Kidney/qc/StandR/standr_logcounts.csv"),
  de.glom_v_tubule.default = file.path(project.root, "test_datasets/Human_Kidney/de/DKD_Glomerulus_vs_Tubule_default_de.results.csv"),
  de.glom_v_tubule.standr  = file.path(project.root, "test_datasets/Human_Kidney/de/DKD_Glomerulus_vs_Tubule_standr_de.results.csv")
))

# Matching reference reference filename
invisible(ref.files <- list(
  q3.normalized.counts     = file.path(ref.dir, "ref_q3normalized_counts.csv.gz"),
  standr.logcounts         = file.path(ref.dir, "ref_standr_logcounts.csv.gz"),
  de.glom_v_tubule.default = file.path(ref.dir, "ref_DKD_Glomerulus_vs_Tubule_default_de.results.csv.gz"),
  de.glom_v_tubule.standr  = file.path(ref.dir, "ref_DKD_Glomerulus_vs_Tubule_standr_de.results.csv.gz")
))

round.digits <- 3

# Helper: round all numeric columns in a data frame to round.digits decimal
# place to make file sizes smaller
round_numeric_cols <- function(df, digits = round.digits) {
  numeric.cols <- sapply(df, is.numeric)
  df[numeric.cols] <- lapply(df[numeric.cols], round, digits = digits)
  df
}

# Format columns as strings to allow for compression
format_numeric_cols <- function(df, digits = round.digits) {
  numeric.cols <- sapply(df, is.numeric)
  df[numeric.cols] <- lapply(df[numeric.cols], function(x) sprintf(paste0("%.", digits, "f"), x))
  df
}



# Helper: read a CSV robustly (consistent settings across all comparisons)
read_csv_for_test <- function(path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path)
    on.exit(close(con))
    read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  round_numeric_cols(data)
}

# Helper: compare two data frames representing a counts or results table.
# Sorts both by an ID column first, since row order in an exported CSV can
# shift (e.g. tie-breaking in a sort-by-p-value) without the underlying
# values actually differing.
expect_table_equal <- function(current, ref, id.col, tolerance = 1e-3) {
  expect_true(id.col %in% names(current),
              info = paste0("ID column '", id.col, "' not found in current file"))
  expect_true(id.col %in% names(ref),
              info = paste0("ID column '", id.col, "' not found in reference file"))
  
  current.sorted <- current[order(current[[id.col]]), , drop = FALSE]
  ref.sorted  <- ref[order(ref[[id.col]]), , drop = FALSE]
  
  expect_equal(names(current.sorted), names(ref.sorted),
               info = "Column names differ between current and reference files")
  
  expect_equal(nrow(current.sorted), nrow(ref.sorted),
               info = "Row counts differ between current and reference files")
  
  expect_equal(
    current.sorted,
    ref.sorted,
    tolerance = tolerance,
    check.attributes = FALSE
  )
}
