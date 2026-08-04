# Shared setup, loaded automatically by testthat before any test-*.R file runs.
#
# Regression tests compare current pipeline output files on disk against
# trusted reference copies. Tests do not re-run the QC/DE reports -- run the
# pipeline first (e.g. via run_DSP_human_kidney_test.sh), then run tests.
#
# Reference files are stored gzip-compressed (.csv.gz) to keep them small
# enough for git/GitHub.
#
# To create or update reference files: run generate_references.R. Do this
# whenever a change is expected to alter results, and commit the updated
# reference files alongside the change that caused the difference.

library(testthat)

reference.dir <- testthat::test_path("references")

# Anchors paths to the project root regardless of the working directory
# test_dir() uses while running tests.
project.root <- rprojroot::find_root(rprojroot::has_file("run_test.sh"))

# Current pipeline output files
current.files <- list(
  q3.normalized.counts = file.path(project.root, "test_datasets/Human_Kidney/qc/human_kidney_test_q3normalized_counts.csv"),
  standr.logcounts      = file.path(project.root, "test_datasets/Human_Kidney/qc/StandR/standr_logcounts.csv"),
  de.glom_v_tubule.default = file.path(project.root, "test_datasets/Human_Kidney/de/DKD_Glomerulus_vs_Tubule_default_de.results.csv"),
  de.glom_v_tubule.standr  = file.path(project.root, "test_datasets/Human_Kidney/de/DKD_Glomerulus_vs_Tubule_standr_de.results.csv")
)

# Matching reference files (gzip-compressed)
reference.files <- list(
  q3.normalized.counts     = file.path(reference.dir, "ref_q3normalized_counts.csv.gz"),
  standr.logcounts         = file.path(reference.dir, "ref_standr_logcounts.csv.gz"),
  de.glom_v_tubule.default = file.path(reference.dir, "ref_DKD_Glomerulus_vs_Tubule_default_de.results.csv.gz"),
  de.glom_v_tubule.standr  = file.path(reference.dir, "ref_DKD_Glomerulus_vs_Tubule_standr_de.results.csv.gz")
)

# Decimal places for standard numeric columns (e.g. logfc).
round.digits <- 3

# Columns compared as fixed-significant-figure scientific notation strings
# rather than rounded decimals. P-values can span many orders of magnitude
# (e.g. 8.56e-49), where rounding to a fixed number of decimal places is
# meaningless and comparing raw doubles is sensitive to floating-point
# formatting differences. Add other small-magnitude columns here as needed.
scientific.cols <- c("pval", "padj")

# Significant figures kept when formatting scientific.cols.
scientific.sig.figs <- 3

# Formats a numeric vector as scientific notation with a fixed number of
# significant figures, e.g. format_scientific(8.564e-49, 3) -> "8.56e-49".
format_scientific <- function(x, sig.figs = scientific.sig.figs) {
  formatC(x, format = "e", digits = sig.figs - 1)
}

# Normalizes numeric columns for comparison: scientific.cols become
# fixed-sig-fig strings, all other numeric columns are rounded to
# round.digits. Applied identically when writing reference files and when
# reading files for comparison, so both sides are always normalized the
# same way.
round_numeric_cols <- function(df, digits = round.digits) {
  for (col_name in names(df)) {
    if (col_name %in% scientific.cols && is.numeric(df[[col_name]])) {
      df[[col_name]] <- format_scientific(df[[col_name]])
    } else if (is.numeric(df[[col_name]])) {
      df[[col_name]] <- round(df[[col_name]], digits)
    }
  }
  df
}

# Reads a CSV, decompressing transparently if gzipped, then normalizes
# numeric columns via round_numeric_cols().
read_csv_for_test <- function(path) {
  if (grepl("\\.gz$", path)) {
    # read.csv() opens and closes an unopened connection itself; do not
    # close it again afterward.
    con <- gzfile(path)
    data <- read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  round_numeric_cols(data)
}

# Compares two tables (counts or results) after sorting both by id.col,
# since row order in an exported CSV can vary independently of the
# underlying values. Column order is not normalized, so a reordered column
# will be flagged as a difference. After round_numeric_cols(), columns in
# scientific.cols are character strings and compared for exact equality;
# other numeric columns are compared using `tolerance`.
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
