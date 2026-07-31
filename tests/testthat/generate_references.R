# Run via:
#   Rscript tests/testthat/generate_references.R
#
# Re-run deliberately whenever a pipeline change is EXPECTED to alter
# results, then commit the updated golden CSVs alongside the change that
# caused the difference.

source("tests/testthat/helper-setup.R")

if (!dir.exists(ref.dir)) {
  dir.create(ref.dir, recursive = TRUE)
}

for (key in names(current.files)) {
  src <- current.files[[key]]
  dest <- ref.files[[key]]
  
  if (!file.exists(src)) {
    warning("Current output file not found, skipping: ", src)
    next
  }
  
  file.copy(src, dest, overwrite = TRUE)
  message("Saved reference: ", dest)
}