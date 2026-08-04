# Run via:
#   Rscript tests/testthat/generate_references.R
#
# Re-run deliberately whenever a pipeline change is EXPECTED to alter
# results, then commit the updated golden CSVs alongside the change that
# caused the difference.

source("tests/testthat/helper-setup.R")

if (!dir.exists(reference.dir)) {
  dir.create(reference.dir, recursive = TRUE)
}

for (key in names(current.files)) {
  src <- current.files[[key]]
  dest <- reference.files[[key]]
  
  if (!file.exists(src)) {
    warning("Current output file not found, skipping: ", src)
    next
  }
  
  data <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)
  data <- round_numeric_cols(data)
  
  con <- gzfile(dest, "w")
  write.csv(data, con, row.names = FALSE)
  close(con)
  
  original.size.mb <- round(file.size(src) / 1024^2, 1)
  compressed.size.mb <- round(file.size(dest) / 1024^2, 1)
  
  message(
    "Saved reference: ", dest,
    " (", original.size.mb, " MB -> ", compressed.size.mb, " MB)"
  )
}