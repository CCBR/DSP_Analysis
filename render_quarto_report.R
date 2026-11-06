# Render the QC or DE Quarto report for a dataset using its YAML configuration.
# Usage: Rscript render_quarto_report.R <qc|de> <config-path>
# Arguments:
#   qc|de       Report to render: quality control (qc) or differential expression (de).
#   config-path Path to the YAML file; it must define paths.output.folder and
#               paths.dataset.name.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("Usage: Rscript render_quarto_report.R <qc|de> <config-path>", call. = FALSE)
}

report_type <- tolower(args[[1L]])
if (!report_type %in% c("qc", "de")) {
  stop("Report type must be 'qc' or 'de'.", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Unable to determine the report helper location.", call. = FALSE)
}
# Resolve paths relative to this script so callers can run it from any directory.
repo_root <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
config_path <- normalizePath(args[[2L]], mustWork = TRUE)

# Use the config to determine where the rendered report should be written.
parameters <- yaml::read_yaml(config_path)
project_folder <- parameters$paths$project.folder
output_folder <- parameters$paths$output.folder
dataset_name <- parameters$paths$dataset.name

# Validate output path components before passing them to Quarto.
if (!is.character(project_folder) || length(project_folder) != 1L ||
    is.na(project_folder) || !nzchar(project_folder)) {
  stop("Config must define a non-empty paths.project.folder.", call. = FALSE)
}
if (!is.character(output_folder) || length(output_folder) != 1L ||
    is.na(output_folder) || !nzchar(output_folder)) {
  stop("Config must define a non-empty paths.output.folder.", call. = FALSE)
}
if (!is.character(dataset_name) || length(dataset_name) != 1L ||
    is.na(dataset_name) || !nzchar(dataset_name) ||
    grepl("[/\\\\]", dataset_name) || dataset_name %in% c(".", "..")) {
  stop("Config must define a valid paths.dataset.name.", call. = FALSE)
}

project_folder <- path.expand(project_folder)
if (!grepl("^/", project_folder)) {
  project_folder <- file.path(repo_root, project_folder)
}
project_folder <- normalizePath(project_folder, winslash = "/", mustWork = FALSE)

output_folder <- path.expand(output_folder)
if (!grepl("^/", output_folder)) {
  output_folder <- file.path(repo_root, output_folder)
}
output_folder <- normalizePath(output_folder, winslash = "/", mustWork = FALSE)
if (identical(project_folder, output_folder) ||
    startsWith(project_folder, paste0(output_folder, "/")) ||
    startsWith(output_folder, paste0(project_folder, "/"))) {
  stop("Report output folder must be separate from the dataset input folder.",
       call. = FALSE)
}
report_output_folder <- file.path(output_folder, "reports", report_type)
if (!dir.exists(report_output_folder) &&
    !dir.create(report_output_folder, recursive = TRUE, showWarnings = FALSE)) {
  stop("Unable to create report output folder: ", report_output_folder, call. = FALSE)
}

# Select the matching report and give the output a dataset-specific name.
report_file <- switch(
  report_type,
  qc = "DSP_Analysis_QC_Report.qmd",
  de = "DSP_DE_Analysis_Report.qmd"
)
output_file <- paste0(dataset_name, "_", toupper(report_type), "_Report.html")

render_args <- c(
  "render",
  file.path(repo_root, report_file),
  "--execute-dir", repo_root,
  "-P", paste0("config.path:", config_path),
  "--output-dir", report_output_folder,
  "--output", output_file
)

# Propagate Quarto failures as R errors so shell callers can stop on failure.
status <- system2("quarto", args = vapply(render_args, shQuote, character(1L)))
if (!identical(status, 0L)) {
  stop("Quarto failed to render ", report_type, " report (exit status ", status, ").",
       call. = FALSE)
}
