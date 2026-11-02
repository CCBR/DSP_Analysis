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
repo_root <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
config_path <- normalizePath(args[[2L]], mustWork = TRUE)

parameters <- yaml::read_yaml(config_path)
project_folder <- parameters$paths$project.folder
dataset_name <- parameters$paths$dataset.name

if (!is.character(project_folder) || length(project_folder) != 1L ||
    is.na(project_folder) || !nzchar(project_folder)) {
  stop("Config must define a non-empty paths.project.folder.", call. = FALSE)
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
if (!dir.exists(project_folder) &&
    !dir.create(project_folder, recursive = TRUE, showWarnings = FALSE)) {
  stop("Unable to create project output folder: ", project_folder, call. = FALSE)
}

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
  "--output-dir", project_folder,
  "--output", output_file
)

status <- system2("quarto", args = vapply(render_args, shQuote, character(1L)))
if (!identical(status, 0L)) {
  stop("Quarto failed to render ", report_type, " report (exit status ", status, ").",
       call. = FALSE)
}
