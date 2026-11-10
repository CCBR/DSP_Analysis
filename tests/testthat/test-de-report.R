library(testthat)

de.model.test.env <- new.env(parent = globalenv())
suppressWarnings(sys.source(
  file.path(
    rprojroot::find_root(rprojroot::has_file("run_test.sh")),
    "DSP_functions.R"
  ),
  envir = de.model.test.env
))

test_that("DKD Glomerulus vs Tubule (default method) DE results match reference", {
  
  skip_if_not(file.exists(reference.files$de.glom_v_tubule.default),
              "Reference not found -- run generate_references.R first")
  skip_if_not(file.exists(current.files$de.glom_v_tubule.default),
              paste("Current output file not found:", current.files$de.glom_v_tubule.default))
  
  current <- read_csv_for_test(current.files$de.glom_v_tubule.default)
  ref  <- read_csv_for_test(reference.files$de.glom_v_tubule.default)
  
  # Adjust "gene" to match the actual gene/target ID column name in this file.
  expect_table_equal(current, ref, id.col = "gene")
})

test_that("default DE runner records the model selected by diffExpr", {
  original.diff_expr <- de.model.test.env$diffExpr
  on.exit(de.model.test.env$diffExpr <- original.diff_expr, add = TRUE)

  de.model.test.env$diffExpr <- function(...) {
    cat("Running DE with linear model (no mixed-effect)\n")
    list(results = data.frame(
      Gene = "gene1",
      contrast_logFC = 1,
      contrast_pval = 0.01,
      contrast_adjpval = 0.02
    ))
  }

  fixed.result <- de.model.test.env$run_default_de(
    lmm.input = list(subset.object = NULL),
    comp.parameters = list(analysis.type = "Within Groups")
  )

  expect_true(is.character(fixed.result$model))
  if (is.character(fixed.result$model)) {
    expect_match(fixed.result$model, "Linear model")
    expect_match(fixed.result$model, "~ sub.group")
  }

  de.model.test.env$diffExpr <- function(...) {
    cat("Running DE with linear mixed-effect model\n")
    list(results = data.frame(
      Gene = "gene1",
      contrast_logFC = 1,
      contrast_pval = 0.01,
      contrast_adjpval = 0.02
    ))
  }

  mixed.result <- de.model.test.env$run_default_de(
    lmm.input = list(subset.object = NULL),
    comp.parameters = list(analysis.type = "Between Groups")
  )

  expect_true(is.character(mixed.result$model))
  if (is.character(mixed.result$model)) {
    expect_match(mixed.result$model, "mixed-effects model")
    expect_match(mixed.result$model, "~ main.group \\+ \\(1 \\| slide\\)")
  }

  within.mixed.result <- de.model.test.env$run_default_de(
    lmm.input = list(subset.object = NULL),
    comp.parameters = list(analysis.type = "Within Groups")
  )

  expect_match(
    within.mixed.result$model,
    "~ sub.group \\+ \\(1 \\+ sub.group \\| slide\\)"
  )
})

test_that("StandR model description uses the design and contrast matrices", {
  design <- matrix(
    c(1, 0, 0.5, 0, 1, -0.5),
    nrow = 2,
    dimnames = list(NULL, c("DKD", "normal", "ruv_W1"))
  )
  contrast <- matrix(
    c(1, -1, 0),
    ncol = 1,
    dimnames = list(colnames(design), "DKD_vs_normal")
  )

  expect_true(is.function(de.model.test.env$describe_standr_model))
  if (is.function(de.model.test.env$describe_standr_model)) {
    model <- de.model.test.env$describe_standr_model(
      design = design,
      contrast = contrast,
      block.field = "slide_name"
    )

    expect_match(model, "DKD")
    expect_match(model, "normal")
    expect_match(model, "DKD_vs_normal")
    expect_match(model, "slide_name")

    unblocked.model <- de.model.test.env$describe_standr_model(
      design = design,
      contrast = contrast
    )
    expect_match(unblocked.model, "no blocking")
  }
})

test_that("model descriptions are printed from cached method output", {
  expect_true(is.function(de.model.test.env$print_de_model))
  if (is.function(de.model.test.env$print_de_model)) {
    method.output <- list(model = "captured model metadata")
    printed <- capture.output(
      de.model.test.env$print_de_model(method.output, "standr")
    )

    expect_true(any(grepl("captured model metadata", printed, fixed = TRUE)))
    expect_false(any(grepl("**", printed, fixed = TRUE)))
    expect_false(any(grepl("`", printed, fixed = TRUE)))

    cached.output <- tempfile()
    on.exit(unlink(cached.output), add = TRUE)
    saveRDS(method.output, cached.output)
    loaded.output <- readRDS(cached.output)

    expect_identical(loaded.output$model, method.output$model)
    expect_true(any(grepl(
      "captured model metadata",
      capture.output(de.model.test.env$print_de_model(loaded.output, "standr")),
      fixed = TRUE
    )))
  }
})

test_that("legacy cached results report their model as unavailable", {
  expect_true(is.function(de.model.test.env$print_de_model))
  if (is.function(de.model.test.env$print_de_model)) {
    printed <- capture.output(
      de.model.test.env$print_de_model(list(results = data.frame()), "default")
    )

    expect_true(any(grepl("unavailable", printed, ignore.case = TRUE)))
    expect_true(any(grepl("overwrite.results = TRUE", printed, fixed = TRUE)))
  }
})

test_that("DKD Glomerulus vs Tubule (standr method) DE results match reference", {
  
  skip_if_not(file.exists(reference.files$de.glom_v_tubule.standr),
              "Reference not found -- run generate_references.R first")
  skip_if_not(file.exists(current.files$de.glom_v_tubule.standr),
              paste("Current output file not found:", current.files$de.glom_v_tubule.standr))
  
  current <- read_csv_for_test(current.files$de.glom_v_tubule.standr)
  ref  <- read_csv_for_test(reference.files$de.glom_v_tubule.standr)
  
  expect_table_equal(current, ref, id.col = "gene")
})