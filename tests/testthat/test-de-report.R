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

test_that("default DE runner accepts a custom formula and comparison variable", {
  original.diff_expr <- de.model.test.env$diffExpr
  original.mixed_model_de <- de.model.test.env$mixedModelDE
  on.exit({
    de.model.test.env$diffExpr <- original.diff_expr
    de.model.test.env$mixedModelDE <- original.mixed_model_de
  }, add = TRUE)

  de.model.test.env$diffExpr <- function(...) {
    stop("diffExpr should not run when model.formula is configured.")
  }

  calls <- new.env(parent = emptyenv())
  de.model.test.env$mixedModelDE <- function(object, elt, modelFormula,
                                              groupVar, ...) {
    calls$object <- object
    calls$elt <- elt
    calls$formula <- modelFormula
    calls$group.var <- groupVar
    mixed.results <- matrix(
      vector("list", 4),
      nrow = 2,
      dimnames = list(c("anova", "lsmeans"), c("gene1", "gene2"))
    )
    mixed.results["lsmeans", 1] <- list(matrix(
      c(1.23456, 0.012345),
      nrow = 1,
      dimnames = list("treatmentB - treatmentA", c("Estimate", "Pr(>|t|)"))
    ))
    mixed.results["lsmeans", 2] <- list(matrix(
      c(-0.76543, 0.02055),
      nrow = 1,
      dimnames = list("treatmentB - treatmentA", c("Estimate", "Pr(>|t|)"))
    ))
    mixed.results
  }

  custom.result <- de.model.test.env$run_default_de(
    lmm.input = list(subset.object = "subset-object"),
    comp.parameters = list(
      model.formula = "treatment * area + (1 | mouse_id)",
      group.var = "treatment"
    )
  )

  expect_identical(calls$object, "subset-object")
  expect_identical(calls$elt, "log_q")
  expect_equal(
    paste(deparse(calls$formula), collapse = ""),
    "~treatment * area + (1 | mouse_id)"
  )
  expect_identical(calls$group.var, "treatment")
  expect_equal(custom.result$results$gene, c("gene1", "gene2"))
  expect_equal(custom.result$results$logfc, c(1.235, -0.765))
  expect_equal(custom.result$results$pval, c(0.0123, 0.0205))
  expect_equal(custom.result$results$padj, c(0.0205, 0.0205))
  expect_match(custom.result$model, "treatment \\* area")
  expect_match(custom.result$model, "pairwise contrasts by treatment")
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

test_that("report DEG selection applies cutoffs without changing DE results", {
  results <- data.frame(
    gene = c("up", "down", "padj_boundary", "fc_boundary", "not_significant"),
    logfc = c(1.2, -1.4, 3, 1, 2),
    padj = c(0.01, 0.02, 0.05, 0.01, 0.06)
  )
  original.results <- results

  degs <- de.model.test.env$select_degs_for_report(
    results.df = results,
    fc.limit = 1,
    padj.limit = 0.05
  )

  expect_identical(degs$gene, c("up", "down"))
  expect_identical(results, original.results)
})

test_that("static and interactive volcano builders preserve their intended layers", {
  results <- data.frame(
    gene = c("up", "down", "not_significant"),
    logfc = c(2, -2, 0.25),
    padj = c(0.01, 0.02, 0.8)
  )
  parameters <- list(
    contrast.name = "test_comparison",
    fc.limit = 1,
    padj.limit = 0.05,
    legend.coordinates = NULL
  )

  static.plot <- de.model.test.env$make_de_volcano(
    results.df = results,
    comp.parameters = parameters
  )
  interactive.plot <- de.model.test.env$make_de_volcano_interactive(
    results.df = results,
    comp.parameters = parameters
  )

  expect_true(any(vapply(
    static.plot$layers,
    function(layer) inherits(layer$geom, "GeomTextRepel"),
    logical(1)
  )))
  expect_false(any(vapply(
    interactive.plot$layers,
    function(layer) inherits(layer$geom, "GeomTextRepel"),
    logical(1)
  )))
  expect_identical(
    rlang::as_name(interactive.plot$mapping$text),
    "gene"
  )
  interactive.widget <- plotly::ggplotly(interactive.plot, tooltip = "text")
  hover.genes <- unlist(
    lapply(interactive.widget$x$data, `[[`, "text"),
    use.names = FALSE
  )
  expect_setequal(hover.genes[nzchar(hover.genes)], results$gene)
  expect_identical(
    interactive.plot$data$de_direction,
    c("UP", "DOWN", "NONE")
  )
  expect_identical(
    ggplot2::ggplot_build(interactive.plot)$data[[3]]$y,
    -log10(results$padj)
  )
  expect_identical(
    ggplot2::ggplot_build(interactive.plot)$data[[3]]$colour,
    c("red", "blue", "gray60")
  )
})