library(testthat)

de.model.test.env <- new.env(parent = globalenv())
suppressWarnings(sys.source(
  file.path(
    rprojroot::find_root(rprojroot::has_file("run_test.sh")),
    "DSP_functions.R"
  ),
  envir = de.model.test.env
))

test_that("DE backups preserve CSV and RDS results separately", {
  backup.dir <- tempfile("de-backup-")
  dir.create(backup.dir)
  on.exit(unlink(backup.dir, recursive = TRUE), add = TRUE)

  csv.source <- file.path(backup.dir, "results.csv")
  rds.source <- file.path(backup.dir, "results.RDS")
  csv.data <- data.frame(gene = c("A", "B"), logfc = c(1.2, -0.5))
  rds.data <- list(results = csv.data, model = "test model")
  write.csv(csv.data, csv.source, row.names = FALSE)
  saveRDS(rds.data, rds.source)

  backups <- de.model.test.env$backup_de_results(
    result.file = csv.source,
    result.RDS = rds.source,
    backup.folder = file.path(backup.dir, "backups"),
    backup.prefix = "comparison_default_20260101_120000"
  )

  expect_true(file.exists(backups[["csv"]]))
  expect_true(file.exists(backups[["rds"]]))
  expect_equal(read.csv(backups[["csv"]]), csv.data)
  expect_identical(readRDS(backups[["rds"]]), rds.data)
})

test_that("DE backup failures identify missing source files", {
  expect_error(
    de.model.test.env$backup_de_results(
      result.file = tempfile("missing-csv-"),
      result.RDS = tempfile("missing-rds-"),
      backup.folder = tempfile("de-backup-"),
      backup.prefix = "comparison_default"
    ),
    "source file\\(s\\) do not exist"
  )
})

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

test_that("default DE runner records the mixed model selected by diffExpr", {
  original.diff_expr <- de.model.test.env$diffExpr
  original.metadata <- de.model.test.env$de_model_metadata
  on.exit({
    de.model.test.env$diffExpr <- original.diff_expr
    de.model.test.env$de_model_metadata <- original.metadata
  }, add = TRUE)
  de.model.test.env$de_model_metadata <- function(object) {
    data <- data.frame(slide_name = rep(c("s1", "s2"), each = 4),
                       region = rep(c("a", "b"), 4),
                       group = rep(c("g", "g", "h", "h"), 2))
    data[rep(seq_len(nrow(data)), each = 2), ]
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
    comp.parameters = list(
      analysis.type = "Between Groups", group.col = "group",
      groups = c("g", "h"), region.col = "region", regions = "a"
    )
  )

  expect_true(is.character(mixed.result$model))
  if (is.character(mixed.result$model)) {
    expect_match(mixed.result$model, "mixed-effects model")
    expect_match(mixed.result$model, "~ main.group \\+ \\(1 \\| slide\\)")
  }

  within.mixed.result <- de.model.test.env$run_default_de(
    lmm.input = list(subset.object = NULL),
    comp.parameters = list(
      analysis.type = "Within Groups", group.col = "group", groups = "g",
      region.col = "region", regions = c("a", "b")
    )
  )

  expect_match(
    within.mixed.result$model,
    "~ sub.group \\+ \\(1 \\+ sub.group \\| slide\\)"
  )
})

test_that("default DE runner accepts a custom formula and comparison variable", {
  original.diff_expr <- de.model.test.env$diffExpr
  original.mixed_model_de <- de.model.test.env$mixedModelDE
  original.metadata <- de.model.test.env$de_model_metadata
  on.exit({
    de.model.test.env$diffExpr <- original.diff_expr
    de.model.test.env$mixedModelDE <- original.mixed_model_de
    de.model.test.env$de_model_metadata <- original.metadata
  }, add = TRUE)
  de.model.test.env$de_model_metadata <- function(object) {
    data.frame(treatment = c("A", "B", "A", "B"),
               area = c("x", "x", "y", "y"),
               mouse_id = c("1", "1", "2", "2"))
  }

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

test_that("default DE reports the underlying worker error for both fit paths", {
  skip_if_not(isTRUE(parallel::detectCores() > 1))
  original.diff_expr <- de.model.test.env$diffExpr
  original.mixed_model_de <- de.model.test.env$mixedModelDE
  original.metadata <- de.model.test.env$de_model_metadata
  on.exit({
    de.model.test.env$diffExpr <- original.diff_expr
    de.model.test.env$mixedModelDE <- original.mixed_model_de
    de.model.test.env$de_model_metadata <- original.metadata
  }, add = TRUE)
  de.model.test.env$de_model_metadata <- function(object) {
    data.frame(class = c("tk", "wt", "tk", "wt"),
               mouse_id = c("1", "1", "2", "2"),
               slide_name = rep(c("s1", "s2"), each = 4),
               region = rep(c("a", "b"), 4),
               group = rep("g", 8))
  }

  calls <- integer()
  fail.fit <- function(n.cores) {
    calls <<- c(calls, n.cores)
    if (n.cores > 1) {
      warning("all scheduled cores encountered errors in user code")
      stop("subscript out of bounds")
    }
    stop("number of levels of each grouping factor must be < number of observations (problems: mouse_id)")
  }
  de.model.test.env$diffExpr <- function(..., n.cores, multi.core) {
    fail.fit(n.cores)
  }
  de.model.test.env$mixedModelDE <- function(..., nCores, multiCore) {
    fail.fit(nCores)
  }

  for (custom in c(FALSE, TRUE)) {
    calls <- integer()
    parameters <- list(
      contrast.name = "ca1", group.col = "group", groups = "g",
      region.col = "region", regions = c("a", "b")
    )
    if (custom) {
      parameters$model.formula <- "class + (1 | mouse_id)"
      parameters$group.var <- "class"
    }
    expect_error(
      suppressWarnings(suppressMessages(de.model.test.env$run_default_de(
        lmm.input = list(subset.object = NULL),
        comp.parameters = parameters
      ))),
      "comparison 'ca1'.*mouse_id"
    )
    expect_equal(length(calls), 2)
    expect_gt(calls[1], 1)
    expect_identical(calls[2], 1L)
  }
})

test_that("explicit fixed models fit the requested contrast and adjustment", {
  skip_if_not_installed("GeomxTools")
  loadNamespace("GeomxTools")
  example <- system.file(
    "extdata/DSP_NGS_Example_Data/demoData.rds", package = "GeomxTools"
  )
  skip_if_not(nzchar(example), "GeomxTools example data not installed")
  object <- readRDS(example)[1:3, 1:10]
  metadata <- Biobase::pData(object)
  metadata$class <- factor(rep(c("tk", "wt"), each = 5))
  metadata$age <- c(1, 2, 3, 4, 5, 1.5, 2.5, 3.5, 4.5, 5.5)
  Biobase::pData(object) <- metadata
  counts <- rbind(
    c(4, 5, 4.2, 6, 5.4, 1, 2.5, 1.5, 3, 2.8),
    c(2, 3.1, 4, 3, 4.2, 3, 2.5, 3.7, 4.2, 5),
    c(1, 2, 1.5, 2.8, 2.6, 2, 2.9, 3.5, 4, 4.1)
  )
  rownames(counts) <- rownames(object)
  colnames(counts) <- colnames(object)
  Biobase::assayDataElement(object, "log_q") <- counts
  params <- list(
    contrast.name = "fixed", model.type = "fixed",
    model.formula = "~ class + age", group.var = "class",
    contrast.numerator = "tk", contrast.denominator = "wt"
  )

  result <- de.model.test.env$run_default_de(list(subset.object = object), params)
  reference <- t(apply(counts, 1, function(y) {
    d <- data.frame(y = y, class = relevel(metadata$class, "wt"),
                    age = metadata$age)
    summary(lm(y ~ class + age, data = d))$coefficients["classtk", c(1, 4)]
  }))
  expect_equal(unname(result$results$logfc), unname(reference[, 1]))
  expect_equal(unname(result$results$pval), unname(reference[, 2]))
  expect_equal(unname(result$results$padj),
               unname(p.adjust(reference[, 2], "fdr")))
  expect_match(result$model, "no random effect")
  expect_match(result$model, "tk - wt", fixed = TRUE)
  expect_match(
    paste(capture.output(de.model.test.env$print_de_model(result, "default")),
          collapse = ""),
    "Default DE \\(fixed-effects\\)"
  )

  params$contrast.numerator <- "wt"
  params$contrast.denominator <- "tk"
  reverse <- de.model.test.env$run_default_de(list(subset.object = object), params)
  expect_equal(reverse$results$logfc, -result$results$logfc)
  expect_equal(reverse$results$pval, result$results$pval)

  params$model.formula <- "~ class * age"
  params$contrast.numerator <- "tk"
  params$contrast.denominator <- "wt"
  interaction <- de.model.test.env$run_default_de(
    list(subset.object = object), params
  )
  expected <- t(apply(counts, 1, function(y) {
    d <- data.frame(y = y, class = relevel(metadata$class, "wt"),
                    age = metadata$age)
    fit <- lm(y ~ class * age, data = d)
    contrast <- c(0, 1, 0, mean(d$age))
    estimate <- sum(coef(fit) * contrast)
    se <- sqrt(drop(t(contrast) %*% vcov(fit) %*% contrast))
    c(estimate, 2 * pt(abs(estimate / se), df.residual(fit), lower.tail = FALSE))
  }))
  expect_equal(unname(interaction$results$logfc), unname(expected[, 1]))
  expect_equal(unname(interaction$results$pval), unname(expected[, 2]))

  params$model.formula <- "~ class * scale(age + as.numeric(class))"
  transformed <- de.model.test.env$run_default_de(
    list(subset.object = object), params
  )
  expected <- t(apply(counts, 1, function(y) {
    d <- data.frame(y = y, class = metadata$class, age = metadata$age)
    fit <- lm(y ~ class * scale(age + as.numeric(class)), data = d)
    high <- low <- d
    high$class <- factor("tk", levels = levels(d$class))
    low$class <- factor("wt", levels = levels(d$class))
    terms <- delete.response(terms(fit))
    contrast <- colMeans(model.matrix(terms, high) - model.matrix(terms, low))
    estimate <- sum(coef(fit) * contrast)
    se <- sqrt(drop(t(contrast) %*% vcov(fit) %*% contrast))
    c(estimate, 2 * pt(abs(estimate / se), df.residual(fit), lower.tail = FALSE))
  }))
  expect_equal(unname(transformed$results$logfc), unname(expected[, 1]))
  expect_equal(unname(transformed$results$pval), unname(expected[, 2]))

  params$contrast.denominator <- "absent"
  expect_error(
    de.model.test.env$run_default_de(list(subset.object = object), params),
    "both contrast levels"
  )
  params$contrast.denominator <- "wt"
  params$model.formula <- "~ class + (1 | age)"
  expect_error(
    de.model.test.env$run_default_de(list(subset.object = object), params),
    "random effects.*fixed"
  )
  params$model.formula <- "~ class + duplicate"
  metadata$duplicate <- metadata$class
  Biobase::pData(object) <- metadata
  expect_error(
    de.model.test.env$run_default_de(list(subset.object = object), params),
    "rank-deficient"
  )
  params$model.formula <- "~ class"
  counts[1, 1] <- NA_real_
  Biobase::assayDataElement(object, "log_q") <- counts
  expect_error(
    de.model.test.env$run_default_de(list(subset.object = object), params),
    "non-finite expression"
  )
})

test_that("single-slide default fallback requires explicit fixed-model choice", {
  original.diff_expr <- de.model.test.env$diffExpr
  original.metadata <- de.model.test.env$de_model_metadata
  on.exit({
    de.model.test.env$diffExpr <- original.diff_expr
    de.model.test.env$de_model_metadata <- original.metadata
  }, add = TRUE)
  de.model.test.env$diffExpr <- function(...) stop("must not fit")
  de.model.test.env$de_model_metadata <- function(object) {
    data.frame(slide_name = c(rep(c("s1", "s2"), each = 4), "s1", "s1"),
               region = rep(c("a", "b"), 5),
               group = c(rep("g1", 8), "g2", "g2"))
  }
  expect_error(
    de.model.test.env$run_default_de(
      list(subset.object = NULL),
      list(contrast.name = "example", analysis.type = "Within Groups",
           group.col = "group", groups = c("g1", "g2"),
           region.col = "region", regions = c("a", "b"))
    ),
    "example.*g2.*model.type.*fixed"
  )

  de.model.test.env$diffExpr <- function(...) {
    cat("Running DE with linear model (no mixed-effect)\n")
    list(results = "must not be returned")
  }
  expect_error(
    de.model.test.env$run_default_de(
      list(subset.object = NULL),
      list(contrast.name = "example", analysis.type = "Within Groups",
           group.col = "group", groups = "g1",
           region.col = "region", regions = c("a", "b"))
    ),
    "example.*implicit fixed-effects"
  )

  de.model.test.env$de_model_metadata <- function(object) {
    data.frame(slide_name = rep("s1", 4), region = rep("a", 4),
               group = rep(c("g1", "g2"), each = 2))
  }
  expect_error(
    de.model.test.env$run_default_de(
      list(subset.object = NULL),
      list(contrast.name = "between", analysis.type = "Between Groups",
           group.col = "group", groups = c("g1", "g2"),
           region.col = "region", regions = "a")
    ),
    "between.*slide.*parameters"
  )
})

test_that("unestimable random effects are rejected before fitting", {
  original.mixed_model_de <- de.model.test.env$mixedModelDE
  original.metadata <- de.model.test.env$de_model_metadata
  on.exit({
    de.model.test.env$mixedModelDE <- original.mixed_model_de
    de.model.test.env$de_model_metadata <- original.metadata
  }, add = TRUE)
  de.model.test.env$mixedModelDE <- function(...) stop("must not fit")
  de.model.test.env$de_model_metadata <- function(object) {
    data.frame(class = c("tk", "tk", "wt", NA),
               mouse_id = c("1", "2", "3", "3"))
  }
  expect_error(
    de.model.test.env$run_default_de(
      list(subset.object = NULL),
      list(contrast.name = "ca1", model.formula = "~ class + (1 | mouse_id)",
           group.var = "class")
    ),
    "ca1.*mouse_id.*parameters"
  )
  de.model.test.env$de_model_metadata <- function(object) {
    data.frame(class = c("tk", "wt", "tk", "wt"),
               mouse_id = c("1", "1", "2", "2"))
  }
  expect_error(
    de.model.test.env$run_default_de(
      list(subset.object = NULL),
      list(contrast.name = "slope",
           model.formula = "~ class + (1 + class | mouse_id)",
           group.var = "class")
    ),
    "slope.*random term.*parameters"
  )
})

test_that("fixed-model choices cannot reuse incompatible cached results", {
  parameters <- list(
    contrast.name = "example", model.type = "fixed", model.formula = "~ class",
    group.var = "class", contrast.numerator = "tk", contrast.denominator = "wt"
  )
  fixed.output <- list(
    model = "Custom linear model (no random effect): ~class",
    model.spec = de.model.test.env$fixed_de_model_spec(parameters)
  )
  expect_invisible(de.model.test.env$validate_cached_default_model(
    fixed.output, parameters, NULL, "example"
  ))
  changed <- parameters
  changed$contrast.numerator <- "wt"
  changed$contrast.denominator <- "tk"
  expect_error(de.model.test.env$validate_cached_default_model(
    fixed.output, changed, NULL, "example"
  ), "example.*overwrite.results")
  expect_error(de.model.test.env$validate_cached_default_model(
    list(model = "Linear model (no random effect): ~ sub.group"),
    list(), NULL, "example"
  ), "example.*overwrite.results")
  expect_error(de.model.test.env$validate_cached_default_model(
    list(model = "Custom linear mixed-effects model"), parameters, NULL, "example"
  ), "example.*overwrite.results")
  expect_error(de.model.test.env$validate_cached_default_model(
    fixed.output, list(model.type = "typo"), NULL, "example"
  ), "model.type must be")
  expect_error(de.model.test.env$validate_cached_default_model(
    fixed.output, list(model.type = "fixed"), NULL, "example"
  ), "example.*requires model.formula")
})

test_that("default DE never substitutes diagnostic results for a failed fit", {
  expect_error(
    suppressWarnings(suppressMessages(de.model.test.env$run_default_de_fit(
      parallel.fit = function() {
        warning("scheduled core 1 encountered error in user code")
        list(results = "incomplete")
      },
      serial.fit = function() list(results = "diagnostic only"),
      comparison = "test"
    ))),
    "comparison 'test'.*without reproducing the error"
  )

  expect_error(
    de.model.test.env$run_default_de_fit(
      parallel.fit = function() stop("invalid model"),
      serial.fit = function() stop("unexpected retry"),
      comparison = "test"
    ),
    "comparison 'test' failed: invalid model"
  )

  expect_error(
    suppressWarnings(suppressMessages(de.model.test.env$run_default_de_fit(
      parallel.fit = function() stop("subscript out of bounds"),
      serial.fit = function() {
        stop("number of levels of each grouping factor must be < number of observations (problems: mouse_id)")
      },
      comparison = "test"
    ))),
    "mouse_id\\).*parameters file"
  )
})

test_that("real parallel GeoMx model errors expose the grouping failure", {
  skip_on_os("windows")
  skip_if_not_installed("GeomxTools")
  loadNamespace("GeomxTools")
  example <- system.file(
    "extdata/DSP_NGS_Example_Data/demoData.rds",
    package = "GeomxTools"
  )
  skip_if_not(nzchar(example), "GeomxTools example data not installed")
  object <- readRDS(example)[1:2, 1:10]
  metadata <- Biobase::pData(object)
  metadata$class <- factor(rep(c("tk", "wt"), each = 5))
  metadata$mouse_id <- factor(seq_len(nrow(metadata)))
  Biobase::pData(object) <- metadata
  Biobase::assayDataElement(object, "log_q") <-
    log2(Biobase::assayDataElement(object, "exprs") + 1)

  fit <- function(cores) {
    GeomxTools::mixedModelDE(
      object, elt = "log_q", modelFormula = ~class + (1 | mouse_id),
      groupVar = "class", nCores = cores, multiCore = cores > 1
    )
  }
  expect_error(
    suppressWarnings(suppressMessages(de.model.test.env$run_default_de_fit(
      parallel.fit = function() fit(2),
      serial.fit = function() fit(1),
      comparison = "ca1"
    ))),
    "comparison 'ca1'.*number of levels.*mouse_id"
  )

  metadata$class <- factor(rep(c("tk", "tk", "tk", "wt", "wt"), each = 2))
  metadata$mouse_id <- factor(rep(seq_len(5), each = 2))
  Biobase::pData(object) <- metadata
  counts <- Biobase::assayDataElement(object, "log_q")
  counts[1, ] <- c(1, 1.3, 1.5, 1.8, 2, 2.2, 3, 3.2, 3.5, 3.8)
  counts[2, ] <- NA_real_
  Biobase::assayDataElement(object, "log_q") <- counts
  expect_error(
    suppressWarnings(suppressMessages(de.model.test.env$run_default_de_fit(
      parallel.fit = function() fit(2),
      serial.fit = function() fit(1),
      comparison = "gene-specific"
    ))),
    "comparison 'gene-specific'.*0 \\(non-NA\\) cases"
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

test_that("Venn reporting handles empty DEG lists and missing methods", {
  skip_if_not_installed("ggVennDiagram")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("openxlsx")

  report <- readLines(file.path(project.root, "DSP_DE_Analysis_Report.qmd"))
  report.chunk <- function(label) {
    start <- match(paste0("```{r ", label, "}"), report)
    stopifnot(!is.na(start))
    end <- start + match("```", report[(start + 1L):length(report)])
    parse(text = report[(start + 1L):(end - 1L)])
  }
  env <- new.env(parent = globalenv())
  env$ggVennDiagram <- ggVennDiagram::ggVennDiagram
  env$plot_layout <- patchwork::plot_layout
  env$plot_annotation <- patchwork::plot_annotation
  env$str_to_title <- stringr::str_to_title
  env$createWorkbook <- openxlsx::createWorkbook
  env$addWorksheet <- openxlsx::addWorksheet
  env$writeData <- openxlsx::writeData
  env$saveWorkbook <- openxlsx::saveWorkbook
  eval(report.chunk("Venn Diagram Function"), envir = env)

  empty <- data.frame(gene = "a", logfc = 0, padj = 0.9)
  up <- data.frame(gene = "b", logfc = 1, padj = 0.01)
  both <- data.frame(gene = c("b", "c"), logfc = c(1, -1),
                     padj = c(0.01, 0.01))
  env$de.outputs <- list(
    no_degs = list(standr = list(results = empty), default = list(results = empty)),
    missing_method = list(default = list(results = empty)),
    one_sided = list(standr = list(results = empty), default = list(results = up)),
    populated = list(standr = list(results = up), default = list(results = both))
  )

  output.dir <- tempfile("venn-report-")
  dir.create(output.dir)
  dir.create(file.path(output.dir, "method_comparison"))
  on.exit(unlink(output.dir, recursive = TRUE), add = TRUE)
  env$parameters <- list(
    de = list(common = list(export.venn = TRUE, write.deg.lists = TRUE)),
    paths = list(de.folder = paste0(output.dir, "/"))
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  invisible(capture.output(eval(report.chunk("Create Venn Diagrams"), envir = env)))

  expect_null(env$comparison.venn.list$no_degs$up.venn)
  expect_null(env$comparison.venn.list$no_degs$down.venn)
  expect_equal(nrow(env$comparison.venn.list$no_degs$venn.degs.df), 0)
  expect_null(env$comparison.venn.list$missing_method$venn.degs.df)
  expect_identical(env$comparison.venn.list$missing_method$missing.methods, "StandR")
  expect_false(is.null(env$comparison.venn.list$one_sided$up.venn))
  expect_null(env$comparison.venn.list$one_sided$down.venn)
  expect_equal(
    env$comparison.venn.list$populated$venn.degs.df,
    data.frame(group = c("up_shared", "down_default"), gene = c("b", "c"))
  )

  workbook <- file.path(output.dir, "method_comparison", "deg_list_method_comparison.xlsx")
  expect_setequal(openxlsx::getSheetNames(workbook), c("no_degs", "one_sided", "populated"))
  expect_true(file.exists(file.path(output.dir, "method_comparison", "no_degs_de_venn.png")))
  expect_false(file.exists(file.path(output.dir, "method_comparison", "missing_method_de_venn.png")))

  rendered <- paste(capture.output(
    eval(report.chunk("Display Venn Diagrams"), envir = env)
  ), collapse = "\n")
  expect_match(rendered, "No upregulated DEGs for StandR or Default", fixed = TRUE)
  expect_match(rendered, "No downregulated DEGs for StandR or Default", fixed = TRUE)
  expect_match(rendered, "Missing: StandR", fixed = TRUE)
  expect_match(rendered, "No upregulated DEGs for Default", fixed = TRUE)
  expect_match(rendered, "### Populated", fixed = TRUE)
})