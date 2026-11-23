# Fixed 13-gene expression scoring; base R only.
# Genes are rows, samples are columns, and expression must not be log-transformed.
load_score_model <- function(path = "score_model.rda") {
  e <- new.env(parent = emptyenv())
  loaded <- load(path, envir = e)
  if (!identical(loaded, "score_model"))
    stop("Model file must contain only the object score_model.")
  model <- e$score_model
  validate_score_model(model)
  model
}

validate_score_model <- function(model) {
  required <- c("gene", "symbol", "coefficient", "cgga693_symbol", "cgga325_symbol")
  if (!is.list(model) || !identical(model$schema_version, "1.0") ||
      !is.data.frame(model$features) ||
      !all(required %in% names(model$features)))
    stop("Unsupported score_model structure.")
  f <- model$features
  if (nrow(f) != 13L || anyNA(f$gene) || any(!nzchar(f$gene)) ||
      anyDuplicated(f$gene) || !is.numeric(f$coefficient) ||
      any(!is.finite(f$coefficient)))
    stop("Model must contain 13 unique genes and finite numeric coefficients.")
  if (!identical(model$preprocessing$transformation, "log2(expression+1)"))
    stop("Unsupported expression transformation.")
  invisible(TRUE)
}

read_expression_rda <- function(path) {
  e <- new.env(parent = emptyenv())
  loaded <- load(path, envir = e)
  if (length(loaded) != 1L)
    stop("Expression file must contain exactly one numeric matrix.")
  x <- e[[loaded]]
  if (!is.matrix(x) || !is.numeric(x))
    stop("Expression file must contain a numeric gene-by-sample matrix.")
  x
}

score_signature <- function(expression, platform, model, reference = NULL,
                            details = FALSE, input_scale = "untransformed") {
  validate_score_model(model)
  platforms <- c("TCGA", "CGGA693", "CGGA325")
  if (length(platform) != 1L || is.na(platform) || !platform %in% platforms)
    stop("platform must be TCGA, CGGA693 or CGGA325.")
  if (!identical(input_scale, "untransformed"))
    stop("Only untransformed expression is supported; do not supply log-transformed values.")
  if (!is.matrix(expression) || !is.numeric(expression) ||
      nrow(expression) < 1L || ncol(expression) < 1L)
    stop("expression must be a nonempty numeric gene-by-sample matrix.")
  ids <- rownames(expression)
  samples <- colnames(expression)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)))
    stop("Every expression row needs a gene identifier.")
  if (is.null(samples) || anyNA(samples) || any(!nzchar(samples)) ||
      anyDuplicated(samples))
    stop("Sample column names must be present, nonempty and unique.")
  if (!is.logical(details) || length(details) != 1L || is.na(details))
    stop("details must be TRUE or FALSE.")

  f <- model$features
  if (platform == "TCGA") ids <- sub("[.].*", "", ids)
  map_col <- switch(platform, TCGA = "gene", CGGA693 = "cgga693_symbol",
                    CGGA325 = "cgga325_symbol")
  target <- f[[map_col]]
  if (anyNA(target) || any(!nzchar(target)) || anyDuplicated(target))
    stop("Ambiguous or missing target mapping in model.")
  selected <- integer(nrow(f))
  matches <- integer(nrow(f))
  for (j in seq_len(nrow(f))) {
    idx <- which(ids == target[j])
    if (!length(idx))
      stop("Missing model gene: ", target[j], ". No imputation or zero replacement is performed.")
    if (length(idx) > 1L) {
      if (platform != "TCGA")
        stop("Duplicate model symbol: ", target[j], ". Resolve explicitly before scoring.")
      values <- expression[idx, , drop = FALSE]
      if (any(!is.finite(values)) || any(values < 0))
        stop("Invalid duplicate-gene expression.")
      # Resolve duplicates BEFORE log2, across ALL supplied sample columns.
      chosen <- idx[order(-rowMeans(values), idx)[1L]]
      warning("TCGA duplicate ", target[j], ": retained input row ", chosen,
              " using highest mean untransformed expression.")
    } else chosen <- idx
    selected[j] <- chosen
    matches[j] <- length(idx)
  }
  e <- expression[selected, , drop = FALSE]
  if (any(!is.finite(e)) || any(e < 0))
    stop("Selected model-gene expression must be finite, non-negative and nonmissing.")
  transformed <- log2(e + 1)
  raw <- as.numeric(t(transformed) %*% f$coefficient)
  scores <- data.frame(sample = samples, raw_score = raw, stringsAsFactors = FALSE)
  if (!is.null(reference)) {
    if (length(reference) != 1L || is.na(reference) ||
        !reference %in% platforms || reference != platform)
      stop("reference must explicitly name the matching original cohort: ", platform, ".")
    refs <- model$reference_parameters
    needed <- c("platform", "mean_raw", "sd_raw", "median_raw")
    if (!is.data.frame(refs) || !all(needed %in% names(refs)))
      stop("Missing reference parameters.")
    ref <- refs[refs$platform == reference, , drop = FALSE]
    if (nrow(ref) != 1L || !is.numeric(ref$mean_raw) || !is.numeric(ref$sd_raw) ||
        !is.numeric(ref$median_raw) || !is.finite(ref$mean_raw) ||
        !is.finite(ref$sd_raw) || ref$sd_raw <= 0 || !is.finite(ref$median_raw))
      stop("Invalid or ambiguous reference parameters.")
    scores$score_z <- (raw - ref$mean_raw) / ref$sd_raw
    scores$risk_group <- ifelse(raw >= ref$median_raw, "High", "Low")
  }
  if (!details) return(scores)
  contributions <- do.call(rbind, lapply(seq_len(ncol(e)), function(k)
    data.frame(sample = samples[k], gene = f$gene, symbol = f$symbol,
               input_gene_id = rownames(expression)[selected], expression = e[, k],
               log2_expression_plus1 = transformed[, k], coefficient = f$coefficient,
               contribution = transformed[, k] * f$coefficient, row.names = NULL,
               stringsAsFactors = FALSE)))
  mapping <- data.frame(gene = f$gene, symbol = f$symbol,
                        input_gene_id = rownames(expression)[selected],
                        input_row = selected, n_matching_rows = matches,
                        coefficient = f$coefficient, stringsAsFactors = FALSE)
  list(scores = scores, contributions = contributions, mapping_audit = mapping)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% c(2L, 3L))
    stop("Usage: Rscript score_signature.R expression.rda TCGA|CGGA693|CGGA325 [matching_reference_cohort]")
  file_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1L]
  root <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  model <- load_score_model(file.path(root, "score_model.rda"))
  expression <- read_expression_rda(args[1L])
  reference <- if (length(args) == 3L) args[3L] else NULL
  scores <- score_signature(expression, args[2L], model, reference = reference)
  print(scores, row.names = FALSE, digits = 17)
}
