auc_score <- function(y, probability) {
  if (any(!is.finite(probability)) || !setequal(unique(y), 0:1)) stop("AUC requires finite predictions and both classes.")
  positive <- sum(y == 1)
  negative <- sum(y == 0)
  (sum(rank(probability, ties.method = "average")[y == 1]) - positive * (positive + 1) / 2) / (positive * negative)
}

roc_points <- function(y, probability) {
  order <- order(probability, decreasing = TRUE)
  p <- probability[order]
  target <- y[order]
  ends <- c(which(diff(p) != 0), length(p))
  data.frame(threshold = c(Inf, p[ends]),
             sensitivity = c(0, cumsum(target == 1)[ends] / sum(target == 1)),
             specificity = c(1, 1 - cumsum(target == 0)[ends] / sum(target == 0)))
}

select_threshold <- function(y, probability) {
  curve <- roc_points(y, probability)
  curve <- curve[is.finite(curve$threshold), ]
  # Deterministic tie break: the highest threshold among equal Youden scores.
  curve$threshold[which.max(curve$sensitivity + curve$specificity - 1)]
}

classification_metrics <- function(y, probability, threshold) {
  predicted <- as.integer(probability >= threshold)
  tp <- sum(predicted == 1 & y == 1)
  tn <- sum(predicted == 0 & y == 0)
  fp <- sum(predicted == 1 & y == 0)
  fn <- sum(predicted == 0 & y == 1)
  divide <- function(a, b) if (b > 0) a / b else NA_real_
  data.frame(auc = auc_score(y, probability), accuracy = (tp + tn) / length(y),
             sensitivity = divide(tp, tp + fn), specificity = divide(tn, tn + fp),
             precision = divide(tp, tp + fp), f1 = divide(2 * tp, 2 * tp + fp + fn),
             tp = tp, tn = tn, fp = fp, fn = fn)
}

candidate_models <- function() {
  data.frame(
    model = c("clinical_ridge", "clinical_lasso", "clinical_elastic_net", "genomic_ridge", "genomic_lasso", "genomic_elastic_net", "genomic_weighted_lasso", "combined_weighted_ridge", "combined_weighted_lasso", "combined_weighted_elastic_net"),
    features = c(rep("clinical", 3), rep("genomic", 4), rep("combined", 3)),
    alpha = c(0, 1, 0.5, 0, 1, 0.5, 1, 0, 1, 0.5),
    weighted = c(rep(FALSE, 6), rep(TRUE, 4)), stringsAsFactors = FALSE)
}

lambda_grid <- function(alpha, n) if (alpha == 0) 10^seq(3, -3, length.out = n) else 10^seq(0, -4, length.out = n)

fit_path <- function(x, y, candidate, lambda) {
  weights <- if (candidate$weighted) ifelse(y == 1, 1 / mean(y == 1), 1 / mean(y == 0)) else rep(1, length(y))
  glmnet::glmnet(x, y, family = "binomial", alpha = candidate$alpha, lambda = lambda,
                 weights = weights, standardize = TRUE, maxit = 100000L, thresh = 1e-7)
}

tune_models <- function(data, train, config) {
  candidates <- candidate_models()
  folds <- stratified_folds(data$y[train], config$folds, config$seed + 1L)
  predictions <- lapply(seq_len(nrow(candidates)), function(i) matrix(NA_real_, length(train), config$lambda_count))
  audit <- vector("list", config$folds)
  for (fold in seq_len(config$folds)) {
    fitting <- train[folds != fold]
    validating <- train[folds == fold]
    stopifnot(length(intersect(fitting, validating)) == 0)
    prepared <- prepare_partitions(data, fitting, validating, config$top_genes)
    audit[[fold]] <- data.frame(fold = fold, gene = prepared$state$selected_genes)
    for (i in seq_len(nrow(candidates))) {
      candidate <- candidates[i, ]
      lambda <- lambda_grid(candidate$alpha, config$lambda_count)
      fitted <- fit_path(prepared$train[[candidate$features]], data$y[fitting], candidate, lambda)
      predictions[[i]][folds == fold, ] <- predict(fitted, newx = prepared$assessment[[candidate$features]], s = lambda, type = "response")
    }
    message("Completed training validation fold ", fold, "/", config$folds)
  }
  selected <- vector("list", nrow(candidates))
  grid <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    candidate <- candidates[i, ]
    scores <- apply(predictions[[i]], 2, function(p) auc_score(data$y[train], p))
    best <- which.max(scores)
    lambda <- lambda_grid(candidate$alpha, config$lambda_count)
    threshold <- select_threshold(data$y[train], predictions[[i]][, best])
    selected[[i]] <- cbind(candidate, lambda = lambda[best], validation_auc = scores[best], threshold = threshold)
    grid[[i]] <- data.frame(model = candidate$model, lambda = lambda, validation_auc = scores)
    predictions[[i]] <- data.frame(patient_id = data$ids[train], fold = folds, model = candidate$model,
                                    actual = data$y[train], probability = predictions[[i]][, best])
  }
  list(selected = do.call(rbind, selected), grid = do.call(rbind, grid), oof = do.call(rbind, predictions),
       folds = folds, screening = do.call(rbind, audit))
}
