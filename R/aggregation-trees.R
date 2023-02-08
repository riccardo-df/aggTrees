#' Aggregation Trees
#'
#' Nonparametric data-driven approach to discovering heterogeneous subgroups in a selection-on-observables framework.
#' The approach constructs a sequence of groupings, one for each level of granularity. Groupings are nested and
#' feature an optimality property. For each grouping, we obtain point estimation and inference about group average
#' treatment effects (GATEs). Additionally, we compute the average characteristics of units in each group.
#'
#' @param y Outcome vector.
#' @param D Treatment vector.
#' @param X Covariate matrix (no intercept).
#' @param honest_frac Fraction of observations to be allocated to honest sample.
#' @param method Either \code{"raw"} or \code{"aipw"}, controls how node predictions are computed.
#' @param scores Optional, vector of scores to be used in computing node predictions. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#' @param cates Estimated CATEs. If not provided by the user, CATEs are estimated internally via a \code{\link[grf]{causal_forest}}.
#' @param is_honest Logical vector denoting which observations belong to the honest sample. Required only if the \code{cates} argument is used.
#' @param ... Further arguments from \code{\link[rpart]{rpart.control}}.
#'
#' @return
#' \code{\link{build_aggtree}} returns an \code{aggTrees} object. \code{\link{analyze_aggtree}} returns the fitted model, as an
#' \code{\link[estimatr]{lm_robust}} object, and the scores (if \code{method == "raw"}, this is \code{NULL}). Additionally,
#' it prints LATEX code in the console if \code{verbose == TRUE} (the default).
#'
#' @examples
#' \donttest{
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Construct sequence of groupings. CATEs estimated internally,
#' groupings <- build_aggtree(y, D, X, method = "aipw")
#'
#' ## We can estimate CATEs and pass them.
#' splits <- sample_split(length(y), training_frac = 0.5)
#' training_idx <- splits$training_idx
#' honest_idx <- splits$honest_idx
#'
#' y_tr <- y[training_idx]
#' D_tr <- D[training_idx]
#' X_tr <- X[training_idx, ]
#'
#' y_hon <- y[honest_idx]
#' D_hon <- D[honest_idx]
#' X_hon <- X[honest_idx, ]
#'
#' library(grf)
#' forest <- causal_forest(X_tr, y_tr, D_tr) # Use training sample.
#' cates <- predict(forest, X)$predictions
#'
#' groupings <- build_aggtree(y, D, X, method = "aipw", cates = cates,
#'                            is_honest = 1:length(y) %in% honest_idx)
#'
#' ## We have compatibility with generic S3-methods.
#' summary(groupings)
#' print(groupings)
#' plot(groupings) # Try also setting 'sequence = TRUE'.
#'
#' ## To predict, do the following.
#' tree <- subtree(groupings$tree, cv = TRUE) # Select by cross-validation.
#' predict(tree, data.frame(X))
#'
#' ## Analyze results with 4 groups.
#' results <- analyze_aggtree(groupings, n_groups = 4, method = "aipw", scores = groupings$scores)
#' summary(results$model)}
#'
#' @details
#' Aggregation trees are a three-step procedure. First, CATEs are estimated using any estimator. Second, a tree is grown
#' to approximate the CATEs. Third, the tree is pruned to derive a nested sequence of optimal groupings, one for each
#' granularity level. For each level of granularity, we can obtain point estimation and inference about GATEs.\cr
#'
#' \code{\link{build_aggtree}} constructs the sequence of groupings and estimate GATEs in each node. GATEs can be estimated
#' in several ways. This is controlled by the \code{method} argument. If \code{method == "raw"}, we compute the difference
#' in mean outcomes between treated and control observations in each node. This is an unbiased estimator in randomized experiment.
#' If \code{method == "aipw"}, we construct doubly-robust scores and average them in each node. This is unbiased also in
#' observational studies. Honest regression forests and 5-fold cross fitting are used to estimate the propensity score and the
#' conditional mean function of the outcome (unless the user specifies the argument \code{scores}).\cr
#'
#' The user can provide a vector of estimated CATEs via the \code{cates} argument. If so, the user needs to specify a logical
#' vector to denote which observations belong to the honest sample. If honesty is not desired, \code{is_honest} must be a
#' vector of \code{FALSE}s. If no vector of CATEs is provided, these are estimated internally via a
#' \code{\link[grf]{causal_forest}}.\cr
#'
#' The tree is grown up to some stopping criteria that can be specified by the user. Please refer to
#' the \code{\link[rpart]{rpart.control}} documentation for this.\cr
#'
#' \code{\link{analyze_aggtree}} takes as input an \code{aggTrees} object constructed by \code{\link{build_aggtree}}. Then, for the
#' desired granularity level, chosen via the \code{n_groups} argument, it provides point estimation and inference about
#' GATEs, together with the average characteristics of the units in each group. As before, the \code{method} argument controls
#' how GATEs are estimated. If \code{method == "raw"}, we estimate via OLS the following linear model:
#'
#' \deqn{Y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th group, and \code{|T|} the
#' number of groups. If the treatment is randomly assigned, one can show that the estimated betas identify the GATEs of
#' each group. Thus, we can interpret the OLS results as usual. However, in observational studies these estimates are biased
#' due to selection into treatment. To get unbiased estimates, we can set \code{method} to \code{"aipw"} to construct
#' doubly-robust scores \code{y_i^*} and use them as a pseudo-outcome in the following regression:
#'
#' \deqn{Y_i^* = \sum_{l = 1}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' This way, we get unbiased GATEs estimates, and we can again interpret OLS results as usual.\cr
#'
#' Regardless of the chosen \code{method}, both functions estimate GATEs using observations in the honest sample. If the honest
#' sample is empty, the same data used to construct the tree are used to estimate GATEs. This is fine for prediction but
#' invalidates the inference obtained by \code{\link{analyze_aggtree}}.
#'
#' @import rpart grf
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @seealso
#' \code{\link{plot.aggTrees}}
#'
#' @export
build_aggtree <- function(y, D, X,
                          honest_frac = 0.5, method = "aipw", scores = NULL,
                          cates = NULL, is_honest = NULL, ...) {
  ## Handling inputs and checks.
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)
  if (honest_frac < 0 | honest_frac >= 1) stop("Invalid 'honest_frac'. This must be in the interval [0, 1).", call. = FALSE)
  if (!is.null(cates)) if(is.null(is_honest)) stop("When 'cates' is used, you need to specify honest observations via the argument 'is_honest'.", call. = FALSE)

  ## If necessary, perform training-honest split.
  if (is.null(cates)) {
    idx <- sample_split(length(y), training_frac = (1 - honest_frac))
    training_idx <- idx$training_idx
    honest_idx <- idx$honest_idx
  } else {
    training_idx <- which(!is_honest)
    if (sum(is_honest) > 0) honest_idx <- which(is_honest) else honest_idx = NULL
  }

  ## Define variables.
  y_tr <- y[training_idx]
  D_tr <- D[training_idx]
  X_tr <- X[training_idx, ]

  y_hon <- y[honest_idx]
  D_hon <- D[honest_idx]
  X_hon <- X[honest_idx, ]

  ## If necessary, estimate the CATEs using training sample (estimation step).
  if (is.null(cates)) {
    forest <- grf::causal_forest(X_tr, y_tr, D_tr, num.trees = 4000, tune.parameters = "all")
    forest_predictions <- stats::predict(forest, X) # Predict on whole sample.
    cates <- forest_predictions$predictions
  }

  ## Grow the tree using training sample (tree-growing step).
  tree <- rpart::rpart(cates ~ ., data = data.frame("cates" = cates[training_idx], X_tr), method = "anova", control = rpart::rpart.control(...), model = TRUE)

  ## If adaptive, replace each node with predictions computed in training sample. Otherwise, honest trees.
  if (honest_frac == 0 | (!is.null(is_honest) & sum(is_honest) == 0)) {
    results <- estimate_rpart(tree, y_tr, D_tr, X_tr, method, scores = scores)
  } else {
    results <- estimate_rpart(tree, y_hon, D_hon, X_hon, method, scores = scores)
  }

  new_tree <- results$tree
  scores <- results$scores

  ## Output.
  if (!is.null(is_honest)) forest <- NULL

  out <- list("tree" = new_tree,
              "forest" = forest,
              "scores" = scores,
              "dta" = data.frame(y, D, X),
              "idx" = list("training_idx" = training_idx, "honest_idx" = honest_idx))
  class(out) <- "aggTrees"
  return(out)
}


#' Analyzing Aggregation Trees
#'
#' @param object An \code{aggTrees} object.
#' @param n_groups Number of desired groups.
#' @param method Either \code{"raw"} or \code{"aipw"}, controls how GATEs are estimated.
#' @param scores Optional, vector of scores to be used in estimating GATEs. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#' @param verbose Logical, whether to print in console.
#'
#' @rdname build_aggtree
#'
#' @export
analyze_aggtree <- function(object, n_groups, method = "aipw", scores = NULL, verbose = TRUE) {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (is.null(object$idx$honest_idx)) warning("Inference is not valid, because the same data have been used to construct the tree and estimate GATEs.")
  if (n_groups <= 1) stop("Invalid 'n_groups'. This must be greater than or equal to 2.", call. = FALSE)
  if (!(verbose %in% c(TRUE, FALSE))) stop("Invalid 'verbose'. This must be either TRUE or FALSE.", call. = FALSE)

  tree <- object$tree

  if (n_groups > get_leaves(tree)) stop("The sequence you provided does not contain any grouping with ", n_groups, " groups.", call. = FALSE)

  ## Select appropriate sample (adaptive/honest) according to the output of build_aggtree.
  if (is.null(object$idx$honest_idx)) {
    y <- object$dta$y
    D <- object$dta$D
    X <- object$dta[, -c(1:2)]
  } else {
    y <- object$dta$y[object$idx$honest_idx]
    D <- object$dta$D[object$idx$honest_idx]
    X <- object$dta[object$idx$honest_idx, -c(1:2)]
  }

  ## Select granularity level.
  groups <- subtree(tree, leaves = n_groups)

  ## GATEs point estimates and standard errors.
  results <- causal_ols_rpart(groups, y, X, D, method = method, scores = scores)

  model <- results$model
  scores <- results$scores

  if (method == "raw") {
    gates_point <- coef(summary(model))[(n_groups+1):(n_groups*2), "Estimate"]
    gates_sd <- coef(summary(model))[(n_groups+1):(n_groups*2), "Std. Error"]
    gates_lower <- coef(summary(model))[(n_groups+1):(n_groups*2), "CI Lower"]
    gates_upper <- coef(summary(model))[(n_groups+1):(n_groups*2), "CI Upper"]
  } else if (method == "aipw") {
    gates_point <- coef(summary(model))[, "Estimate"]
    gates_sd <- coef(summary(model))[, "Std. Error"]
    gates_lower <- coef(summary(model))[, "CI Lower"]
    gates_upper <- coef(summary(model))[, "CI Upper"]
  }

  ## Print table.
  if (verbose) avg_characteristics_rpart(groups, X, gates_point, gates_sd)

  ## Output.
  return(list("model" = model, "scores" = scores))
}
