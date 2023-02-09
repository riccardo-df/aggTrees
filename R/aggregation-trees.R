#' Aggregation Trees
#'
#' Nonparametric data-driven approach to discovering heterogeneous subgroups in a selection-on-observables framework.
#' The approach constructs a sequence of groupings, one for each level of granularity. Groupings are nested and
#' feature an optimality property. For each grouping, we obtain point estimation and standard errors for group average
#' treatment effects (GATEs). Additionally, we assess whether systematic heterogeneity is found by testing the hypothesis
#' that GATEs are the same across groups and that each GATE is different from the smallest. Finally, we investigate the driving
#' factors of effect heterogeneity by computing the average characteristics of units in each group.
#'
#' @param y Outcome vector.
#' @param D Treatment vector.
#' @param X Covariate matrix (no intercept).
#' @param honest_frac Fraction of observations to be allocated to honest sample.
#' @param method Either \code{"raw"} or \code{"aipw"}, controls how node predictions are computed.
#' @param scores Optional, vector of scores to be used in computing node predictions. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#' @param cates Optional, estimated CATEs. If not provided by the user, CATEs are estimated internally via a \code{\link[grf]{causal_forest}}.
#' @param is_honest Logical vector denoting which observations belong to the honest sample. Required only if the \code{cates} argument is used.
#' @param ... Further arguments from \code{\link[rpart]{rpart.control}}.
#'
#' @return
#' \code{\link{build_aggtree}} returns an \code{aggTrees} object.\cr
#'
#' \code{\link{inference_aggtree}} returns an \code{aggTrees.inference} object, which in turn contains the \code{aggTrees} object used
#' in the call.
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
#' ## Inference with 4 groups.
#' results <- inference_aggtree(groupings, n_groups = 4)
#'
#' summary(results$model_gates) # Coefficient of leafk is GATE in k-th leaf.
#' results$gates_all_equal # We reject the null that all GATEs are equal.
#'
#' summary(results$model_diff) # leafk is difference between k-th GATE and smallest GATE.
#' results$gates_diff_smallest # GATEs are significantly different from GATE in leaf 1.
#' print(results, table = "diff")
#'
#' print(results, table = "avg_char")}
#'
#' @md
#' @details
#' Aggregation trees are a three-step procedure. First, conditional average treatment effects (CATEs) are estimated using any
#' estimator. Second, a tree is grown to approximate the CATEs. Third, the tree is pruned to derive a nested sequence of optimal
#' groupings, one for each granularity level. For each level of granularity, we can obtain point estimation and inference about
#' GATEs.\cr
#'
#' To implement this methodology, the user can rely on two core functions that handle the various steps.\cr
#'
#' ## Constructing the Sequence of Groupings
#' \code{\link{build_aggtree}} constructs the sequence of groupings (i.e., the tree) and estimate GATEs in each node. GATEs can be estimated
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
#' ## GATEs Estimation and Inference
#' \code{\link{inference_aggtree}} takes as input an \code{aggTrees} object constructed by \code{\link{build_aggtree}}. Then, for the
#' desired granularity level, chosen via the \code{n_groups} argument, it provides point estimation and standard errors for
#' GATEs. Additionally, it performs some hypothesis testing to assess whether we find systematic heterogeneity and computes
#' the average characteristics of the units in each group to investigate the driving mechanism.
#'
#' ### Point estimates and standard errors for GATEs
#' GATEs and their standard errors are obtained by fitting an appropriate linear model. If \code{method == "raw"}, we estimate
#' via OLS the following:
#'
#' \deqn{Y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th group, and |T| the
#' number of groups. If the treatment is randomly assigned, one can show that the betas identify the GATE of
#' each group. However, this is not true in observational studies due to selection into treatment. In this case, the user is expected
#' to use \code{method == "aipw"} when calling \code{\link{build_aggtree}}. In this case, \code{\link{inference_aggtree}} uses the scores
#' in the following regression:
#'
#' \deqn{score_i = \sum_{l = 1}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' This way, betas again identify GATEs. Regardless of \code{method}, standard errors are estimated via the Eicker-Huber-White
#' estimator.
#'
#' ### Hypothesis testing
#' \code{\link{inference_aggtree}} performs two types of hypothesis testing. First, it uses standard errors from the above models
#' to test whether GATEs in each leaf are the same by constructing a finite-sample F statistic for carrying out a Wald-test-based
#' comparison between a model and a linearly restricted model. Second, it fits a new linear model to estimate and make inference
#' about the difference between each GATE and the smallest GATE. If \code{method == "raw"}, the following model is used:
#'
#' \deqn{Y_i = \delta + \lambda D_i + \sum_{l = 2}^{|T|} L_{i, l} \gamma_l + \sum_{l = 2}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' where leaves 2, ..., |T| are ordered in increasing order of their GATEs. One can show that each beta_l identifies the difference
#' between the l-th GATE and the smallest GATE. Similarly, if \code{method == "aipw"} we fit the following model:
#'
#' \deqn{score_i = \delta + \sum_{l = 2}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' to obtain the same result. We can then test separately the usual null hypotheses that each beta is zero. We adjust p-values to
#' account for multiple hypothesis testing using the Holm's procedure. As before, standard errors are always estimated via the
#' Eicker-Huber-White estimator.\cr
#'
#' ### Average Characteristics
#' \code{\link{inference_aggtree}} regresses each covariate on a set of dummies denoting group membership. This way, we get the
#' average characteristics of units in each leaf, together with a standard error. Leaves are ordered in increasing order of their
#' predictions (from most negative to most positive). Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' ## Caution on Inference
#' Regardless of the chosen \code{method}, both functions estimate GATEs, linear models, and average characteristics of units in each
#' group using only observations in the honest sample. If the honest sample is empty (this happens because the user either sets
#' \code{honest_frac = 0} or passes a vector of \code{FALSE}s as \code{is_honest} when calling \code{\link{build_aggtree}}), the
#' same data used to construct the tree are used to estimate the above quantities. This is fine for prediction but invalidates
#' inference.
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
#' \code{\link{plot.aggTrees}} \code{\link{print.aggTrees.inference}}
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
              "method" = method,
              "dta" = data.frame(y, D, X),
              "idx" = list("training_idx" = training_idx, "honest_idx" = honest_idx))
  class(out) <- "aggTrees"
  return(out)
}


#' Analyzing Aggregation Trees
#'
#' @param object An \code{aggTrees} object.
#' @param n_groups Number of desired groups.
#'
#' @rdname build_aggtree
#'
#' @import car
#'
#' @export
inference_aggtree <- function(object, n_groups) {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (is.null(object$idx$honest_idx)) warning("Inference is not valid, because the same data have been used to construct the tree and estimate GATEs.")
  if (n_groups <= 1) stop("Invalid 'n_groups'. This must be greater than or equal to 2.", call. = FALSE)

  tree <- object$tree

  if (n_groups > get_leaves(tree)) stop("The sequence you provided does not contain any grouping with ", n_groups, " groups.", call. = FALSE)

  method <- object$method
  scores <- object$scores

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

  ## GATEs point estimates and standard errors, and hypotheses testing.
  results <- causal_ols_rpart(groups, y, D, X, method = method, scores = scores)

  model_gates <- results$model_gates
  model_diff <- results$model_diff
  gates_all_equal <- results$gates_all_equal
  gates_diff_smallest <- results$gates_diff_smallest

  ## Compute average characteristics of units in each leaf.
  avg_characteristics <- avg_characteristics_rpart(groups, X)

  ## Output.
  output <- list("aggTree" = object,
                 "groups" = groups,
                 "model_gates" = model_gates,
                 "model_diff" = model_diff,
                 "gates_all_equal" = gates_all_equal,
                 "gates_diff_smallest" = gates_diff_smallest,
                 "avg_characteristics" = avg_characteristics)
  class(output) <- "aggTrees.inference"
  return(output)
}
