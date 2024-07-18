#' Aggregation Trees
#'
#' Nonparametric data-driven approach to discovering heterogeneous subgroups in a selection-on-observables framework.
#' The approach constructs a sequence of groupings, one for each level of granularity. Groupings are nested and
#' feature an optimality property. For each grouping, we obtain point estimation and standard errors for the group average
#' treatment effects (GATEs) using debiased machine learning procedures. Additionally, we assess whether systematic heterogeneity
#' is found by testing the hypotheses that the differences in the GATEs across all pairs of groups are zero. Finally, we investigate
#' the driving mechanisms of effect heterogeneity by computing the average characteristics of units in each group.
#'
#' @param Y_tr Outcome vector for training sample.
#' @param Y_hon Outcome vector for honest sample.
#' @param D_tr Treatment vector for training sample.
#' @param D_hon Treatment vector for honest sample.
#' @param X_tr Covariate matrix (no intercept) for training sample.
#' @param X_hon Covariate matrix (no intercept) for honest sample.
#' @param cates_tr Optional, predicted CATEs for training sample. If not provided by the user, CATEs are estimated internally via a \code{\link[grf]{causal_forest}}.
#' @param cates_hon Optional, predicted CATEs for honest sample. If not provided by the user, CATEs are estimated internally via a \code{\link[grf]{causal_forest}}.
#' @param method Either \code{"raw"} or \code{"aipw"}, controls how node predictions are computed.
#' @param scores Optional, vector of scores to be used in computing node predictions. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#' @param ... Further arguments from \code{\link[rpart]{rpart.control}}.
#'
#' @return
#' \code{\link{build_aggtree}} returns an \code{aggTrees} object.\cr
#'
#' \code{\link{inference_aggtree}} returns an \code{aggTrees.inference} object, which in turn contains the \code{aggTrees} object used
#' in the call.
#'
#' @examples
#' \donttest{## Generate data.
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
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Training-honest sample split.
#' honest_frac <- 0.5
#' splits <- sample_split(length(Y), training_frac = (1 - honest_frac))
#' training_idx <- splits$training_idx
#' honest_idx <- splits$honest_idx
#'
#' Y_tr <- Y[training_idx]
#' D_tr <- D[training_idx]
#' X_tr <- X[training_idx, ]
#'
#' Y_hon <- Y[honest_idx]
#' D_hon <- D[honest_idx]
#' X_hon <- X[honest_idx, ]
#'
#' ## Construct sequence of groupings. CATEs estimated internally.
#' groupings <- build_aggtree(Y_tr, D_tr, X_tr, # Training sample.
#'                            Y_hon, D_hon, X_hon) # Honest sample.
#'
#' ## Alternatively, we can estimate the CATEs and pass them.
#' library(grf)
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # Use training sample.
#' cates_tr <- predict(forest, X_tr)$predictions
#' cates_hon <- predict(forest, X_hon)$predictions
#'
#' groupings <- build_aggtree(Y_tr, D_tr, X_tr, # Training sample.
#'                            Y_hon, D_hon, X_hon, # Honest sample.
#'                            cates_tr, cates_hon) # Predicted CATEs.
#'
#' ## We have compatibility with generic S3-methods.
#' summary(groupings)
#' print(groupings)
#' plot(groupings) # Try also setting 'sequence = TRUE'.
#'
#' ## To predict, do the following.
#' tree <- subtree(groupings$tree, cv = TRUE) # Select by cross-validation.
#' head(predict(tree, data.frame(X_hon)))
#'
#' ## Inference with 4 groups.
#' results <- inference_aggtree(groupings, n_groups = 4)
#'
#' summary(results$model) # Coefficient of leafk is GATE in k-th leaf.
#'
#' results$gates_diff_pairs$gates_diff # GATEs differences.
#' results$gates_diff_pairs$holm_pvalues # leaves 1-2 not statistically different.
#'
#' ## LATEX.
#' print(results, table = "diff")
#' print(results, table = "avg_char")}
#'
#' @md
#' @details
#' Aggregation trees are a three-step procedure. First, the conditional average treatment effects (CATEs) are estimated using any
#' estimator. Second, a tree is grown to approximate the CATEs. Third, the tree is pruned to derive a nested sequence of optimal
#' groupings, one for each granularity level. For each level of granularity, we can obtain point estimation and inference about
#' the GATEs.\cr
#'
#' To implement this methodology, the user can rely on two core functions that handle the various steps.\cr
#'
#' ## Constructing the Sequence of Groupings
#' \code{\link{build_aggtree}} constructs the sequence of groupings (i.e., the tree) and estimate the GATEs in each node. The
#' GATEs can be estimated in several ways. This is controlled by the \code{method} argument. If \code{method == "raw"}, we
#' compute the difference in mean outcomes between treated and control observations in each node. This is an unbiased estimator
#' in randomized experiment. If \code{method == "aipw"}, we construct doubly-robust scores and average them in each node. This
#' is unbiased also in observational studies. Honest regression forests and 5-fold cross fitting are used to estimate the
#' propensity score and the conditional mean function of the outcome (unless the user specifies the argument \code{scores}).\cr
#'
#' The user can provide a vector of the estimated CATEs via the \code{cates_tr} and \code{cates_hon} arguments. If no CATEs are provided,
#' these are estimated internally via a \code{\link[grf]{causal_forest}} using only the training sample, that is, \code{Y_tr}, \code{D_tr},
#' and \code{X_tr}.\cr
#'
#' ## GATEs Estimation and Inference
#' \code{\link{inference_aggtree}} takes as input an \code{aggTrees} object constructed by \code{\link{build_aggtree}}. Then, for
#' the desired granularity level, chosen via the \code{n_groups} argument, it provides point estimation and standard errors for
#' the GATEs. Additionally, it performs some hypothesis testing to assess whether we find systematic heterogeneity and computes
#' the average characteristics of the units in each group to investigate the driving mechanisms.
#'
#' ### Point estimates and standard errors for the GATEs
#' GATEs and their standard errors are obtained by fitting an appropriate linear model. If \code{method == "raw"}, we estimate
#' via OLS the following:
#'
#' \deqn{Y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th group, and |T| the
#' number of groups. If the treatment is randomly assigned, one can show that the betas identify the GATE of
#' each group. However, this is not true in observational studies due to selection into treatment. In this case, the user is
#' expected to use \code{method == "aipw"} when calling \code{\link{build_aggtree}}. In this case,
#' \code{\link{inference_aggtree}} uses the scores in the following regression:
#'
#' \deqn{score_i = \sum_{l = 1}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' This way, betas again identify the GATEs.\cr
#'
#' Regardless of \code{method}, standard errors are estimated via the Eicker-Huber-White estimator.\cr
#'
#' If \code{boot_ci == TRUE}, the routine also computes asymmetric bias-corrected and accelerated 95% confidence intervals using 2000 bootstrap
#' samples. Particularly useful when the honest sample is small-ish.
#'
#' ### Hypothesis testing
#' \code{\link{inference_aggtree}} uses the standard errors obtained by fitting the linear models above to test the hypotheses
#' that the GATEs are different across all pairs of leaves. Here, we adjust p-values to account for multiple hypotheses testing
#' using Holm's procedure.
#'
#' ### Average Characteristics
#' \code{\link{inference_aggtree}} regresses each covariate on a set of dummies denoting group membership. This way, we get the
#' average characteristics of units in each leaf, together with a standard error. Leaves are ordered in increasing order of their
#' predictions (from most negative to most positive). Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' ## Caution on Inference
#' Regardless of the chosen \code{method}, both functions estimate the GATEs, the linear models, and the average characteristics
#' of units in each group using only observations in the honest sample. If the honest sample is empty (this happens when the
#' user either does not provide \code{Y_hon}, \code{D_hon}, and \code{X_hon} or sets them to \code{NULL}), the same data used to
#' construct the tree are used to estimate the above quantities. This is fine for prediction but invalidates inference.
#'
#' @import rpart grf
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R. (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso \code{\link{plot.aggTrees}} \code{\link{print.aggTrees.inference}}
#'
#' @export
build_aggtree <- function(Y_tr, D_tr, X_tr,
                          Y_hon = NULL, D_hon = NULL, X_hon = NULL,
                          cates_tr = NULL, cates_hon = NULL,
                          method = "aipw", scores = NULL,
                          ...) {
  ## Handling inputs and checks.
  if (any(!(D_tr %in% c(0, 1)))) stop("Invalid 'D_tr'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.null(D_hon) & any(!(D_hon %in% c(0, 1)))) stop("Invalid 'D_hon'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.null(X_tr) &!is.matrix(X_tr) & !is.data.frame(X_tr)) stop("Invalid 'X_tr'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.null(X_hon) & !is.matrix(X_hon) & !is.data.frame(X_hon)) stop("Invalid 'X_hon'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.null(Y_hon) & (is.null(D_hon) | is.null(X_hon))) stop("Either you provide all 'Y_hon', 'D_hon', and 'X_hon', or set all to NULL.", call. = FALSE)
  if (!is.null(D_hon) & (is.null(Y_hon) | is.null(X_hon))) stop("Either you provide all 'Y_hon', 'D_hon', and 'X_hon', or set all to NULL.", call. = FALSE)
  if (!is.null(X_hon) & (is.null(Y_hon) | is.null(D_hon))) stop("Either you provide all 'Y_hon', 'D_hon', and 'X_hon', or set all to NULL.", call. = FALSE)
  if (!is.null(cates_tr) & is.null(cates_hon)) stop("If you provide 'cates_tr' you must also provide 'cates_hon'. Alternatively, set both to NULL and CATEs will be estimated internally.", call. = FALSE)

  ## If necessary, estimate the CATEs using training sample (estimation step).
  if (is.null(cates_tr)) {
    forest <- grf::causal_forest(X_tr, Y_tr, D_tr, num.trees = 4000, tune.parameters = "all")
    cates_tr <- stats::predict(forest, X_tr)$predictions
    cates_hon <- stats::predict(forest, X_hon)$predictions
  }

  ## Grow the tree using training sample (tree-growing step).
  tree <- rpart::rpart(cates ~ ., data = data.frame("cates" = cates_tr, X_tr), method = "anova", control = rpart::rpart.control(...), model = TRUE)

  ## If adaptive, replace each node with predictions computed in training sample. Otherwise, honest trees.
  if (is.null(Y_hon)) {
    results <- estimate_rpart(tree, Y_tr, D_tr, X_tr, method, scores = scores)
  } else {
    results <- estimate_rpart(tree, Y_hon, D_hon, X_hon, method, scores = scores)
  }

  new_tree <- results$tree
  scores <- results$scores

  ## Output.
  if (!is.null(Y_hon)) forest <- NULL

  training_sample <- data.frame(cates_tr, Y_tr, D_tr, X_tr)
  colnames(training_sample) <- c("cates", "Y", "D", colnames(X_tr))

  if (!is.null(Y_hon)) {
    honest_sample <- data.frame(cates_hon, Y_hon, D_hon, X_hon)
    colnames(honest_sample) <- colnames(training_sample)
  } else {
    honest_sample <- NULL
  }

  out <- list("tree" = new_tree,
              "forest" = forest,
              "scores" = scores,
              "method" = method,
              "training_sample" = training_sample,
              "honest_sample" = honest_sample)
  class(out) <- "aggTrees"
  return(out)
}


#' Analyzing Aggregation Trees
#'
#' @param object An \code{aggTrees} object.
#' @param n_groups Number of desired groups.
#' @param boot_ci Logical, whether to compute bootstrap confidence intervals.
#' @param boot_R Number of bootstrap replications. Ignored if \code{boot_ci == FALSE}.
#'
#' @rdname build_aggtree
#'
#' @importFrom broom tidy
#' @importFrom boot boot
#'
#' @export
inference_aggtree <- function(object, n_groups,
                              boot_ci = FALSE, boot_R = 2000) {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("Invalid 'object'. This must be an aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("Invalid 'object'. This must be a valid aggTrees object.", call. = FALSE)
  if (is.null(object$honest_sample)) warning("Inference is not valid, because the same data have been used to construct the tree and estimate GATEs.")
  if (n_groups <= 1) stop("Invalid 'n_groups'. This must be greater than or equal to 2.", call. = FALSE)
  if (!(boot_ci %in% c(FALSE, TRUE))) stop("Invalid 'boot_ci'. This must be either FALSE or TRUE.", call. = FALSE)
  if (boot_R < 0 | boot_R %% 1 != 0) stop("Invalid 'boot_R'. This must be a positive integer.", call. = FALSE)

  tree <- object$tree

  if (n_groups > get_leaves(tree)) stop("The sequence you provided does not contain any grouping with ", n_groups, " groups.", call. = FALSE)

  method <- object$method
  scores <- object$scores

  ## Select appropriate sample (full/honest) according to the output of build_aggtree.
  if (is.null(object$honest_sample)) {
    Y <- object$training_sample$Y
    D <- object$training_sample$D
    X <- object$training_sample[, -c(1:3)]
  } else {
    Y <- object$honest_sample$Y
    D <- object$honest_sample$D
    X <- object$honest_sample[, -c(1:3)]
  }

  ## Select granularity level.
  groups <- subtree(tree, leaves = n_groups)

  ## GATEs point estimates and standard errors, hypotheses testing, and boostrap CI if required.
  results <- causal_ols_rpart(groups, Y, D, X, method = method, scores = scores, boot_ci = boot_ci, boot_R = boot_R)

  model <- results$model
  gates_diff_pairs <- results$gates_diff_pairs
  boot_ci <- results$boot_ci

  ## Compute average characteristics of units in each leaf.
  avg_characteristics <- avg_characteristics_rpart(groups, X)

  ## Output.
  output <- list("aggTree" = object,
                 "groups" = groups,
                 "model" = model,
                 "boot_ci" = boot_ci,
                 "gates_diff_pairs" = gates_diff_pairs,
                 "avg_characteristics" = avg_characteristics)
  class(output) <- "aggTrees.inference"
  return(output)
}
