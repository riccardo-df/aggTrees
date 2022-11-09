#' Estimating Leaf-Effects via Linear Models
#'
#' Uses the leaves of a tree stored in an \code{aggTrees} object to estimate a linear model via OLS. The estimated
#' coefficients correspond to the GATEs in each leaf. If the data used in the OLS estimation have not been used to
#' grow the tree (a condition called "honesty"), then one can use the standard errors for the tree's estimates to
#' conduct valid inference using conventional approaches.
#'
#' @param object An \code{aggTrees} object.
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' The fitted model, as a \code{\link[estimatr]{lm_robust}} object.
#'
#' @details
#' In randomized experiments, taking the difference between the sample average of the observed outcomes of treated units
#' and the sample average of the observed outcomes of control units is an unbiased estimators of Average Treatment Effect
#' (ATE). One can use this estimator in different subpopulations to estimate the Group Average Treatment Effects (GATEs).
#' Notice that this is equivalent to regress the outcome on a set of dummies denoting group membership and the interactions
#' of these dummies on the binary treatment indicator. One of the advantage of using the regression is that one gets standard
#' errors for the GATEs.\cr
#'
#' \code{causal_ols_aggtree} fits via OLS the following model:
#'
#' \deqn{y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th leaf of the tree, and \code{|T|} the
#' number of leaves. It is immediate to notice that, in randomized experiments, the estimated betas correspond to the GATEs
#' of the l-th group:
#'
#' \deqn{E[Y | D = 1, L_l = 1] - E[Y | D = 0, L_l = 1] = \gamma_l + \beta_l - \gamma_l = \beta_l}
#'
#' Thus, standard errors for the estimated betas are standard errors for the estimated GATEs. Standard errors are estimated
#' via the Eicker-Huber-White estimator.\cr
#'
#' Notice that "honesty" is a necessary requirement to get valid inference. Thus, observations in \code{y}, \code{X}, and
#' \code{D} must not have been used to grow the tree.\cr
#'
#' If the tree consists of a root only, \code{causal_ols_aggtree} regresses \code{y} on a constant and \code{D}, thus
#' estimating the ATE.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item S Athey, G Imbens (2016). Recursive partitioning for heterogeneous causal effects. Proceedings of the National Academy of Sciences. \doi{10.1073/pnas.1510489113}.
#' }
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{estimate_aggtree}}
#'
#' @export
causal_ols_aggtree <- function(object, y, X, D) {
  ## Handling inputs and checks
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- object$tree

  ## Fit the model.
  model <- causal_ols_rpart(tree, y, X, D)

  ## Output
  return(model)
}


#' Estimating Leaf-Effects via Linear Models
#'
#' Uses the leaves of a tree stored in an \code{\link[rpart]{rpart}} object to estimate a linear model via OLS. The
#' estimated coefficients correspond to the GATEs in each leaf. If the data used in the OLS estimation have not been
#' used to grow the tree (a condition called "honesty"), then one can use the standard errors for the tree's estimates
#' to conduct valid inference using conventional approaches.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' The fitted model, as a \code{\link[estimatr]{lm_robust}} object.
#'
#' @details
#' In randomized experiments, taking the difference between the sample average of the observed outcomes of treated units
#' and the sample average of the observed outcomes of control units is an unbiased estimators of Average Treatment Effect
#' (ATE). One can use this estimator in different subpopulations to estimate the Group Average Treatment Effects (GATEs).
#' Notice that this is equivalent to regress the outcome on a set of dummies denoting group membership and the interactions
#' of these dummies on the binary treatment indicator. One of the advantage of using the regression is that one gets standard
#' errors for the GATEs.\cr
#'
#' \code{causal_ols_rpart} fits via OLS the following model:
#'
#' \deqn{y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th leaf of the tree, and \code{|T|} the
#' number of leaves. It is immediate to notice that, in randomized experiments, the estimated betas correspond to the GATEs
#' of the l-th group:
#'
#' \deqn{E[Y | D = 1, L_l = 1] - E[Y | D = 0, L_l = 1] = \gamma_l + \beta_l - \gamma_l = \beta_l}
#'
#' Thus, standard errors for the estimated betas are standard errors for the estimated GATEs. Standard errors are estimated
#' via the Eicker-Huber-White estimator.\cr
#'
#' Notice that "honesty" is a necessary requirement to get valid inference. Thus, observations in \code{y}, \code{X}, and
#' \code{D} must not have been used to grow the tree.\cr
#'
#' If the tree consists of a root only, \code{causal_ols_aggtree} regresses \code{y} on a constant and \code{D}, thus
#' estimating the ATE.
#'
#' @import rpart estimatr grf
#' @importFrom stats predict
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item S Athey, G Imbens (2016). Recursive partitioning for heterogeneous causal effects. Proceedings of the National Academy of Sciences. \doi{10.1073/pnas.1510489113}.
#' }
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{estimate_rpart}}
#'
#' @export
causal_ols_rpart <- function(tree, y, X, D) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  ## Generate leaves indicators.
  leaves <- leaf_membership(tree, X)

  if (length(unique(leaves)) < get_leaves(tree)) warning("One or more leaves are empty: No observations in X fall there.")

  ## OLS estimation.
  if (length(unique(leaves)) == 1) {
    model <- estimatr::lm_robust(y ~ D, data = data.frame("y" = y, "D" = D), se_type = "HC1")
  } else {
    model <- estimatr::lm_robust(y ~ 0 + leaf + D:leaf, data = data.frame("y" = y, "leaf" = leaves, "D" = D), se_type = "HC1")
  }

  ## Output.
  return(model)
}


#' Estimation of aggTrees Objects
#'
#' Replace leaf predictions of a tree stored in an \code{aggTrees} object using external data.
#'
#' @param object An \code{aggTrees} object.
#' @param X Covariate matrix (no intercept).
#' @param y Outcome vector.
#' @param D Treatment assignment vector.
#' @param cates Estimated CATEs.
#' @param method Either \code{"raw"} or \code{"cates"}, defines how leaf predictions are replaced.
#'
#' @return
#' An \code{aggTrees} object.
#'
#' @details
#' \code{estimate_aggtree} replaces the leaf estimates of a tree as follows. First, it "pushes" all observations in
#' \code{X} down the tree and finds the leaves where they fall. Then, it replaces the predictions in each leaf according
#' to the user-specified \code{method}.
#'
#' If \code{method = "raw"}, \code{estimate_aggtree} replaces leaf predictions with the difference between the sample average
#' of the observed outcomes of treated units and the sample average of the observed outcomes of control units in each leaf,
#' which is an unbiased estimator of the GATEs if the assignment to treatment is randomized. This requires the user to
#' specify \code{y} and \code{D}.\cr
#'
#' If \code{method = "cates"}, \code{estimate_aggtree} replaces leaf predictions with the sample average of \code{cates}
#' in each leaf. This is a valid estimator of the GATEs in observational studies.\cr
#'
#' \code{estimate_aggtree} allows the user to implement "honest" estimation. If observations in \code{X} have not been
#' used to construct the tree, then the new predictions are honest in the sense of Athey and Imbens (2016). This allows
#' the user to conduct valid inference about the estimated GATEs with standard approaches, e.g., by constructing conventional
#' confidence intervals. To get standard errors for the tree's estimates, please use \code{\link{causal_ols_aggtree}}.\cr
#'
#' Due to coding limitations, \code{estimate_aggtree} replaces the estimates only in the leaves of a tree. The
#' internal nodes' estimates are not replaced, so the user should ignore them and focus on the leaves of the tree.\cr
#'
#' @import treeClust Rcpp
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item S Athey, G Imbens (2016). Recursive partitioning for heterogeneous causal effects. Proceedings of the National Academy of Sciences. \doi{10.1073/pnas.1510489113}.
#' }
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{causal_ols_aggtree}}
#'
#' @export
estimate_aggtree <- function(object, X, y = NULL, D = NULL, cates = NULL, method = "raw") {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  if (!(method %in% c("raw", "cates"))) stop("You must provide a valid method.", call. = FALSE)
  if (method == "raw" & (is.null(y) | is.null(D))) stop("'raw' method requires the user to specify both 'y' and 'D'.", call. = FALSE)
  if (method == "cates" & (is.null(cates))) stop("'cates' method requires the user to specify 'cates'.", call. = FALSE)

  tree <- object$tree

  ## Honest re-estimation.
  new_tree <- estimate_rpart(tree, X, y, D, cates, method)

  ## Output.
  out <- list("tree" = new_tree, "honesty" = TRUE)
  class(out) <- "aggTrees"
  return(out)
}


#' Estimation of rpart Objects
#'
#' Replace leaf predictions of a \code{\link[rpart]{rpart}} object using external data.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param X Covariate matrix (no intercept).
#' @param y Outcome vector.
#' @param D Treatment assignment vector.
#' @param cates Estimated CATEs.
#' @param method Either \code{"raw"} or \code{"cates"}, defines how leaf predictions are replaced.
#'
#' @return
#' A tree with leaf predictions replaced, as a \code{\link[rpart]{rpart}} object.
#'
#' @details
#' \code{estimate_rpart} replaces the leaf estimates of a tree as follows. First, it "pushes" all observations in
#' \code{X} down the tree and finds the leaves where they fall. Then, it replaces the predictions in each leaf according
#' to the user-specified \code{method}.
#'
#' If \code{method = "raw"}, \code{estimate_rpart} replaces leaf predictions with the difference between the sample average
#' of the observed outcomes of treated units and the sample average of the observed outcomes of control units in each leaf,
#' which is an unbiased estimator of the GATEs if the assignment to treatment is randomized. This requires the user to
#' specify \code{y} and \code{D}.\cr
#'
#' If \code{method = "cates"}, \code{estimate_rpart} replaces leaf predictions with the sample average of \code{cates}
#' in each leaf. This is a valid estimator of the GATEs in observational studies.\cr
#'
#' \code{estimate_rpart} allows the user to implement "honest" estimation. If observations in \code{X} have not been used
#' to construct the tree, then the new predictions are honest in the sense of Athey and Imbens (2016). This allows the user
#' to conduct valid inference about the estimated GATEs with standard approaches, e.g., by constructing conventional
#' confidence intervals. To get standard errors for the tree's estimates, please use \code{\link{causal_ols_rpart}}.\cr
#'
#' Due to coding limitations, \code{estimate_rpart} replaces the estimates only in the leaves of a tree. The
#' internal nodes' estimates are not replaced, so the user should ignore them and focus on the leaves of the tree.\cr
#'
#' @import treeClust Rcpp causalTree
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item S Athey, G Imbens (2016). Recursive partitioning for heterogeneous causal effects. Proceedings of the National Academy of Sciences. \doi{10.1073/pnas.1510489113}.
#' }
#'
#' @seealso \code{\link{causal_ols_rpart}}
#'
#' @export
estimate_rpart <- function(tree, X, y = NULL, D = NULL, cates = NULL, method = "raw") {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  if (!(method %in% c("raw", "cates"))) stop("You must provide a valid method.", call. = FALSE)
  if (method == "raw" & (is.null(y) | is.null(D))) stop("'raw' method requires the user to specify both 'y' and 'D'.", call. = FALSE)
  if (method == "cates" & (is.null(cates))) stop("'cates' method requires the user to specify 'cates'.", call. = FALSE)

  new_tree <- tree

  if (method == "raw") {
    new_tree <- causalTree::estimate.causalTree(tree, data.frame(X, y), treatment = D)
  } else if (method == "cates") {
    ## Extract leaves where observations in X fall.
    leaves <- treeClust::rpart.predict.leaves(tree, data.frame(X), type = "where") # Row numbers of tree$frame.
    names(leaves) <- 1:length(y)
    unique_leaves <- unique(leaves)

    ## Call cpp to compute honest leaf estimates.
    honest_estimates <- as.matrix(honest_rpart_cpp(unique_leaves, y, leaves))

    ## Replace leaf estimates.
    for (leaf in unique(leaves)) {
      new_tree$frame$yval[leaf] <- honest_estimates[as.numeric(names(leaves)[leaves == leaf][1])]
      new_tree$frame$n[leaf] <- sum(leaves == leaf)
    }
  }

  ## Output.
  return(new_tree)
}
