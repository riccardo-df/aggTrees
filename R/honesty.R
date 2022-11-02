#' Estimation of aggTrees Objects
#'
#' Replace leaf predictions of a tree stored in an \code{aggTrees} object using external data.
#'
#' @param object An \code{aggTrees} object.
#' @param cates Estimated CATEs.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' An \code{aggTrees} object.
#'
#' @details
#' \code{estimate_aggtree} replaces the leaf estimates of a tree as follows. First, it "pushes" all observations in
#' \code{X} down the tree and finds the leaves where they fall. Then, it computes the average value of \code{cates} in
#' each leaf. These values are the new estimates.\cr
#'
#' If observations in \code{X} have not been used to construct the tree, then the new predictions are "honest".\cr
#'
#' Due to coding limitations, \code{estimate_aggtree} replaces the estimates only in the leaves of a tree. The
#' internal nodes' estimates are not replaced, so the user should ignore them and focus on the leaves of the tree.\cr
#'
#' To get standard errors on the tree's estimates, please use \code{\link{ols_aggtree}}.
#'
#' @import treeClust Rcpp
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{ols_aggtree}}, \code{\link{aggregation_tree}}
#'
#' @export
estimate_aggtree <- function(object, cates, X) {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- object$tree

  ## Honest re-estimation.
  new_tree <- estimate_rpart(tree, cates, X)

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
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' A tree with leaf predictions replaced, as a \code{\link[rpart]{rpart}} object.
#'
#' @details
#' Due to coding limitations, \code{estimate_rpart} replaces the estimates only in the leaves of a tree. The
#' internal nodes' estimates are not replaced, so the user should ignore them and focus on the leaves of the tree.\cr
#'
#' To get standard errors on the tree's estimates, please use \code{\link{ols_rpart}}.
#'
#' @import treeClust Rcpp
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{ols_rpart}}
#'
#' @export
estimate_rpart <- function(tree, y, X) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  new_tree <- tree

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

  ## Output.
  return(new_tree)
}


#' OLS Estimation of aggTrees Objects
#'
#' Uses the leaves of a tree stored in an \code{aggTrees} object to estimate a linear model via OLS. If the data used in
#' the OLS estimation have not been used to grow the tree, then standard errors for the tree's estimates are valid.
#'
#' @param object An \code{aggTrees} object.
#' @param cates Estimated CATEs.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' The fitted model, as a \code{\link[estimatr]{lm_robust}} object.
#'
#' @details
#' To get standard errors on the tree's estimates, \code{ols_aggtree} fits via OLS the following model:
#'
#' \deqn{cates_i = \sum_{l = 1}^{|T|} L_l \beta_l + \epsilon_i}
#'
#' with \code{L_l} the l-th leaf of the tree, and \code{|T|} the number of leaves. It is immediate to notice that the l-th
#' coefficient corresponds to the GATEs of the l-th group. Thus, standard errors on the estimated coefficients are standard
#' errors on the estimated GATEs.\cr
#'
#' Notice that honesty is a necessary requirement to get valid standard errors. Thus, observations in \code{X} must not have
#' been used to grow the tree or estimate the cates.\cr
#'
#' Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{estimate_aggtree}}
#'
#' @export
ols_aggtree <- function(object, cates, X) {
  ## Handling inputs and checks
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- object$tree

  ## Fit the model.
  model <- ols_rpart(tree, cates, X)

  ## Output
  return(model)
}


#' OLS Estimation of rpart Objects
#'
#' Uses the leaves of a tree stored in an \code{\link[rpart]{rpart}} object to estimate a linear model via OLS.
#' If the data used in the OLS estimation have not been used to grow the tree, then standard errors for the tree's
#' estimates are valid.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' The fitted model, as a \code{\link[estimatr]{lm_robust}} object.
#'
#' @details
#' To get standard errors on the tree's estimates, \code{ols_rpart} fits via OLS the following model:
#'
#' \deqn{y_i = \sum_{l = 1}^{|T|} L_l \beta_l + \epsilon_i}
#'
#' with \code{L_l} the l-th leaf of the tree, and \code{|T|} the number of leaves. It is immediate to notice that the l-th
#' coefficient corresponds to the estimated conditional expectation of \code{y_i} in the l-th leaf. Thus, standard errors
#' on the estimated coefficients are standard errors on these estimates.\cr
#'
#' Notice that honesty is a necessary requirement to get valid standard errors. Thus, \code{y_honest} and {X_honest} must
#' not have been used to grow the tree.\cr
#'
#' Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' @import rpart estimatr
#' @importFrom stats predict
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{estimate_rpart}}
#'
#' @export
ols_rpart <- function(tree, y, X) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  ## Generate leaves indicators.
  leaves <- leaf_membership(tree, X)

  ## OLS estimation.
  model <- estimatr::lm_robust(y ~ 0 + leaf, data = data.frame("y" = y, "leaf" = leaves), se_type = "HC1")

  ## Output.
  return(model)
}
