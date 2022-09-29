#' Honest aggTrees
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
#' \code{honest_estimation} replaces the leaf estimates of a tree as follows. First, it "pushes" all observations in
#' \code{X} down the tree and finds the leaves where they fall. Then, it computes the average value of \code{cates} in
#' each leaf. These values are the new estimates, and are honest if \code{X} has not been used to construct the tree.\cr
#'
#' Due to coding limitations, honest trees from \code{aggTrees} objects show honest estimates only in their leaves. The
#' internal nodes show non-honest estimates, so the user should ignore them and focus on the leaves of tree.\cr
#'
#' To get standard errors on the tree's estimates, please use \code{\link{honest_ols}}.
#'
#' @import treeClust Rcpp
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{honest_ols}},\code{\link{aggregation_tree}}
#'
#' @export
honest_estimation <- function(object, cates, X) {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- object$tree

  ## Honest re-estimation.
  honest_tree <- honest_estimation_rpart(tree, cates, X)

  ## Output.
  out <- list("tree" = honest_tree, "honesty" = TRUE)
  class(out) <- "aggTrees"
  return(out)
}


#' Honest rpart
#'
#' Replace leaf predictions of a \code{\link[rpart]{rpart}} object using external data.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param y_honest Outcome vector of the honest sample.
#' @param X_honest Covariate matrix (no intercept) of the honest sample.
#'
#' @return
#' An honest tree, as a \code{\link[rpart]{rpart}} object.
#'
#' @details
#' Due to coding limitations, the honest tree shows honest estimates only in its leaves. The internal nodes show
#' non-honest estimates, so the user should ignore them and focus on the leaves of tree.\cr
#'
#' To get standard errors on the tree's estimates, please use \code{\link{honest_ols_rpart}}.
#'
#' @import treeClust Rcpp
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{honest_ols_rpart}}
#'
#' @export
honest_estimation_rpart <- function(tree, y_honest, X_honest) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  honest_tree <- tree

  ## Extract leaves where honest observations fall.
  leaf_honest <- treeClust::rpart.predict.leaves(tree, data.frame(X_honest), type = "where") # Row numbers of tree$frame.
  names(leaf_honest) <- 1:length(y_honest)
  unique_leaves_honest <- unique(leaf_honest)

  ## Call cpp to compute honest leaf estimates.
  honest_estimates <- as.matrix(honest_rpart_cpp(unique_leaves_honest, y_honest, leaf_honest))

  ## Replace leaf estimates.
  for (leaf in unique(leaf_honest)) {
    honest_tree$frame$yval[leaf] <- honest_estimates[as.numeric(names(leaf_honest)[leaf_honest == leaf][1])]
    honest_tree$frame$n[leaf] <- sum(leaf_honest == leaf)
  }

  ## Output.
  return(honest_tree)
}


#' Honest OLS Estimation aggTrees
#'
#' Uses the leaves of a tree stored in an \code{aggTrees} object to estimate a linear model via OLS. If the data used in
#' the OLS estimation have not been used to grow the tree, then standard errors for the tree's estimates are valid.
#'
#' @param object An \code{aggTrees} object.
#' @param cates Estimated CATEs.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' The fitted model, as a \code{\link[stats]{lm}} object.
#'
#' @details
#' To get standard errors on the tree's estimates, \code{honest_ols} fits via OLS the following model:
#'
#' \deqn{cates_i = \sum_{l = 1}^{|T|} L_l \beta_l + \epsilon_i}
#'
#' with \code{L_l} the l-th leaf of the tree, and \code{|T|} the number of leaves. It is immediate to notice that the l-th
#' coefficient corresponds to the GATEs of the l-th group. Thus, standard errors on the estimated coefficients are standard
#' errors on the estimated GATEs.\cr
#'
#' Notice that honesty is a necessary requirement to get valid standard errors. Thus, \code{X} must not have been used to
#' grow the tree or estimate the cates.\cr
#'
#' Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{honest_estimation}}
#'
#' @export
honest_ols <- function(object, cates, X) {
  ## Handling inputs and checks
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- object$tree

  ## Fit the model.
  model <- honest_ols_rpart(tree, cates, X)

  ## Output
  return(model)
}


#' Honest OLS Estimation rpart
#'
#' Uses the leaves of a tree stored in an \code{\link[rpart]{rpart}} object to estimate a linear model via OLS.
#' If the data used in the OLS estimation have not been used to grow the tree, then standard errors for the tree's
#' estimates are valid.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param y_honest Outcome vector of the honest sample.
#' @param X_honest Covariate matrix (no intercept) of the honest sample.
#'
#' @return
#' The fitted model, as a \code{\link[stats]{lm}} object.
#'
#' @details
#' To get standard errors on the tree's estimates, \code{honest_ols_rpart} fits via OLS the following model:
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
#' @seealso \code{\link{aggregation_tree}}, \code{\link{honest_estimation}}
#'
#' @export
honest_ols_rpart <- function(tree, y_honest, X_honest) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  ## Generate leaves indicators. Inspired by https://bookdown.org/halflearned/tutorial/hte1.html.
  tree_predictions <- predict(tree, data.frame(X_honest))
  n_leaves <- length(unique(tree_predictions))
  leaves <- factor(tree_predictions, levels = sort(unique(tree_predictions)), labels = seq(n_leaves))

  ## OLS estimation.
  model <- estimatr::lm_robust(y ~ 0 + leaf, data = data.frame("y" = y_honest, "leaf" = leaves), se_type = "HC1")

  ## Output.
  return(model)
}
