#' Aggregation Trees
#'
#' Grows an aggregation tree.
#'
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept).
#' @param maxdepth The maximum depth (i.e, the maximum number of nodes connecting the root to the leaves) of the tree.
#' @param cp Minimum MSE increase to accept splits.
#' @param honest Logical. Whether honest estimation should be implemented.
#' @param honest_frac Number in the (0, 1] interval. The fraction of observations to be used in estimating the tree structure. Ignored if \code{honesty} is \code{FALSE}.
#'
#' @return
#' A \code{\link[rpart]{rpart}} object.
#'
#' @details
#' For details about aggregation trees and their use, please see Di Francesco (2022). \cr
#'
#' @examples
#' ## Loading data (using random subsample provided by Matias Cattaneo).
#' dta <- haven::read_dta("http://www.stata-press.com/data/r13/cattaneo2.dta")
#'
#' Y <- as.matrix(dta[, "bweight"])
#' D <- as.matrix(dta[, "mbsmoke"])
#' X_names <- c("bweight", "mbsmoke", "deadkids", "monthslb", "lbweight")
#' X <- as.matrix(dta[, !(colnames(dta) %in% X_names)])
#'
#' ## Splitting sample.
#' set.seed(1986)
#' n <- dim(dta)[1]
#' est_idx <- sample(1:n, n / 2, replace = FALSE)
#'
#' X_est <- X[est_idx, ]
#' Y_est <- Y[est_idx]
#' D_est <- D[est_idx]
#'
#' X_agg <- X[-est_idx, ]
#' Y_agg <- Y[-est_idx]
#' D_agg <- D[-est_idx]
#'
#' ## Estimating CATEs using only estimation sample.
#' cates_forest <- grf::causal_forest(X = X_est, Y = Y_est, W = D_est)
#'
#' ## Growing tree using only aggregation sample.
#' cates <- predict(cates_forest, newdata = X_agg)$predictions
#' tree <- aggregation_tree(cates, X_agg, maxdepth = 3, cp = 0.01)
#'
#' ## Plotting.
#' plot_aggregation_tree(tree)
#'
#' ## It can be annoying selecting a palette each time.
#' palette <- colorspace::choose_palette()
#' plot_aggregation_tree(tree, palette)
#' plot_aggregation_tree(tree, palette)
#'
#' ## Displaying the whole sequence.
#' plot_aggregation_tree(tree, palette, sequence = TRUE)
#'
#' @export
aggregation_tree <- function(cates, X, maxdepth, cp, honesty = FALSE, honest_frac = 0.5) {
  ## Growing tree.
  if (honesty) {
    honest_idx <- sample(1:dim(X)[1], floor(dim(X)[1] * honest_frac), replace = FALSE)

    tree <- rpart::rpart(cates ~ .,
                         data = data.frame("cates" = cates[honest_idx], X[honest_idx, ]), method = "anova",
                         control = rpart::rpart.control(maxdepth = maxdepth, cp = cp), model = TRUE)
  } else {
    tree <- rpart::rpart(cates ~ .,
                         data = data.frame("cates" = cates, X), method = "anova",
                         control = rpart::rpart.control(maxdepth = maxdepth, cp = cp), model = TRUE)
  }

  ## Output.
  return(tree)
}


#' Cross-Validated Tree
#'
#' Uses the cross-validation criterion to select the "best" subtree.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#'
#' @return
#' The cross-validated tree, as a \code{\link[rpart]{rpart}} object.
#'
#' @details
#' \code{tree} should be deep enough to allow the criterion to explore more trees. \cr
#'
#' @examples
#' ## Loading data (using random subsample provided by Matias Cattaneo).
#' dta <- haven::read_dta("http://www.stata-press.com/data/r13/cattaneo2.dta")
#'
#' Y <- as.matrix(dta[, "bweight"])
#' D <- as.matrix(dta[, "mbsmoke"])
#' X_names <- c("bweight", "mbsmoke", "deadkids", "monthslb", "lbweight")
#' X <- as.matrix(dta[, !(colnames(dta) %in% X_names)])
#'
#' ## Splitting sample.
#' set.seed(1986)
#' n <- dim(dta)[1]
#' est_idx <- sample(1:n, n / 2, replace = FALSE)
#'
#' X_est <- X[est_idx, ]
#' Y_est <- Y[est_idx]
#' D_est <- D[est_idx]
#'
#' X_agg <- X[-est_idx, ]
#' Y_agg <- Y[-est_idx]
#' D_agg <- D[-est_idx]
#'
#' ## Estimating CATEs using only estimation sample.
#' cates_forest <- grf::causal_forest(X = X_est, Y = Y_est, W = D_est)
#'
#' ## Growing tree using only aggregation sample.
#' cates <- predict(cates_forest, newdata = X_agg)$predictions
#' tree <- aggregation_tree(cates, X_agg, maxdepth = 3, cp = 0.01)
#'
#' ## Plotting.
#' plot_aggregation_tree(tree)
#'
#' ## Cross-validation.
#' cv_tree <- cross_validated_tree(tree)
#' plot_aggregation_tree(cv_tree)
#'
#' @export
cross_validated_tree <- function(tree) {
  if (!inherits(tree, "rpart")) stop("The input must be a rpart object.")

  return(rpart::prune(tree, tree$cptable[, 1][which.min(tree$cptable[, 4])]))
}


#' Number of Leaves in Tree
#'
#' Extracts the number of leaves of a given tree.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#'
#' @return
#' The number of leaves.
#'
#' @details
#' If the tree is just a root, then \code{get_leaves} returns 1.
#'
#' @export
get_leaves <- function(tree) {
  if (!inherits(tree, "rpart")) stop("The input must be a rpart object.")

  return(dim(tree$frame[tree$frame$var == "<leaf>", ])[1])
}

