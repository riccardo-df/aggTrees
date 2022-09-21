#' Aggregation Trees
#'
#' Grows an aggregation tree.
#'
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept).
#' @param maxdepth The maximum depth of the tree (i.e, the maximum number of nodes connecting the root to the leaves).
#' @param cp Any split that does not decrease the overall lack of fit by a factor of cp is not attempted. Default is zero, meaning that trees are fully grown up to \code{maxdepth}.
#'
#' @return
#' The tree, as a \code{\link[rpart]{rpart}} object.
#'
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{subtree}}, \code{\link{plot_tree}}
#'
#' @export
aggregation_tree <- function(cates, X, maxdepth, cp = 0) {
  ## Checks.
  if (!is.numeric(cates)) stop("'cates' must be a numeric vector.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)
  if (maxdepth < 1 | maxdepth > 30) stop("'maxdepth' must be in the interval [1, 30].", call. = FALSE)

  ## Growing tree.
  tree <- rpart::rpart(cates ~ ., data = data.frame("cates" = cates, X), method = "anova", control = rpart::rpart.control(maxdepth = maxdepth, cp = cp), model = TRUE)

  ## Output.
  return(tree)
}


#' Subtree
#'
#' Extracts a subtree with a user-specified number of leaves.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param leaves Number of leaves of the desired subtree. Default is number of leaves of full tree.
#' @param cv Logical. If \code{TRUE}, \code{leaves} is ignored and a cross-validation criterion is used to select a partition.
#'
#' @return
#' The subtree, as a \code{\link[rpart]{rpart}} object.
#'
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree}}, \code{\link{plot_tree}}
#'
#' @export
subtree <- function(tree, leaves = NULL, cv = FALSE) {
  ## Checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)
  if (is.null(leaves) & cv == FALSE) stop("Invalid combination of 'leaves' and 'cv'. Please specify a number of leaves or select the cross-validation option.", call. = FALSE)
  if (!is.null(leaves)) if (leaves < 1) stop("'leaves' must be a positive number.", call. = FALSE)
  if (!is.null(leaves) & leaves > get_leaves(tree)) stop("'leaves' is greater than the number of leaves of 'tree'. Please provide a deeper 'tree'.", call. = FALSE)

  ## Output.
  if (cv) return(rpart::prune(tree, tree$cptable[, 1][which.min(tree$cptable[, 4])])) else return(rpart::prune(tree, tree$cptable[tree$cptable[, "nsplit"] == leaves - 1, "CP"])) # Number of leaves in cptable is nsplit + 1.
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
#' If \code{tree} is just a root, then \code{get_leaves} returns 1.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @export
get_leaves <- function(tree) {
  if (!inherits(tree, "rpart")) stop("'tree must be a rpart object.")

  return(dim(tree$frame[tree$frame$var == "<leaf>", ])[1])
}
