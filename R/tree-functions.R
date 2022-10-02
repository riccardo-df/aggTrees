#' Aggregation Trees
#'
#' Grows an aggregation tree for discovering heterogeneous sub-populations that differ in the magnitude
#' of their treatment effects.
#'
#' @param cates Estimated CATEs.
#' @param X Covariate matrix (no intercept).
#' @param maxdepth The maximum depth of the tree.
#' @param ... Further arguments from \code{\link[rpart]{rpart.control}}.
#'
#' @return
#' An \code{aggTrees} object.
#'
#' @details
#' Aggregation trees are a three-step procedure. First, CATEs are estimated using any estimator. Second, a tree is grown
#' to approximate the CATEs. Third, the tree is pruned to derive a nested sequence of optimal partitions. Different
#' subsamples are used: an "estimation sample" for the CATEs estimation, and an "aggregation sample" for the tree-growing
#' and the tree-pruning.\cr
#'
#' \code{aggregation_tree} uses \code{X} to grow a tree approximating the CATEs provided by the user. For this, the user
#' passes in as an argument a vector of CATEs \code{cates} obtained by fitting a model in the estimation sample and
#' predicting in the aggregation sample \code{X}. \code{X} must consists of a covariate matrix with variables that are a
#' subset of those used to construct \code{cates}. Importantly, observations in \code{X} must not have been used to
#' construct \code{cates}.\cr
#'
#' The tree is grown up to some stopping criteria that can be specified by the user. Please refer to
#' the \code{\link[rpart]{rpart.control}} documentation for this.\cr
#'
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso
#' \code{\link{honest_estimation}}, \code{\link{honest_ols}}, \code{\link{plot.aggTrees}}, \code{\link{subtree}}
#'
#' @export
aggregation_tree <- function(cates, X, maxdepth = 3, ...) {
  ## Handling inputs and checks.
  if (length(cates) == 0) stop("'cates' has length zero.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)
  if (maxdepth < 1 | maxdepth > 30) stop("'maxdepth' must be in the interval [1, 30].", call. = FALSE)

  ## Grow the tree using only aggregation sample.
  tree <- rpart::rpart(cates ~ ., data = data.frame("cates" = cates, X), method = "anova", control = rpart::rpart.control(maxdepth = maxdepth, ...), model = TRUE)

  ## Output.
  out <- list("tree" = tree, "honesty" = FALSE)
  class(out) <- "aggTrees"
  return(out)
}


#' aggTrees Subtree
#'
#' Extracts a subtree with a user-specified number of leaves from an \code{aggTrees} object.
#'
#' @param x An \code{aggTrees} object.
#' @param leaves Number of leaves of the desired subtree. Default is number of leaves of full tree.
#' @param cv If \code{TRUE}, \code{leaves} is ignored and a cross-validation criterion is used to select a partition.
#'
#' @return
#' An \code{aggTrees} object.
#'
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree_rpart}}, \code{\link{plot.aggTrees}}
#'
#' @export
subtree <- function(x, leaves = NULL, cv = FALSE) {
  ## Handling inputs and checks.
  if (!(inherits(x, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if(!(inherits(x$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  if (!(cv %in% c(TRUE, FALSE))) stop("'cv' must be either TRUE or FALSE.", call. = FALSE)
  if (is.null(leaves) & cv == FALSE) stop("Invalid combination of 'leaves' and 'cv'. Please specify a number of leaves or select the cross-validation option.", call. = FALSE)

  tree <- x$tree

  if (!is.null(leaves)) {
    if (leaves < 1) stop("'leaves' must be a positive number.", call. = FALSE)
    if (leaves > get_leaves(tree)) stop("'leaves' is greater than the number of leaves of 'tree'. Please provide a deeper 'tree'.", call. = FALSE)
  }

  ## Selecting subtree.
  subtree <- subtree_rpart(tree, leaves, cv)

  ## Output.
  out <- list("tree" = subtree, "honesty" = x$honesty)
  class(out) <- "aggTrees"
  return(out)
}


#' rpart Subtree
#'
#' Extracts a subtree with a user-specified number of leaves from an \code{\link[rpart]{rpart}} object.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param leaves Number of leaves of the desired subtree. Default is number of leaves of full tree.
#' @param cv If \code{TRUE}, \code{leaves} is ignored and a cross-validation criterion is used to select a partition.
#'
#' @return
#' The subtree, as a \code{\link[rpart]{rpart}} object.
#'
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree}}, \code{\link{plot.aggTrees}}
#'
#' @export
subtree_rpart <- function(tree, leaves = NULL, cv = FALSE) {
  ## Handling inputs and checks.
  if(!(inherits(tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(cv %in% c(TRUE, FALSE))) stop("'cv' must be either TRUE or FALSE.", call. = FALSE)
  if (is.null(leaves) & cv == FALSE) stop("Invalid combination of 'leaves' and 'cv'. Please specify a number of leaves or select the cross-validation option.", call. = FALSE)
  if (!is.null(leaves)) {
    if (leaves < 1) stop("'leaves' must be a positive number.", call. = FALSE)
    if (leaves > get_leaves(tree)) stop("'leaves' is greater than the number of leaves of 'tree'. Please provide a deeper 'tree'.", call. = FALSE)
  }

  ## Output.
  if (cv) return(rpart::prune(tree, tree$cptable[, 1][which.min(tree$cptable[, 4])])) else return(rpart::prune(tree, tree$cptable[tree$cptable[, "nsplit"] == leaves - 1, "CP"]))
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
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @export
get_leaves <- function(tree) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.")

  return(dim(tree$frame[tree$frame$var == "<leaf>", ])[1])
}
