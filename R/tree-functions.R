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
#' @seealso \code{\link{aggregation_tree}}
#'
#' @export
subtree <- function(tree, leaves = NULL, cv = FALSE) {
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
#' Extracts the number of leaves of an \code{\link[rpart]{rpart}} object.
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


#' Leaf Membership
#'
#' Constructs a variable that encodes in which leaf of an \code{\link[rpart]{rpart}} object the units in a given data frame fall.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' A factor whose levels denote in which leaf each unit falls.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @export
leaf_membership <- function(tree, X) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.")
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)

  ## Inspired by https://bookdown.org/halflearned/tutorial/hte1.html.
  tree_predictions <- predict(tree, data.frame(X))
  n_leaves <- length(unique(tree_predictions))
  leaves <- factor(tree_predictions, levels = sort(unique(tree_predictions)), labels = seq(n_leaves))

  return(leaves)
}
