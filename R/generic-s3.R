#' Plot Method for aggTrees Objects
#'
#' Plots an \code{aggTrees} object.
#'
#' @param x An \code{aggTrees} object.
#' @param leaves Number of leaves of the desired tree. This can be used to plot subtrees.
#' @param type Plotting style.
#' @param palette Palette to color the nodes.
#' @param sequence If \code{TRUE}, the whole sequence of optimal partitions is displayed in a short animation.
#' @param ... Further arguments from \code{\link[rpart.plot]{prp}}.
#'
#' @details
#' Due to coding limitations, honest trees from \code{aggTrees} objects show honest estimates only in their leaves. Thus, if
#' the user wants to plot an honest tree, \code{type} 3 or 5 is recommended to avoid misleading plots.\cr
#'
#' \code{palette} can be either a vector of colors, or a function that takes as an argument the number of nodes of a tree and
#' returns a vector of colors. Please refer to \code{\link[colorspace]{choose_palette}} for high-quality palettes.\cr
#'
#' @import rpart.plot
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso
#' \code{\link{aggregation_tree}}, \code{\link{subtree_aggtree}}, \code{\link{recursive_partitioning_plot}}
#'
#' @export
plot.aggTrees <- function(x, leaves = get_leaves(x$tree),
                          type = 2, palette = c("BuRd"),
                          sequence = FALSE, ...) {
  ## Handling inputs and checks.
  if (!(inherits(x, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(x$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- x$tree

  # if (!is.null(x$idx$honest_idx) & !(type %in% c(3, 5))) warning("Only the leaf estimates are honest. Internal nodes show non-honest estimates. Consider using 'type' 3 or 5.")

  ## Plotting.
  if (sequence) {
    prefix <- c("ATE = ", rep("GATE = ", times = (length(tree$frame$n) - 1)))
    sizes <- tree$frame$n
    percentages <- round(tree$frame$n / tree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    full_nodes <- rownames(tree$frame) # Storing rownames of full tree (they serve as node labels).
    alpha_values <- rev(tree$cptable[, "CP"]) # Threshold values of cost-complexity parameter.

    for (alpha in alpha_values) {
      temp_nodes <- rownames(rpart::prune.rpart(tree, alpha)$frame) # Node labels of pruned tree.
      colors <- ifelse(full_nodes %in% temp_nodes, 1, "white") # Graying out collapsed nodes.

      grDevices::dev.hold()

      rpart.plot::prp(tree,
                      type = type,
                      extra = 0,
                      under = FALSE,
                      fallen.leaves = TRUE,
                      round = 0,
                      leaf.round = 0,
                      prefix = prefix,
                      suffix = suffix,
                      pal.thresh = tree$frame$yval[2],
                      box.palette = if (is.function(palette)) palette(nrow(tree$frame)) else palette,
                      branch = 0.3,
                      branch.col = colors,
                      split.col = colors,
                      col = colors,
                      ...)

      grDevices::dev.flush()
      Sys.sleep(1)
    }
  } else {
    subtree <- subtree(tree, leaves)
    prefix <- c("ATE = ", rep("GATE = ", times = (length(subtree$frame$n) - 1)))
    sizes <- subtree$frame$n
    percentages <- round(subtree$frame$n / subtree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    rpart.plot::prp(subtree,
                    type = type,
                    extra = 0,
                    under = FALSE,
                    fallen.leaves = TRUE,
                    round = 0,
                    leaf.round = 0,
                    prefix = prefix,
                    suffix = suffix,
                    pal.thresh = subtree$frame$yval[2],
                    box.palette = if (!is.function(palette)) palette else palette(nrow(subtree$frame)),
                    branch = 0.3,
                    ...)
  }
}


#' Summary Method for aggTrees Objects
#'
#' Summarizes an \code{aggTrees} object.
#'
#' @param object \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
summary.aggTrees <- function(object, ...) {
  cat("Honest estimates:", object$honest_estimates, "\n")
  summary(object$tree)
}


#' Print Method for aggTrees Objects
#'
#' Prints an \code{aggTrees} object.
#'
#' @param x \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
print.aggTrees <- function(x, ...) {
  cat("Honest estimates:", x$honest_estimates, "\n")
  print(x$tree)
}
