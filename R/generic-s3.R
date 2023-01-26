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
#' \code{\link{build_aggtree}}, \code{\link{analyze_aggtree}}
#'
#' @export
plot.aggTrees <- function(x, leaves = get_leaves(x$tree),
                          type = 2, palette = c("BuRd"),
                          sequence = FALSE, ...) {
  ## Handling inputs and checks.
  if (!(inherits(x, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(sequence %in% c(TRUE, FALSE))) stop("Invalid 'sequence'. This must be either TRUE or FALSE.", call. = FALSE)

  tree <- x$tree

  ## Plotting.
  if (sequence) {
    prefix <- c("ATE = ", rep("GATE = ", times = (length(tree$frame$n) - 1)))
    sizes <- tree$frame$n
    percentages <- round(tree$frame$n / tree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    nodes <- as.numeric(rownames(tree$frame)) # Node numbers of full tree.
    # parents <- full_nodes %/% 2 # Parent of each node.
    alpha_values <- rev(tree$cptable[, "CP"]) # Threshold values of cost-complexity parameter.

    for (alpha in alpha_values) {
      temp_nodes <- as.numeric(rownames(rpart::prune.rpart(tree, alpha)$frame)) # Node numbers of subtree.
      # temp_parents <- temp_nodes %/% 2 # Parent numbers of subtree.
      # split_colors <- ifelse(parents %in% c(0, temp_nodes), 1, "white")
      colors <- ifelse(nodes %in% temp_nodes, 1, "white") # Graying out pruned leaves.

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
                      pal.thresh = tree$frame$yval[1],
                      box.palette = if (is.function(palette)) palette(nrow(tree$frame)) else palette,
                      branch = 0.3,
                      branch.col = colors,
                      split.col = "black",
                      col = colors)
                      # ...)

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
                    pal.thresh = tree$frame$yval[1],
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
#' \code{\link{build_aggtree}}, \code{\link{analyze_aggtree}}
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
  if (is.null(object$idx$honest_idx)) cat("Honest estimates:", FALSE, "\n") else cat("Honest estimates:", TRUE, "\n")
  summary(object$tree)
}


#' Print Method for aggTrees Objects
#'
#' Prints an \code{aggTrees} object.
#'
#' @param x \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' \code{\link{build_aggtree}}, \code{\link{analyze_aggtree}}
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
  if (is.null(x$idx$honest_idx)) cat("Honest estimates:", FALSE, "\n") else cat("Honest estimates:", TRUE, "\n")
  print(x$tree)
}
