#' Tree-Plot
#'
#' Plots a \code{\link[rpart]{rpart}} object.
#'
#' @param tree A  \code{\link[rpart]{rpart}} object.
#' @param leaves Number of leaves of the desired tree. This can be used to plot subtrees of \code{tree}.
#' @param palette Palette to color the nodes. This must be a vector of colors.
#' @param ... Further arguments for \code{\link[rpart.plot]{prp}}. See the details section for suggested parameters.
#'
#' @return
#' None. It plots a tree.
#'
#' @details
#' Please refer to \code{\link[colorspace]{choose_palette}} for high-quality palettes.\cr
#'
#' Among the most useful parameters to be passed in \code{...} there is:
#'   \item{\code{main}}{String for the title of the plot.}
#'   \item{\code{cex.main}}{Size of the title.}
#'   \item{\code{split.cex}}{Size of splitting labels.}
#'
#' @import rpart.plot
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree}}, \code{\link{plot_tree_sequence}}, \code{\link{recursive_partitioning_plot}}
#'
#' @export
plot_tree <- function(tree, leaves = get_leaves(tree), palette = "auto", ...) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  subtree <- subtree(tree, leaves)
  labels <- c("ATE \n", rep("GATE \n", times = (length(subtree$frame$n) - 1)))

  ## Plotting.
  rpart.plot::prp(subtree,
    type = 2,
    extra = 101,
    under = FALSE,
    fallen.leaves = TRUE,
    round = 0,
    leaf.round = 0,
    prefix = labels,
    box.palette = palette(nrow(subtree$frame)),
    branch = 0.3,
    tweak = 1,
    ...)
}


#' Tree-Plot sequence
#'
#' Plots a sequence of \code{\link[rpart]{rpart}} subtrees.
#'
#' @param tree A  \code{\link[rpart]{rpart}} object.
#' @param palette Palette to color the nodes. this must be a vector of colors.
#' @param ... Further arguments for \code{\link[rpart.plot]{prp}}. See the details section for suggested parameters.
#'
#' @return
#' None. It plots a sequence of subtree.
#'
#' @details
#' Please refer to \code{\link[colorspace]{choose_palette}} for high-quality palettes.\cr
#'
#' Among the most useful parameters to be passed in \code{...} there is:
#' \describe{
#'   \item{\code{main}}{String for the title of the plot.}
#'   \item{\code{cex.main}}{Size of the title.}
#'   \item{\code{split.cex}}{Size of splitting labels.}
#' }
#'
#' @import rpart rpart.plot
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree}}, \code{\link{plot_tree}}, \code{\link{recursive_partitioning_plot}}
#'
#' @export
plot_tree_sequence <- function(tree, palette = "auto", ...) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)
  labels <- c("ATE \n", rep("GATE \n", times = (length(tree$frame$n) - 1)))

  full_nodes <- rownames(tree$frame) # Storing rownames of full tree (they serve as node labels).
  alpha_values <- rev(tree$cptable[, "CP"]) # Threshold values of cost-complexity parameter.

  ## Plotting the sequence.
  for (alpha in alpha_values) {
    temp_nodes <- rownames(rpart::prune.rpart(tree, alpha)$frame) # Node labels of pruned tree.
    colors <- ifelse(full_nodes %in% temp_nodes, 1, "white") # Graying out collapsed nodes.

    grDevices::dev.hold()

    rpart.plot::prp(tree,
      type = 2,
      extra = 101,
      under = FALSE,
      fallen.leaves = TRUE,
      round = 0,
      leaf.round = 0,
      prefix = labels,
      box.palette = palette(nrow(tree$frame)),
      branch = 0.3,
      tweak = 1,
      branch.col = colors,
      split.col = colors,
      col = colors,
      ...)

    grDevices::dev.flush()
    Sys.sleep(0.05)
  }
}


#' Recursive Partitioning Plots
#'
#' Plots the recursive partitioning of \code{\link[rpart]{rpart}} objects for a two-dimensional covariate space.
#'
#' @param tree A \code{\link[rpart]{rpart}} object. The tree must have been built using only two covariates.
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept), the same used to construct \code{tree}.
#' @param low String. Color to represent more negative \code{cates}.
#' @param high String. Color to represent more positive \code{cates}.
#' @param size Size of points in the scatter plot.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot2}} object.
#'
#' @details
#' The plot is built as follows. First, a scatter plot of the two covariates is displayed. Each point is colored
#' according to the associated value of \code{cates}, so that it is immediate noticing where treatment effects
#' are stronger and lighter (colors can be specified by the user). Second, the axis-aligned splits of \code{tree}
#' are overimposed.
#'
#' @import ggplot2 parttree
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree}}, \code{\link{plot_tree}}
#'
#' @export
recursive_partitioning_plot <- function(tree, cates, X, low = "yellow", high = "red", size = 1) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)
  if (dim(X)[2] > 2) stop("You are using more than two covariates. You may consider a tree-plot (see ??plot_tree).", call. = FALSE)
  if (length((tree$ordered)) > 2) stop("'tree' has been grown with more than two covariates. You may consider a tree-plot (see ??plot_tree).", call. = FALSE)

  ## Plotting.
  plot <- ggplot2::ggplot(data.frame("x1" = X[, 1], "x2" = X[, 2], "cates" = cates), ggplot2::aes(x = x1, y = x2)) +
    ggplot2::geom_point(ggplot2::aes(color = cates), size = size) +
    ggplot2::scale_color_gradient(low = low, high = high) +
    ggplot2::xlab(colnames(X)[1]) +
    ggplot2::ylab(colnames(X)[2]) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")

  tree_plot <- plot +
    parttree::geom_parttree(data = tree, show.legend = FALSE, flipaxes = FALSE, ggplot2::aes(fill = cates), alpha = 0, colour = "black")

  ## Output.
  return(tree_plot)
}
