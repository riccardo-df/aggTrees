#' Cluster Plot
#'
#' Plots the average values of the covariates for each cluster in a two-dimensional covariate space. Useful to compare
#' clustering methods with aggregation trees.
#'
#' @param clusters The output of \code{\link{kmeans}} or \code{\link[T4cluster]{kmeanspp}}.
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept). The same used in \code{tree}.
#' @param plot Optional, the output of \code{\link{recursive_partitioning_plot}}. If provided, the cluster plot is overimposed to \code{plot}.
#' @param low String. Color to represent more negative \code{cates}. Ignored if \code{plot} is provided.
#' @param high String. Color to represent more positive \code{cates}. Ignored if \code{plot} is provided.
#' @param size Size of points.
#'
#' @return
#' The desired plot.
#'
#' @details
#' \code{cluster_plot} helps visualizing how clustering methods compress the information of the \code{cates}. \cr\cr
#' The plot is built as follows. First, a scatter plot of the two covariates is displayed. Each point is colored
#' according to the associated value of \code{cates}, so that it is immediate noticing where treatment effects
#' are stronger and lighter (colors can be specified by the user). Second, additional points are overimposed, one
#' for each cluster. The coordinates of these points correspond to the average value of \code{X[, 1]} and \code{X[, 2]}
#' for the associated cluster.
#'
#' @import rpart parttree
#'
#' @author Riccardo Di Francesco
#'
#' @seealso aggregation_tree, subtree, plot_tree
#'
#' @export
cluster_plot <- function(clusters, cates, X,
                         plot = NULL, low = "yellow", high = "red", size = 1) {
  ## Handling inputs and checks.
  if (dim(X)[2] > 2) stop("You are using more than two covariates. \n When more than two covariates are available, the recursive partitioning plot is not avaialable. \n You may consider a tree-plot (see ??plot_tree).")

  n_cluster <- length(unique(clusters$cluster))

  clusters_x1_coord <- numeric(n_cluster)
  clusters_x2_coord <- numeric(n_cluster)
  for (i in seq_len(n_cluster)) {
    clusters_x1_coord[i] <- mean(X[, 1][clusters$cluster == i])
    clusters_x2_coord[i] <- mean(X[, 2][clusters$cluster == i])
  }

  ## Plotting.
  if (is.null(plot)) {
    plot <- ggplot2::ggplot(data.frame("x1" = X[, 1], "x2" = X[, 2], "cates" = cates), ggplot2::aes(x = x1, y = x2)) +
      ggplot2::geom_point(ggplot2::aes(color = cates), size = size) +
      ggplot2::scale_color_gradient(low = low, high = high) +
      ggplot2::xlab(colnames(X)[1]) +
      ggplot2::ylab(colnames(X)[2]) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")
  }

  cluster_plot <- plot +
    ggplot2::geom_point(
      data = data.frame(x1 = clusters_x1_coord, x2 = clusters_x2_coord),
      colour = "black", fill = "black", alpha = 1, pch = 21, size = 8
    )

  ## Output.
  return(cluster_plot)
}
