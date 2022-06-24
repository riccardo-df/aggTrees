#' Tree-Plot for Aggregation Trees
#'
#' Plots an aggregation tree via the tree representation.
#'
#' @param aggregation_tree The output of \code{\link{aggregation_tree}}.
#' @param palette Palette to color the nodes.
#' @param sequence Logical. If \code{FALSE} (default value), the entire tree is displayed. Otherwise, an animation showing the whole sequence of subtrees is provided.
#' @param ... Additional arguments for \code{\link[rpart.plot]{prp}}. See the details section for suggested parameters.
#'
#'
#' @return
#' None. It plots an aggregation tree.
#'
#' @details
#' \code{palette} helps visualizing the information from the \code{aggregation_tree}. As suggested in Di Francesco (2022),
#' using diverging colors (e.g., red and blue) to indicate effects stronger and lighter than the average treatment effect
#' gives an immediate indication of the most and least impacted groups. \cr\cr
#' If palette is not passed, \code{plot_aggregation_tree} asks the user to select one. See the example to understand how
#' to avoid being asked a palette each time. \cr\cr
#' Among the most useful parameters to be passed in \code{...} there is:
#' \describe{
#'   \item{\code{main}}{String for the title of the plot.}
#'   \item{\code{cex.main}}{Size of the title.}
#'   \item{\code{split.cex}}{Size of splitting labels.}
#'   }
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
plot_aggregation_tree <- function(aggregation_tree,
                                  palette = NULL, sequence = FALSE, ...) {
  ## Handling inputs and checks.
  if (is.null(palette)) palette <- colorspace::choose_palette()
  labels <- c("ATE \n", rep("GATE \n", times = (length(aggregation_tree$frame$n) - 1)))

  ## Plotting.
  if (!sequence) {
    rpart.plot::prp(aggregation_tree,
      type = 2,
      extra = 101,
      under = FALSE,
      fallen.leaves = TRUE,
      round = 0,
      leaf.round = 0,
      prefix = labels,
      box.palette = palette(nrow(aggregation_tree$frame)),
      branch = 0.3,
      tweak = 1,
      ...
    )
  } else {
    full_nodes <- rownames(aggregation_tree$frame) # Storing rownames of full tree (they serve as node labels).
    alpha_values <- rev(aggregation_tree$cptable[, "CP"]) # Threshold values of cost-complexity parameter.

    for (alpha in alpha_values) {
      temp_nodes <- rownames(rpart::prune.rpart(aggregation_tree, alpha)$frame) # Node labels of pruned tree.
      colors <- ifelse(full_nodes %in% temp_nodes, 1, "gray90") # Graying out collapsed nodes.

      grDevices::dev.hold()
      rpart.plot::prp(aggregation_tree,
        type = 2,
        extra = 101,
        under = FALSE,
        fallen.leaves = TRUE,
        round = 0,
        leaf.round = 0,
        prefix = labels,
        box.palette = palette(nrow(aggregation_tree$frame)),
        branch = 0.3,
        tweak = 1,
        branch.col = colors,
        split.col = colors,
        col = colors,
        ...
      )
      grDevices::dev.flush()
      Sys.sleep(0.05)
    }
  }
}


#' Recursive Partitioning Plot for Aggregation Trees
#'
#' Plots the recursive partitioning of an aggregation tree for a two-dimensional covariate space.
#'
#' @param aggregation_tree The output of \code{\link{aggregation_tree}}. The tree must have been built using two covariates only.
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept). The same used in \code{aggregation_tree}.
#' @param low String. Color to represent more negative \code{cates}.
#' @param high String. Color to represent more positive \code{cates}.
#'
#' @return
#' The desired plot.
#'
#' @details
#' \code{recursive_partitioning_plot} helps visualizing the recursive partitioning of \code{aggregation_tree} when two
#' covariates only are used to build the tree. \cr\cr
#' The plot is built as follows. First, a scatter plot of the two covariates is displayed. Each point is colored
#' according to the associated value of \code{cates}, so that it is immediate noticing where treatment effects
#' are stronger and lighter (colors can be specified by the user). Second, the axis-aligned splits of \code{aggregation_tree}
#' are overimposed, thereby showing the different groups.
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
#' ## Growing tree with two covariates only using only the aggregation sample.
#' X_subsample <- as.data.frame(X_agg[, c("mage", "nprenatal")])
#' cates <- predict(cates_forest, newdata = X_agg)$predictions
#'
#' tree <- aggregation_tree(cates, X_subsample, maxdepth = 3, cp = 0.01)
#'
#' ## Plotting.
#' plot <- recursive_partitioning_plot(tree, cates, X_subsample)
#' plot
#'
#' @export
recursive_partitioning_plot <- function(aggregation_tree, cates, X,
                                        low = "yellow", high = "red") {
  ## Handling inputs and checks.
  if (dim(X)[2] > 2) stop("You are using more than two covariates. \n When more than two covariates are available, the recursive partitioning plot is not avaialable. \n You may consider a tree-plot (see ??plot_aggregation_tree).")
  if (length((aggregation_tree$ordered)) > 2) stop("The tree has been grown with more than two covariates.")

  ## Plotting.
  plot <- ggplot2::ggplot(data.frame("x1" = X[, 1], "x2" = X[, 2], "cates" = cates), ggplot2::aes(x = x1, y = x2)) +
    ggplot2::geom_point(ggplot2::aes(color = cates)) +
    ggplot2::scale_color_gradient(low = low, high = high) +
    ggplot2::xlab(colnames(X)[1]) +
    ggplot2::ylab(colnames(X)[2]) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")

  tree_plot <- plot +
    parttree::geom_parttree(
      data = aggregation_tree, show.legend = FALSE, flipaxes = FALSE,
      ggplot2::aes(fill = cates), alpha = 0, colour = "black"
    )

  ## Output.
  return(tree_plot)
}


#' Cluster Plot
#'
#' Plots the average values of the covariates for each cluster in a two-dimensional covariate space. Useful to compare
#' clustering methods with aggregation trees.
#'
#' @param clusters The output of \code{\link{kmeans}} or \code{\link[T4cluster]{kmeanspp}}.
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept). The same used in \code{aggregation_tree}.
#' @param plot Optional, the output of \code{\link{recursive_partitioning_plot}}. If provided, the cluster plot is overimposed to \code{plot}.
#' @param low String. Color to represent more negative \code{cates}. Ignored if \code{plot} is provided.
#' @param high String. Color to represent more positive \code{cates}. Ignored if \code{plot} is provided.
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
#' ## Clustering using only aggregation sample.
#' cates <- predict(cates_forest, newdata = X_agg)$predictions
#' clusters <- kmeans(cates, 3)
#'
#' ## Plotting.
#' X_subsample <- as.data.frame(X_agg[, c("mage", "nprenatal")])
#' plot <- cluster_plot(clusters, cates, X_subsample)
#' plot
#'
#' ## Overlaying aggregation tree.
#' tree <- aggregation_tree(cates, X_subsample, maxdepth = 3, cp = 0.01)
#' tree_plot <- recursive_partitioning_plot(tree, cates, X_subsample)
#' cluster_plot <- cluster_plot(clusters, cates, X_subsample, plot = tree_plot)
#' cluster_plot
#'
#' @export
cluster_plot <- function(clusters, cates, X,
                         plot = NULL, low = "yellow", high = "red") {
  ## Handling inputs and checks.
  if (dim(X)[2] > 2) stop("You are using more than two covariates. \n When more than two covariates are available, the recursive partitioning plot is not avaialable. \n You may consider a tree-plot (see ??plot_aggregation_tree).")

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
      ggplot2::geom_point(ggplot2::aes(color = cates)) +
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
