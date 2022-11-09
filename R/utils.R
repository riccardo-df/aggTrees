#' Sample Splitting
#'
#' Splits the sample in three different parts: estimation sample, aggregation sample, and honest sample.
#'
#' @param n Size of the sample to be split.
#' @param estimation_frac Fraction of units for the estimation sample.
#' @param aggregation_frac Fraction of units for the aggregation sample.
#' @param honest_frac Fraction of units for the honest sample.
#'
#' @return
#' A list storing the indexes for the three different subsamples.
#'
#' @details
#' The estimation and the aggregation sample are mandatory to use \code{\link{aggregation_tree}}, meaning that both
#' \code{estimation_frac} and \code{aggregation_frac} must be non-zero.\cr
#'
#' The honest sample is required only if the user wants to conduct valid inference, which generally comes at the price of
#' a larger mean squared error. If honesty is not desired, please set \code{honest_frac} equal to zero.
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @author Riccardo Di Francesco
#'
#' @export
sample_split <- function(n,
                         estimation_frac = 0.5, aggregation_frac = 0.25, honest_frac = 0.25) {
  ## Handling inputs and checks.
  if (n <= 0) stop("'n' cannot be equal or lower than zero.", call. = FALSE)
  if (estimation_frac <= 0 | estimation_frac >= 1) stop("'estimation_frac' must lie in the interval (0, 1).", call. = FALSE)
  if (aggregation_frac <= 0 | aggregation_frac >= 1) stop("'aggregation_frac' must lie in the interval (0, 1).", call. = FALSE)
  if (honest_frac < 0 | honest_frac > 1) stop("'honest_frac' must lie in the interval [0, 1].", call. = FALSE)
  if (estimation_frac + aggregation_frac + honest_frac != 1) stop("Fractions for the different samples do not sum up to one.", call. = FALSE)

  ## Split the sample.
  estimation_idx <- sample(1:n, floor(estimation_frac * n), replace = FALSE)
  aggregation_idx <- sample(setdiff(1:n, estimation_idx), floor(aggregation_frac * n), replace = FALSE)
  honest_idx <- setdiff(1:n, union(estimation_idx, aggregation_idx))

  ## Output.
  return(list("estimation_idx" = estimation_idx, "aggregation_idx" = aggregation_idx, "honest_idx" = honest_idx))
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
#' @seealso \code{\link{aggregation_tree}}, \code{\link{subtree_aggtree}}, \code{\link{plot.aggTrees}}
#'
#' @export
recursive_partitioning_plot <- function(tree, cates, X, low = "yellow", high = "red", size = 1) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)
  if (dim(X)[2] > 2) stop("You are using more than two covariates. You may consider a tree-plot (see ??plot_tree).", call. = FALSE)
  if (length((tree$ordered)) > 2) stop("'tree' has been grown with more than two covariates. You may consider a tree-plot (see ??plot_tree).", call. = FALSE)

  x1 <- NULL
  x2 <- NULL

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


#' Renaming Variables for LATEX Usage (Internal Use)
#'
#' Renames variables where the character "_" is used, which causes clashes in LATEX. Useful for the \code{phased} print method.
#'
#' @param names string vector.
#'
#' @return
#' The renamed string vector. Strings where "_" is not found are not modified by \code{rename_latex}.
rename_latex <- function(names) {
  ## Locating variables that need renaming.
  idx <- grepl("_", names, fixed = TRUE)

  ## Renaming variables.
  split_names <- stringr::str_split(string = names[idx], pattern = "_", simplify = TRUE)
  attach_names <- paste(split_names[, 1], split_names[, 2], sep = "\\_")

  ## Replacing.
  names[idx] <- attach_names

  ## Output.
  return(names)
}


#' Covariate Matrix Expansion
#'
#' Expands the covariate matrix, adding interactions and polynomials. This is particularly useful for penalized regressions.
#'
#' @param X Covariate matrix (no intercept).
#' @param int_order Order of interactions to be added. Set equal to one if no interactions are desired.
#' @param poly_order Order of the polynomials to be added. Set equal to one if no polynomials are desired.
#' @param threshold Drop binary variables representing less than \code{threshold}\% of the population. Useful to speed up computation.
#'
#' @return
#' The expanded covariate matrix, as a data frame.
#'
#' @details
#' \code{expand_df} assumes that categorical variables are coded as \code{factors}. Also, no missing values are allowed.
#'
#' @export
expand_df <- function(X, int_order = 2, poly_order = 4, threshold = 0) {
  ## Handling inputs and checks.
  if (int_order < 1 | int_order > 4) stop("Wrong order of interactions! Must be either 1, 2, 3 or 4.")
  if (poly_order < 1) stop("Wrong order of polynomials! Must be greater than zero.")
  if (threshold < 0 | threshold > 1) stop("Wrong threshold! Must lie in the open interval (0, 1).")
  if (threshold == 1) stop("Wrong threshold! Cannot drop all the units.")

  X <- as.data.frame(X)
  X <- stats::model.matrix(~., data = data.frame(X))[, -1]

  X_continuous <- X[, !apply(X, MARGIN = 2, function(x) all(x %in% 0:1))]

  ## Adding int_order-way interactions.
  if (int_order == 1) {
    expanded_X <- X
  } else if (int_order == 2) {
    expanded_X <- stats::model.matrix(~ .^2, data = data.frame(X))[, -1]
  } else if (int_order == 3) {
    expanded_X <- stats::model.matrix(~ .^3, data = data.frame(X))[, -1]
  } else if (int_order == 4) {
    expanded_X <- stats::model.matrix(~ .^4, data = data.frame(X))[, -1]
  }

  ## Adding polynomials for continuous variables (works on original continuous covariates).
  if (poly_order > 1) {
    for (i in seq_len(dim(X_continuous)[2]))
    {
      temp.poly <- stats::poly(X_continuous[, i], degree = poly_order, raw = TRUE)[, -1]
      expanded_X <- data.frame(expanded_X, temp.poly)

      for (j in 2:poly_order)
      {
        colnames(expanded_X)[(dim(expanded_X)[2]) - poly_order + j] <- paste(paste(colnames(X_continuous)[i], "..", sep = ""), j, sep = "")
      }
    }
  }

  ## Dropping binary variables with low variability.
  # Bit tricky: the following is an index which equals TRUE iff the column of expanded_X is binary and with low variability.
  temp_idx <- apply(expanded_X, MARGIN = 2, function(x) all(x %in% 0:1)) * (apply(expanded_X, MARGIN = 2, mean) < threshold)
  expanded_X <- expanded_X[, !temp_idx]

  ## Handling output.
  out <- data.frame(expanded_X)

  ## Output.
  return(out)
}
