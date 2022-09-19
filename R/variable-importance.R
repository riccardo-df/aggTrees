#' Variable Importance for Estimated CATEs
#'
#' Grows a regression forest using the estimated CATEs as the dependent variable,
#' and provides the covariates' relative importance.
#'
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' A data frame storing the relative importance of each covariate.
#'
#' @details
#' \code{var_importance} fits a regression forest using the covariates in \code{X} and \code{cates} as the
#' dependent variable. This allows to gain some knowledge on how covariates are related to treatment effect
#' heterogeneity.\cr
#' Variable importance is computed as a weighted sum of how many times a particular covariate was split on
#' at each depth in the forest. For more details, please see \code{\link[grf]{variable_importance}}.
#'
#' @import grf magrittr
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{plot_importance}}
#'
#' @export
var_importance <- function(cates, X) {
  ## Checks.
  if (!is.numeric(cates)) stop("'cates' must be a numeric vector.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)

  `%>%` <- magrittr::`%>%`

  ## Fitting forest.
  forest <- grf::regression_forest(X, cates, tune.parameters = "all")

  ## Extracting variable importance.
  var_df <- data.frame(varImp = grf::variable_importance(forest), name = colnames(X))

  ## Handling output.
  var_df <- var_df[order(var_df$varImp, decreasing = TRUE), ]
  var_df <- dplyr::tibble(Variable = var_df$name, Importance = as.numeric(var_df$varImp)) %>%
    dplyr::arrange(Importance) %>%
    dplyr::mutate(Variable = factor(Variable, levels = unique(Variable)))

  ## Output.
  return(var_df)
}


#' Variable Importance Plot
#'
#' Plots the relative variable importance.
#'
#' @param importance The output of \code{\link{var_importance}}.
#' @param k Display only the k most important covariates. All covariates are displayed by default.
#' @param title String. The title of the plot.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot2}} object.
#'
#' @import ggplot2
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{var_importance}}
#'
#' @export
plot_importance <- function(importance, k = 0, title = "") {
  ## Checks.
  if (k < 1 | k > dim(importance)[1]) stop("'k' must be in the interval [1, k_max], with k_max the number of covariates.", call. = FALSE)

  ## Selecting desired covariates.
  if (k == 0) idx <- 1:dim(importance)[1] else idx <- dim(importance)[1]:(dim(importance)[1]-k+1)

  ## Generating plot.
  plot <- ggplot2::ggplot(importance[idx, ], ggplot2::aes(x = Variable, y = Importance)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = .6, width = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ## Plotting.
  plot
}
