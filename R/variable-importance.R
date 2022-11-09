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
#'
#' Variable importance is computed as a weighted sum of how many times a particular covariate was split on
#' at each depth in the forest. For more details, please see \code{\link[grf]{variable_importance}}.\cr
#'
#' Notice that variable importance is generally a rough diagnostic of how treatment effects relate to the
#' covariates. If two covariates are highly correlated, trees generally split on only one of those covariates, thus
#' assigning low importance to the other. Therefore, one should not conclude that covariates with low importance are
#' not related to heterogeneity. A more systematic way to assess how treatment effects relate to the covariates consists
#' of investigating how the average characteristics of the units vary across subpopulations that differ in the magnitude
#' of their treatment effects (see \code{\link{avg_characteristics_aggtree}}).
#'
#' @import grf magrittr
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{avg_characteristics_aggtree}}
#'
#' @export
var_importance <- function(cates, X) {
  ## Checks.
  if (!is.numeric(cates)) stop("'cates' must be a numeric vector.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)

  `%>%` <- magrittr::`%>%`

  Importance <- NULL
  Variable <- NULL

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
