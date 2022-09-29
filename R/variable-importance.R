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
#' at each depth in the forest. For more details, please see \code{\link[grf]{variable_importance}}.
#'
#' @import grf magrittr
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}
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
