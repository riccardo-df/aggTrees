#' Variable Importance for Estimated CATEs
#'
#' Grows a random forest using the estimated CATEs as dependent variable, and provides the covariates' relative 
#' variable importance.
#'
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' A data frame storing the relative importance of each covariate.
#'
#' @details
#' \code{cates} are used to grow a random forest using the the covariates in \code{X}. This allows to gain some knowledge
#' on how covariates are related to treatment effect heterogeneity. \cr\cr
#' Variable importance is computed as a weighted sum of how many times a particular covariate was split on
#' at each depth in the forest. For more details, please see \code{\link[grf]{variable_importance}}.
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
#' ## Computing importance using only aggregation sample.
#' cates <- predict(cates_forest, newdata = X_agg)$predictions
#'
#' var_imp <- var_importance(cates, X_agg)
#' plot_importance(var_imp)
#' 
#' ## Looking only at the ten most important covariates.
#' plot_importance(var_imp, k = 10)
#'
#' @export
var_importance <- function(cates, X) {
  ## Fitting forest.
  forest <- grf::regression_forest(X, cates, tune.parameters = "all")
  
  ## Extracting variable importance.
  var_df <- data.frame(varImp = grf::variable_importance(forest), name = colnames(X))
  
  ## Handling output.
  var_df <- var_df[order(var_df$varImp, decreasing = TRUE), ]
  
  `%>%` <- magrittr::`%>%`
  
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
#'
#' @return
#' None. It plots the relative importance of the covariates.
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
#' ## Computing importance using only aggregation sample.
#' cates <- predict(cates_forest, newdata = X_agg)$predictions
#'
#' var_imp <- var_importance(cates, X_agg)
#' plot_importance(var_imp)
#' 
#' ## Looking only at the ten most important covariates.
#' plot_importance(var_imp, k = 10)
#' 
#' @export
plot_importance <- function(importance, k = 0) {
  ## Selecting desired covariates.
  if (k == 0) idx <- 1:dim(importance)[1] else idx <- dim(importance)[1]:(dim(importance)[1]-k+1)
  
  ## Generating plot.
  plot <- ggplot2::ggplot(importance[idx, ], ggplot2::aes(x = Variable, y = Importance)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = .6, width = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::xlab("Labels") + ggplot2::ylab("Variable importance") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  ## Plotting.
  plot
}

