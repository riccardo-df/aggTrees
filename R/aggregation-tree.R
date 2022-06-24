#' Aggregation Trees
#'
#' Grows an aggregation tree.
#'
#' @param cates CATEs vector.
#' @param X Covariate matrix (no intercept).
#' @param maxdepth The maximum depth (i.e, the maximum number of nodes connecting the root to the leaves) of the tree.
#' @param cp Minimum MSE increase to accept splits.
#'
#' @return
#' A \code{\link[rpart]{rpart}} object.
#'
#' @details
#' For details about aggregation trees and their use, please see Di Francesco (2022). \cr
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
aggregation_tree <- function(cates, X, maxdepth, cp) {
  ## Growing tree.
  tree <- rpart::rpart(cates ~ .,
    data = data.frame("cates" = cates, X), method = "anova",
    control = rpart::rpart.control(maxdepth = maxdepth, cp = cp), model = TRUE
  )

  ## Output.
  return(tree)
}
