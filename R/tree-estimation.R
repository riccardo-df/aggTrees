#' Estimating Leaf-Effects via Linear Models
#'
#' Uses the leaves of a tree stored in an \code{\link[rpart]{rpart}} object to estimate a linear model via OLS. The
#' estimated coefficients correspond to the GATEs in each leaf. If the data used in the OLS estimation have not been
#' used to grow the tree (a condition called "honesty"), then one can use the standard errors for the tree's estimates
#' to construct valid confidence intervals.
#'
#' @param tree A \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param D Treatment assignment vector
#' @param X Covariate matrix (no intercept).
#' @param method Either \code{"raw"} or \code{"aipw"}, defines the outcome used in the regression.
#'
#' @return
#' The fitted model, as an \code{\link[estimatr]{lm_robust}} object.
#'
#' @details
#' The \code{method} argument controls how GATEs are estimated. If \code{"method" == "raw"}, we estimate via OLS the following
#' linear model:
#'
#' \deqn{Y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th leaf of the tree, and \code{|T|} the
#' number of groups. If the treatment is randomly assigned, one can show that the estimated betas identify the GATEs of
#' each group. Thus, we can interpret the OLS results as usual. However, in observational studies these estimates are biased
#' due to selection into treatment. To get unbiased estimates, we can set \code{"method"} to \code{"aipw"} to construct
#' doubly-robust scores \code{y_i^*} and use them as a pseudo-outcome in the following regression:
#'
#' \deqn{Y_i^* = \sum_{l = 1}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' This way, we get unbiased GATEs estimates, and we can again interpret OLS results as usual. Scores are
#' constructed via 5-fold cross fitting. Honest regression forests are used to estimate the propensity score and the
#' conditional mean function of the outcome.\cr
#'
#' Notice that "honesty" is a necessary requirement to get valid inference. Thus, observations in \code{y}, \code{X}, and
#' \code{D} must not have been used to grow the \code{tree}.\cr
#'
#' Regardless of \code{method}, standard errors are estimated via the Eicker-Huber-White estimator.\cr
#'
#' If \code{tree} consists of a root only, \code{causal_ols_rpart} regresses \code{y} or the doubly-robust scores on a constant
#' and \code{D}, thus estimating the ATE.
#'
#' @import rpart estimatr causalDML
#' @importFrom stats predict
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item S Athey, G Imbens (2016). Recursive partitioning for heterogeneous causal effects. Proceedings of the National Academy of Sciences. \doi{10.1073/pnas.1510489113}.
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#'   \item J Robins, A Rotnitzky (1995). Semiparametric efficiency in multivariate regression models with missing data. \doi{10.2307/2291135}.
#'   \item V Semenova, V Chernozhukov (2021). Debiased machine learning of conditional average treatment effects and other causal functions. \doi{10.1093/ectj/utaa027}.
#' }
#'
#' @seealso \code{\link{estimate_rpart}}
#'
#' @export
causal_ols_rpart <- function(tree, y, X, D, method = "aipw") {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)
  if(!(method %in% c("raw", "aipw"))) stop("Invalid 'method'. It must be either 'raw' or 'aipw'.", call. = FALSE)

  ## Generate leaves indicators.
  leaves <- leaf_membership(tree, X)

  if (length(unique(leaves)) < get_leaves(tree)) warning("One or more leaves are empty: No observations in X fall there.")

  ## OLS estimation.
  if (method == "raw") {
    if (length(unique(leaves)) == 1) {
      model <- estimatr::lm_robust(y ~ D, data = data.frame("y" = y, "D" = D), se_type = "HC1")
    } else {
      model <- estimatr::lm_robust(y ~ 0 + leaf + D:leaf, data = data.frame("y" = y, "leaf" = leaves, "D" = D), se_type = "HC1")
    }
  } else if (method == "aipw") {
    aipw <- causalDML::DML_aipw(y, D, X)
    scores <- aipw$ATE$delta

    if (length(unique(leaves)) == 1) {
      model <- estimatr::lm_robust(scores ~ 1, data = data.frame("scores" = scores), se_type = "HC1")
    } else {
      model <- estimatr::lm_robust(scores ~ 0 + leaf, data = data.frame("scores" = scores, "leaf" = leaves), se_type = "HC1")
    }
  }

  ## Output.
  return(model)
}


#' Estimation of rpart Objects
#'
#' Replace leaf predictions of an \code{\link[rpart]{rpart}} object using external data.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param D Treatment assignment vector.
#' @param X Covariate matrix (no intercept).
#' @param method Either \code{"raw"} or \code{"cates"}, controls how leaf predictions are replaced.
#'
#' @return
#' A tree with leaf predictions replaced, as an \code{\link[rpart]{rpart}} object.
#'
#' @details
#' \code{estimate_rpart} replaces the leaf estimates of a tree as follows. First, it "pushes" all observations in
#' \code{X} down the tree and finds the leaves where they fall. Then, it replaces the predictions in each leaf according
#' to the user-specified \code{method}.
#'
#' If \code{method = "raw"}, \code{estimate_rpart} replaces leaf predictions with the difference between the sample average
#' of the observed outcomes of treated units and the sample average of the observed outcomes of control units in each leaf,
#' which is an unbiased estimator of the GATEs if the assignment to treatment is randomized.\cr
#'
#' If \code{method = "aipw"}, \code{estimate_rpart} replaces leaf predictions with the sample average of doubly-robust
#' scores in each leaf. This is a valid estimator of the GATEs in observational studies. Scores are constructed via
#' 5-fold cross fitting. Honest regression forests are used to estimate the propensity score and the conditional mean
#' function of the outcome.\cr
#'
#' \code{estimate_rpart} allows the user to implement "honest" estimation. If observations in \code{X} have not been
#' used to construct the tree, then the new predictions are honest in the sense of Athey and Imbens (2016). This allows
#' the user to conduct valid inference about the estimated GATEs with standard approaches, e.g., by constructing conventional
#' confidence intervals. To get standard errors for the tree's estimates, please use \code{\link{causal_ols_rpart}}.\cr
#'
#' Due to coding limitations, \code{estimate_rpart} replaces the estimates only in the leaves of a tree. The
#' internal nodes' estimates are not replaced, so the user should ignore them and focus on the leaves of the tree.\cr
#'
#' @import treeClust Rcpp causalTree causalDML
#' @useDynLib aggTrees
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item S Athey, G Imbens (2016). Recursive partitioning for heterogeneous causal effects. Proceedings of the National Academy of Sciences. \doi{10.1073/pnas.1510489113}.
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#'   \item J Robins, A Rotnitzky (1995). Semiparametric efficiency in multivariate regression models with missing data. \doi{10.2307/2291135}.
#'   \item V Semenova, V Chernozhukov (2021). Debiased machine learning of conditional average treatment effects and other causal functions. \doi{10.1093/ectj/utaa027}.
#' }
#'
#' @seealso \code{\link{causal_ols_rpart}}
#'
#' @export
estimate_rpart <- function(tree, X, y, D, method = "aipw") {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  if (!(method %in% c("raw", "aipw"))) stop("You must provide a valid method.", call. = FALSE)

  new_tree <- tree

  if (method == "raw") {
    new_tree <- causalTree::estimate.causalTree(tree, data.frame(X, y), treatment = D)
  } else if (method == "aipw") {
    ## Extract leaves where observations in X fall.
    leaves <- treeClust::rpart.predict.leaves(tree, data.frame(X), type = "where") # Row numbers of tree$frame.
    names(leaves) <- 1:length(y)
    unique_leaves <- unique(leaves)

    ## Construct scores.
    aipw <- causalDML::DML_aipw(y, D, X)
    scores <- aipw$ATE$delta

    ## Call cpp to compute honest leaf estimates (honest only if y, and X belong to honest sample).
    honest_estimates <- as.matrix(honest_rpart_cpp(unique_leaves, scores, leaves))

    ## Replace leaf estimates.
    for (leaf in unique(leaves)) {
      new_tree$frame$yval[leaf] <- honest_estimates[as.numeric(names(leaves)[leaves == leaf][1])]
      new_tree$frame$n[leaf] <- sum(leaves == leaf)
    }
  }

  ## Output.
  return(new_tree)
}
