#' Subtree
#'
#' Extracts a subtree with a user-specified number of leaves from an \code{\link[rpart]{rpart}} object.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param leaves Number of leaves of the desired subtree.
#' @param cv If \code{TRUE}, \code{leaves} is ignored and a cross-validation criterion is used to select a partition.
#'
#' @return
#' The subtree, as an \code{\link[rpart]{rpart}} object.
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 3000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#'
#' y <- exp(X[, 1]) + 2 * X[, 2] * X[, 2] > 0 + rnorm(n)
#'
#' ## Construct tree.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame(y, X), cp = 0)
#'
#' ## Extract subtree.
#' sub_tree <- subtree(tree, leaves = 4)
#' sub_tree_cv <- subtree(tree, cv = TRUE)
#'
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{get_leaves}} \code{\link{node_membership}} \code{\link{leaf_membership}}
#'
#' @export
subtree <- function(tree, leaves = NULL, cv = FALSE) {
  ## Handling inputs and checks.
  if(!(inherits(tree, "rpart"))) stop("Invalid 'tree'. This must be an aggTrees object.", call. = FALSE)
  if (!(cv %in% c(TRUE, FALSE))) stop("Invalid 'cv'. must be either TRUE or FALSE.", call. = FALSE)
  if (is.null(leaves) & cv == FALSE) stop("Invalid combination of 'leaves' and 'cv'. Please specify a number of leaves or select the cross-validation option.", call. = FALSE)
  if (!is.null(leaves)) {
    if (leaves < 1) stop("Invalid 'leaves'. This must be a positive number.", call. = FALSE)
    if (leaves > get_leaves(tree)) stop("'leaves' is greater than the number of leaves of 'tree'. Please provide a deeper 'tree'.", call. = FALSE)
  }

  ## Output.
  if (cv) return(rpart::prune(tree, tree$cptable[, 1][which.min(tree$cptable[, 4])])) else return(rpart::prune(tree, tree$cptable[tree$cptable[, "nsplit"] == leaves - 1, "CP"]))
}


#' Number of Leaves
#'
#' Extracts the number of leaves of an \code{\link[rpart]{rpart}} object.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#'
#' @return
#' The number of leaves.
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 3000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#'
#' y <- exp(X[, 1]) + 2 * X[, 2] * X[, 2] > 0 + rnorm(n)
#'
#' ## Construct tree.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame(y, X))
#'
#' ## Extract number of leaves.
#' n_leaves <- get_leaves(tree)
#' n_leaves
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{subtree}} \code{\link{node_membership}} \code{\link{leaf_membership}}
#'
#' @export
get_leaves <- function(tree) {
  if (!inherits(tree, "rpart")) stop("Invalid 'tree'. This must be an rpart object.")

  return(dim(tree$frame[tree$frame$var == "<leaf>", ])[1])
}


#' Node Membership
#'
#' Constructs a binary variable that encodes whether each observation falls into a particular node of an
#' \code{\link[rpart]{rpart}} object.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param X Covariate matrix (no intercept).
#' @param node Number of node.
#'
#' @return
#' Logical vector denoting whether each observation in \code{X} falls into \code{node}.
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 3000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#'
#' y <- exp(X[, 1]) + 2 * X[, 2] * X[, 2] > 0 + rnorm(n)
#'
#' ## Construct tree.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame(y, X))
#'
#' ## Extract number of leaves.
#' is_in_third_node <- node_membership(tree, X, 3)
#' head(is_in_third_node)
#'
#' @importFrom utils capture.output
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{subtree}} \code{\link{leaf_membership}} \code{\link{get_leaves}}
#'
#' @export
node_membership <- function(tree, X, node) {                 # Taken from https://stackoverflow.com/questions/36748531/getting-the-observations-in-a-rparts-node-i-e-cart.
  invisible(capture.output(path <- path.rpart(tree, node)))  # Asked by Tal Galili, answered by DatamineR.
  rule <- sapply(path[[1]][-1], function(x) strsplit(x, '(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])', perl = TRUE))
  idx <- apply(do.call(cbind, lapply(rule, function(x) eval(call(x[2], X[, x[1]], as.numeric(x[3]))))), 1, all)
  return(idx)
}


#' Leaf Membership
#'
#' Constructs a variable that encodes in which leaf of an \code{\link[rpart]{rpart}} object the units in a given data frame fall.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' A factor whose levels denote in which leaf each unit falls. Leaves are ordered in increasing order of their predictions
#' (from most negative to most positive).
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 3000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#'
#' y <- exp(X[, 1]) + 2 * X[, 2] * X[, 2] > 0 + rnorm(n)
#'
#' ## Construct tree.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame(y, X))
#'
#' ## Extract number of leaves.
#' leaves_factor <- leaf_membership(tree, X)
#' head(leaves_factor)
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{subtree}} \code{\link{node_membership}} \code{\link{get_leaves}}
#'
#' @export
leaf_membership <- function(tree, X) {
  if (!inherits(tree, "rpart")) stop("Invalid 'tree'. This must be an rpart object.")
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)

  ## Inspired by https://bookdown.org/halflearned/tutorial/hte1.html.
  tree_predictions <- predict(tree, data.frame(X))
  n_leaves <- length(unique(tree_predictions))
  leaves <- factor(tree_predictions, levels = sort(unique(tree_predictions)), labels = seq(n_leaves))

  return(leaves)
}


#' GATE Estimation with rpart Objects
#'
#' Replaces node predictions of an \code{\link[rpart]{rpart}} object using external data to estimate the group average treatment
#' effects (GATEs).
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param D Treatment assignment vector.
#' @param X Covariate matrix (no intercept).
#' @param method Either \code{"raw"} or \code{"aipw"}, controls how node predictions are replaced.
#' @param scores Optional, vector of scores to be used in replacing node predictions. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#'
#' @return
#' A tree with node predictions replaced, as an \code{\link[rpart]{rpart}} object, and the scores (if \code{method == "raw"},
#' this is \code{NULL}).
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Split the sample.
#' splits <- sample_split(length(y), training_frac = 0.5)
#' training_idx <- splits$training_idx
#' honest_idx <- splits$honest_idx
#'
#' y_tr <- y[training_idx]
#' D_tr <- D[training_idx]
#' X_tr <- X[training_idx, ]
#'
#' y_hon <- y[honest_idx]
#' D_hon <- D[honest_idx]
#' X_hon <- X[honest_idx, ]
#'
#' ## Construct a tree using training sample.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame("y" = y_tr, X_tr), maxdepth = 2)
#'
#' ## Estimate GATEs in each node (internal and terminal) using honest sample.
#' new_tree <- estimate_rpart(tree, y_hon, D_hon, X_hon, method = "raw")
#' new_tree$tree
#'
#' @details
#' If \code{method == "raw"}, \code{estimate_rpart} replaces node predictions with the differences between the sample average
#' of the observed outcomes of treated units and the sample average of the observed outcomes of control units in each node,
#' which is an unbiased estimator of the GATEs if the assignment to treatment is randomized.\cr
#'
#' If \code{method == "aipw"}, \code{estimate_rpart} replaces node predictions with sample averages of doubly-robust
#' scores in each node. This is a valid estimator of the GATEs in observational studies. Honest regression forests
#' and 5-fold cross fitting are used to estimate the propensity score and the conditional mean function of the outcome
#' (unless the user specifies the argument \code{scores}).\cr
#'
#' \code{estimate_rpart} allows the user to implement "honest" estimation. If observations in \code{y}, \code{D} and \code{X}
#' have not been used to construct the \code{tree}, then the new predictions are honest in the sense of Athey and Imbens (2016).
#' To get standard errors for the tree's estimates, please use \code{\link{causal_ols_rpart}}.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso \code{\link{causal_ols_rpart}} \code{\link{avg_characteristics_rpart}}
#' @export
estimate_rpart <- function(tree, y, D, X, method = "aipw", scores = NULL) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'Invalid 'tree'. This must be an rpart object.", call. = FALSE)
  if (!(method %in% c("raw", "aipw"))) stop("Invalid 'method'. This must be either 'raw' or 'aipw'.", call. = FALSE)

  leaves <- leaf_membership(tree, X)
  n_leaves <- get_leaves(tree)

  if (length(unique(leaves)) < n_leaves) warning("One or more leaves are empty: No observations in 'X' falls there.")

  for (leaf in seq_len(n_leaves)) {
    if (method == "raw" & length(unique(D[leaves == leaf])) == 1) stop("One or more leaves contain only treated or control observations. This is incompatible with 'method = raw'.", call. = FALSE)
  }

  ## Replace node predictions.
  new_tree <- tree

  if (method == "raw") {
    new_tree$frame$yval[1] <- mean(y[D == 1]) - mean(y[D == 0])

    nodes <- as.numeric(rownames(new_tree$frame))[-1]
    counter <- 2
    for (node in nodes) {
      idx <- node_membership(new_tree, X, node)
      new_tree$frame$yval[counter] <- mean(y[idx][D[idx] == 1]) - mean(y[idx][D[idx] == 0])
      counter <- counter + 1
    }
  } else if (method == "aipw") {
    if (is.null(scores)) {
      scores <- dr_scores(y, D, X)
    }

    new_tree$frame$yval[1] <- mean(scores)

    nodes <- as.numeric(rownames(new_tree$frame))[-1]
    counter <- 2
    for (node in nodes) {
      idx <- node_membership(new_tree, X, node)
      new_tree$frame$yval[counter] <- mean(scores[idx])
      counter <- counter + 1
    }
  }

  ## Output.
  return(list("tree" = new_tree, "scores" = scores))
}


#' Estimation and Inference about the GATEs with rpart Objects
#'
#' Obtains point estimates and standard errors for the group average treatment effects (GATEs), where groups correspond to the
#' leaves of an \code{\link[rpart]{rpart}} object. Additionally, performs some hypothesis testing.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param D Treatment assignment vector
#' @param X Covariate matrix (no intercept).
#' @param method Either \code{"raw"} or \code{"aipw"}, defines the outcome used in the regression.
#' @param scores Optional, vector of scores to be used in the regression. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#' @param boot_ci Logical, whether to compute bootstrap confidence intervals.
#' @param boot_R Number of bootstrap replications. Ignored if \code{boot_ci == FALSE}.
#'
#' @return
#' A list storing:
#'   \item{\code{model}}{The model fitted to get point estimates and standard errors for the GATEs, as an \code{\link[estimatr]{lm_robust}} object.}
#'   \item{\code{gates_diff_pairs}}{Results of testing whether GATEs differ across all pairs of leaves. This is a list storing GATEs differences and p-values adjusted using Holm's procedure (check \code{\link[stats]{p.adjust}}). \code{NULL} if the tree consists of a root only.}
#'   \item{\code{boot_ci}}{Bootstrap confidence intervals (this is an empty list if \code{boot_ci == FALSE}.}
#'   \item{\code{scores}}{Vector of doubly robust scores. \code{NULL} if \code{method == 'raw'}.}
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Split the sample.
#' splits <- sample_split(length(y), training_frac = 0.5)
#' training_idx <- splits$training_idx
#' honest_idx <- splits$honest_idx
#'
#' y_tr <- y[training_idx]
#' D_tr <- D[training_idx]
#' X_tr <- X[training_idx, ]
#'
#' y_hon <- y[honest_idx]
#' D_hon <- D[honest_idx]
#' X_hon <- X[honest_idx, ]
#'
#' ## Construct a tree using training sample.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame("y" = y_tr, X_tr), maxdepth = 2)
#'
#' ## Estimate GATEs in each node (internal and terminal) using honest sample.
#' results <- causal_ols_rpart(tree, y_hon, D_hon, X_hon, method = "raw")
#'
#' summary(results$model) # Coefficient of leafk:D is GATE in k-th leaf.
#'
#' results$gates_diff_pair$gates_diff # GATEs differences.
#' results$gates_diff_pair$holm_pvalues # leaves 1-2 and 3-4 not statistically different.
#'
#' @md
#' @details
#' ## Point estimates and standard errors for the GATEs
#' The GATEs and their standard errors are obtained by fitting an appropriate linear model. If \code{method == "raw"}, we
#' estimate via OLS the following:
#'
#' \deqn{Y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th leaf of \code{tree}, and |T| the number of
#' groups. If the treatment is randomly assigned, one can show that the betas identify the GATE in each leaf. However, this is not true
#' in observational studies due to selection into treatment. In this case, the user is expected to use \code{method == "aipw"} to run
#' the following regression:
#'
#' \deqn{score_i = \sum_{l = 1}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' where score_i are doubly-robust scores constructed via honest regression forests and 5-fold cross fitting (unless the user specifies
#' the argument \code{scores}). This way, betas again identify the GATEs.\cr
#'
#' Regardless of \code{method}, standard errors are estimated via the Eicker-Huber-White estimator.\cr
#'
#' If \code{boot_ci == TRUE}, the routine also computes asymmetric bias-corrected and accelerated 95% confidence intervals using 2000 bootstrap
#' samples.\cr
#'
#' If \code{tree} consists of a root only, \code{causal_ols_rpart} regresses \code{y} on a constant and \code{D} if
#' \code{method == "raw"}, or regresses the doubly-robust scores on a constant if \code{method == "aipw"}. This way,
#' we get an estimate of the overall average treatment effect.
#'
#' ## Hypothesis testing
#' \code{\link{causal_ols_rpart}} uses the standard errors obtained by fitting the linear models above to test the hypotheses
#' that the GATEs are different across all pairs of leaves. Here, we adjust p-values to account for multiple hypotheses testing
#' using Holm's procedure.
#'
#' ## Caution on Inference
#' "honesty" is a necessary requirement to get valid inference. Thus, observations in \code{y}, \code{D}, and
#' \code{X} must not have been used to construct the \code{tree} and the \code{scores}.\cr
#'
#' @import rpart estimatr car stats
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso \code{\link{estimate_rpart}} \code{\link{avg_characteristics_rpart}}
#'
#' @export
causal_ols_rpart <- function(tree, y, D, X, method = "aipw", scores = NULL,
                             boot_ci = FALSE, boot_R = 2000) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("Invalid 'tree'. This must be an rpart object.", call. = FALSE)
  if(!(method %in% c("raw", "aipw"))) stop("Invalid 'method'. This must be either 'raw' or 'aipw'.", call. = FALSE)
  if (!(boot_ci %in% c(FALSE, TRUE))) stop("Invalid 'boot_ci'. This must be either FALSE or TRUE.", call. = FALSE)
  if (boot_R < 0) stop("Invalid 'boot_R'. This must be a positive integer.", call. = FALSE)

  leaves <- leaf_membership(tree, X)
  n_leaves <- get_leaves(tree)

  if (length(unique(leaves)) < n_leaves) warning("One or more leaves are empty: No observations in 'X' falls there.")

  for (leaf in seq_len(n_leaves)) {
    if (sum(leaves == leaf) == 1) stop("One or more leaves contain only one observation of 'X', thus we cannot get standard errors there. \nTry with a lower number of groups or a bigger 'X'.", call. = FALSE)
    if (method == "raw" & length(unique(D[leaves == leaf])) == 1) stop("One or more leaves contain only treated or control observations. This is incompatible with 'method = raw'.", call. = FALSE)
  }

  ## GATEs point estimates and standard errors.
  if (method == "raw") {
    if (length(unique(leaves)) == 1) {
      model <- estimatr::lm_robust(y ~ D, data = data.frame("y" = y, "D" = D), se_type = "HC1")
    } else {
      model <- estimatr::lm_robust(y ~ 0 + leaf + D:leaf, data = data.frame("y" = y, "leaf" = leaves, "D" = D), se_type = "HC1")
    }
  } else if (method == "aipw") {
    if (is.null(scores)) {
      scores <- dr_scores(y, D, X)
    }

    if (length(unique(leaves)) == 1) {
      model <- estimatr::lm_robust(scores ~ 1, data = data.frame("scores" = scores), se_type = "HC1")
    } else {
      model <- estimatr::lm_robust(scores ~ 0 + leaf, data = data.frame("scores" = scores, "leaf" = leaves), se_type = "HC1")
    }
  }

  ## Test if GATEs are different across all pairs of leaves. Adjust p-values by Holm's procedure.
  if (n_leaves > 1) { # Inspired by https://gsbdbi.github.io/ml_tutorial/hte_tutorial/hte_tutorial.html#hte_1:_causal_trees.
    gates_point <- model$coefficients

    differences <- matrix(NA, nrow = n_leaves, ncol = n_leaves)
    rownames(differences) <- paste0("leaf", seq_len(n_leaves))
    colnames(differences) <- paste0("leaf", seq_len(n_leaves))

    p_values <- matrix(NA, nrow = n_leaves, ncol = n_leaves)
    rownames(p_values) <- paste0("leaf", seq_len(n_leaves))
    colnames(p_values) <- paste0("leaf", seq_len(n_leaves))

    if (method == "raw") {
      for (i in seq_len(n_leaves-1)) {
        for (j in seq_len(n_leaves)[-c(1:i)]) {
         leaf_i <- paste0("leaf", i, ":D")
         leaf_j <- paste0("leaf", j, ":D")

         differences[j, i] <- gates_point[leaf_j] - gates_point[leaf_i]
         p_values[j, i] <- car::linearHypothesis(model, paste(leaf_j, "=", leaf_i), test = "F")$`Pr(>F)`[2]
       }
     }
    } else if (method == "aipw") {
      for (i in seq_len(n_leaves-1)) {
        for (j in seq_len(n_leaves)[-c(1:i)]) {
          leaf_i <- paste0("leaf", i)
          leaf_j <- paste0("leaf", j)

          differences[j, i] <- gates_point[leaf_j] - gates_point[leaf_i]
          p_values[j, i] <- car::linearHypothesis(model, paste(leaf_j, "=", leaf_i), test = "F")$`Pr(>F)`[2]
        }
      }
    }

    p_values_vec <- c(p_values) # First column, then second columns, then third columns ...
    p_values_holm_vec <- stats::p.adjust(p_values_vec, method = "holm")
    p_values_holm <- matrix(p_values_holm_vec, nrow = n_leaves, ncol = n_leaves)

    gates_diff_pairs <- list("gates_diff" = differences, "holm_pvalues" = p_values_holm)
  } else {
    gates_diff_pairs <- NULL
  }

  ## Bootstrap confidence intervals, if required.
  if (boot_ci) {
    # Extract GATEs indexes according to estimation strategy.
    if (method == "raw") {
      if (length(unique(leaves)) == 1) {
        gates_idx <- 2
      } else {
        gates_idx <- which(sapply(names(model$coefficients), function(x) grepl(":D", x)))
      }
    } else if (method == "aipw") {
      if (length(unique(leaves)) == 1) {
        gates_idx <- 1
      } else {
        gates_idx <- which(sapply(names(model$coefficients), function(x) grepl("leaf", x)))
      }
    }

    # Define function input for boot::boot().
    boot_fun <- function(data, idx) {
      data_star <- data[idx, ] # It must contain y, D, scores, and leaves.
      if (method == "raw") {
        if (length(unique(leaves)) == 1) {
          model <- estimatr::lm_robust(y ~ D, data = data.frame("y" = data_star$y, "D" = data_star$D), se_type = "HC1")
        } else {
          model <- estimatr::lm_robust(y ~ 0 + leaf + D:leaf, data = data.frame("y" = data_star$y, "leaf" = data_star$leaves, "D" = data_star$D), se_type = "HC1")
        }
      } else if (method == "aipw") {
        if (length(unique(leaves)) == 1) {
          model <- estimatr::lm_robust(scores ~ 1, data = data.frame("scores" = scores[idx]), se_type = "HC1")
        } else {
          model <- estimatr::lm_robust(scores ~ 0 + leaf, data = data.frame("scores" = scores[idx], "leaf" = data_star$leaves), se_type = "HC1")
        }
      }

      return(coef(model)[gates_idx])
    }

    # Run bootstrap and compute confidence intervals.
    boot_out <- boot::boot(data.frame(y, D, X, leaves), boot_fun, R = boot_R)

    boot_ci_lower <- broom::tidy(boot_out, conf.int = TRUE, conf.method = "bca")$conf.low
    boot_ci_upper <- broom::tidy(boot_out, conf.int = TRUE, conf.method = "bca")$conf.high

    names(boot_ci_lower) <- names(gates_idx)
    names(boot_ci_upper) <- names(gates_idx)
  }

  ## Output.
  return(list("model" = model,
              "gates_diff_pairs" = gates_diff_pairs,
              "boot_ci" = if (boot_ci) list("lower" = boot_ci_lower, "upper" = boot_ci_upper) else list(),
              "scores" = scores))
}


#' Leaves Average Characteristics
#'
#' Computes the average characteristics of units in each leaf of an \code{\link[rpart]{rpart}} object.
#'
#' @param tree An \code{rpart} object.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' A list storing each regression as an \code{\link[estimatr]{lm_robust}} object.
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 3
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#'
#' ## Construct a tree.
#' library(rpart)
#' tree <- rpart(y ~ ., data = data.frame("y" = y, X), maxdepth = 2)
#'
#' ## Compute average characteristics in each leaf.
#' results <- avg_characteristics_rpart(tree, X)
#' results
#'
#' @details
#' \code{\link{avg_characteristics_rpart}} regresses each covariate on a set of dummies denoting leaf membership.
#' This way, we get the average characteristics of units in each leaf, together with a standard error.\cr
#'
#' Leaves are ordered in increasing order of their predictions (from most negative to most positive).\cr
#'
#' Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' @import estimatr stats
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso \code{\link{causal_ols_rpart}}, \code{\link{estimate_rpart}}
#'
#' @export
avg_characteristics_rpart <- function(tree, X) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("Invalid 'tree'. This must be an rpart object.", call. = FALSE)

  ## Generate leaves indicators.
  leaves <- leaf_membership(tree, X)

  ## Regress each leaf on the leaf indicator.
  regressions <- apply(X, MARGIN = 2, function(x) {estimatr::lm_robust(x ~ 0 + leaf, data = data.frame("x" = x, "leaf" = leaves), se_type = "HC1")})

  ## Output.
  return(regressions)
}
