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
#' @import rpart
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{get_leaves}} \code{\link{node_membership}} \code{\link{leaf_membership}}
#'
#' @export
subtree <- function(tree, leaves = NULL, cv = FALSE) {
  ## Handling inputs and checks.
  if(!(inherits(tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(cv %in% c(TRUE, FALSE))) stop("'cv' must be either TRUE or FALSE.", call. = FALSE)
  if (is.null(leaves) & cv == FALSE) stop("Invalid combination of 'leaves' and 'cv'. Please specify a number of leaves or select the cross-validation option.", call. = FALSE)
  if (!is.null(leaves)) {
    if (leaves < 1) stop("'leaves' must be a positive number.", call. = FALSE)
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
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{subtree}} \code{\link{node_membership}} \code{\link{leaf_membership}}
#'
#' @export
get_leaves <- function(tree) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.")

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
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{subtree}} \code{\link{node_membership}} \code{\link{get_leaves}}
#'
#' @export
leaf_membership <- function(tree, X) {
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.")
  if (!is.matrix(X) & !is.data.frame(X)) stop("'X' must be either a matrix or a data frame.", call. = FALSE)

  ## Inspired by https://bookdown.org/halflearned/tutorial/hte1.html.
  tree_predictions <- predict(tree, data.frame(X))
  n_leaves <- length(unique(tree_predictions))
  leaves <- factor(tree_predictions, levels = sort(unique(tree_predictions)), labels = seq(n_leaves))

  return(leaves)
}


#' Estimation of rpart Objects
#'
#' Replace node predictions of an \code{\link[rpart]{rpart}} object using external data.
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
#' @details
#' If \code{method == "raw"}, \code{estimate_rpart} replaces node predictions with the differences between the sample average
#' of the observed outcomes of treated units and the sample average of the observed outcomes of control units in each node,
#' which is an unbiased estimator of GATEs if the assignment to treatment is randomized.\cr
#'
#' If \code{method == "aipw"}, \code{estimate_rpart} replaces node predictions with sample averages of doubly-robust
#' scores in each node. This is a valid estimator of GATEs in observational studies. Honest regression forests
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
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)
  if (!(method %in% c("raw", "aipw"))) stop("You must provide a valid method.", call. = FALSE)

  new_tree <- tree

  ## Replace node predictions.
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


#' Estimating Leaf-Effects via Linear Models
#'
#' Uses the leaves of a tree stored in an \code{\link[rpart]{rpart}} object to estimate a linear model via OLS. The
#' estimated coefficients identify the average treatment effect in each leaf under unconfoundedness. If the data used
#' in the OLS estimation have not been used to grow the tree (a condition called "honesty"), then one can use the
#' standard errors for the tree's estimates to construct valid confidence intervals.
#'
#' @param tree An \code{\link[rpart]{rpart}} object.
#' @param y Outcome vector.
#' @param D Treatment assignment vector
#' @param X Covariate matrix (no intercept).
#' @param method Either \code{"raw"} or \code{"aipw"}, defines the outcome used in the regression.
#' @param scores Optional, vector of scores to be used in the regression. Useful to save computational time if scores have already been estimated. Ignored if \code{method == "raw"}.
#'
#' @return
#' A list storing:
#'   \item{\code{model}}{The fitted model, as an \code{\link[estimatr]{lm_robust}} object.}
#'   \item{\code{gates_all_equal}}{Results of testing simultaneously whether all average effects are the same, as an \code{anova} object (check \code{\link[car]{linearHypothesis}}).}
#'   \item{\code{gates_diff_smallest}}{Results of testing individually whether each average effect differs from smallest effect. p-values are adjusted using Holm's procedure (check \code{\link[stats]{p.adjust}}). \code{NULL} if the tree consists of a root only.}
#'   \item{\code{scores}}{Vector of doubly robust scores. \code{NULL} if \code{method == 'raw'}.}
#'
#' @details
#' The \code{method} argument controls how average effects in each leaf are estimated. If \code{"method" == "raw"}, we estimate via
#' OLS the following linear model:
#'
#' \deqn{Y_i = \sum_{l = 1}^{|T|} L_{i, l} \gamma_l + \sum_{l = 1}^{|T|} L_{i, l} D_i \beta_l + \epsilon_i}
#'
#' with \code{L_{i, l}} a dummy variable equal to one if the i-th unit falls in the l-th leaf of the tree, and \code{|T|} the
#' number of leaves. If the treatment is randomly assigned, one can show that the betas identify the average treatment effect in
#' each group. However, in observational studies these estimates are biased due to selection into treatment. To get unbiased
#' estimates, we can set \code{"method"} to \code{"aipw"} to construct doubly-robust scores \code{y_i^*} and use them as a
#' pseudo-outcome in the following regression:
#'
#' \deqn{Y_i^* = \sum_{l = 1}^{|T|} L_{i, l} \beta_l + \epsilon_i}
#'
#' Honest regression forests and 5-fold cross fitting are used to estimate the propensity
#' score and the conditional mean function of the outcome (unless the user specifies the argument \code{scores}).\cr
#'
#' Additionally, \code{causal_ols_rpart} performs two types of hypothesis testing. First, it uses standard
#' errors from the above models to test whether average treatment effects in each leaf are the same by constructing a finite-sample
#' F statistic for carrying out a Wald-test-based comparison between a model and a linearly restricted model. Second, it fits a new
#' linear model by omitting the leaf with the smallest (i.e., more negative) average effect (and its interaction with the treatment
#' variable if \code{method == "raw"}) and adding an intercept (and the treatment variable if \code{method == "raw"}). This way,
#' coefficients give how much each average effect is larger than the smallest effect, and we can test separately the usual null
#' hypotheses that these differences are zero. We correct for multiple hypothesis testing using the Holm's procedure.\cr
#'
#' Notice that "honesty" is a necessary requirement to get valid inference. Thus, observations in \code{y}, \code{D}, and
#' \code{X} must not have been used to construct the \code{tree}.\cr
#'
#' Regardless of \code{method}, standard errors are estimated via the Eicker-Huber-White estimator.\cr
#'
#' If \code{tree} consists of a root only, \code{causal_ols_rpart} regresses \code{y} on a constant and \code{D} if
#' \code{method == "raw"}, or regresses the doubly-robust scores on a constant if \code{method == "aipw"}. This way,
#' we get an estimate of the overall average treatment effect.
#'
#' @examples
#'
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
causal_ols_rpart <- function(tree, y, X, D, method = "aipw", scores = NULL) {
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
    if (is.null(scores)) {
      scores <- dr_scores(y, D, X)
    }

    if (length(unique(leaves)) == 1) {
      model <- estimatr::lm_robust(scores ~ 1, data = data.frame("scores" = scores), se_type = "HC1")
    } else {
      model <- estimatr::lm_robust(scores ~ 0 + leaf, data = data.frame("scores" = scores, "leaf" = leaves), se_type = "HC1")
    }
  }

  ## Test whether GATEs are the same across leaves.
  if (method == "raw") {
    null <- paste0("leaf1:D = leaf", seq(2, get_leaves(tree)), ":D")
    gates_all_equal <- car::linearHypothesis(model, null, test = "F")
  } else if (method == "aipw") {
    null <- paste0("leaf1 = leaf", seq(2, get_leaves(tree)))
    gates_all_equal <- car::linearHypothesis(model, null, test = "F")
  }

  ## Test whether each GATE is different from smallest (more negative) GATE. Fit a new model so that coefficients give difference
  ## of each GATE from smallest GATE. Use Holm's correction to test if differences are zero.
  if (get_leaves(tree) > 1) {
    if (method == "raw") {
      new_model <- estimatr::lm_robust(y ~ leaf*D, data = data.frame("y" = y, "leaf" = leaves, "D" = D), se_type = "HC1")
      parms_idx <- which(sapply(names(coef(new_model)), function(x) grepl(":D", x)))
    } else if (method == "aipw") {
      if (is.null(scores)) scores <- dr_scores(y, D, X)
      new_model <- estimatr::lm_robust(scores ~ leaf, data = data.frame("scores" = scores, "leaf" = leaves), se_type = "HC1")
      parms_idx <- which(sapply(names(coef(new_model)), function(x) grepl("leaf", x)))
    }

    p_values <- new_model$p.value
    p_values_holm <- stats::p.adjust(p_values, method = "holm")

    gates_diff_smallest <- data.frame("GATE_increase" = new_model$coefficients[parms_idx],
                         "se" = new_model$std.error[parms_idx],
                         "adj_pvalue" = round(p_values_holm[parms_idx], 3))
    rownames(gates_diff_smallest) <- paste0("leaf", seq_len(get_leaves(tree))[-1])
  } else {
    gates_diff_smallest <- NULL
  }

  ## Output.
  return(list("model" = model, "gates_all_equal" = gates_all_equal, "gates_diff_smallest" = gates_diff_smallest, "scores" = scores))
}


#' Leaves Average Characteristics
#'
#' Compute the average characteristics of units in each leaf of an \code{\link[rpart]{rpart}} object.
#'
#' @param tree An \code{rpart} object.
#' @param X Covariate matrix (no intercept).
#' @param gates_point Estimated GATEs, one for each leaf, in increasing order (from most negative to most positive).
#' @param gates_sd Standard errors for estimated GATEs, one for each leaf.
#'
#' @return
#' Prints LATEX code in the console.
#'
#' @details
#' \code{avg_characteristics_rpart} regresses each covariate on a set of dummies denoting leaf membership.
#' This way, we get the average characteristics of units in each leaf, together with a standard error. Leaves are
#' ordered in increasing order of their predictions (from most negative to most positive). Standard errors are
#' estimated via the Eicker-Huber-White estimator.\cr
#'
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox},
#' \code{multirow}.
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
avg_characteristics_rpart <- function(tree, X, gates_point = NULL, gates_sd = NULL) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  ## Generate leaves indicators.
  leaves <- leaf_membership(tree, X)

  ## Regress each leaf on the leaf indicator.
  regressions <- apply(X, MARGIN = 2, function(x) {estimatr::lm_robust(x ~ 0 + leaf, data = data.frame("x" = x, "leaf" = leaves), se_type = "HC1")})

  ## Extract information.
  parms <- lapply(regressions, function(x) {stats::coef(summary(x))[, c("Estimate", "Std. Error")]})
  if (is.null(gates_point)) gates_point <- rep("NA", lenght = length(unique(leaves)))
  if (is.null(gates_sd)) gates_sd <- rep("NA", lenght = length(unique(leaves)))

  gates_point <- round(gates_point, 3)
  gates_sd <- round(gates_sd, 3)
  gates_ci_lower <- round(gates_point - 1.96 * gates_sd, 3)
  gates_ci_upper <- round(gates_point + 1.96 * gates_sd, 3)

  ## Write table.
  table_names <- rename_latex(colnames(X))

  cat("\\begingroup
  \\setlength{\\tabcolsep}{8pt}
  \\renewcommand{\\arraystretch}{1.1}
  \\begin{table}[b!]
    \\centering
    \\begin{adjustbox}{width = 1\\textwidth}
    \\begin{tabular}{@{\\extracolsep{5pt}}l ", rep("c ", length(unique(leaves))), "}
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      & ", c(paste("\\textit{Leaf ", 1:(length(unique(leaves))-1), "} & ", sep = ""), paste("\\textit{Leaf ", length(unique(leaves)), "}", sep = "")) ," \\\\
      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n", sep = "")

  cat("      \\multirow{2}{*}{GATEs} & ", paste(gates_point[1:(length(unique(leaves))-1)], " & ", sep = ""), gates_point[length(unique(leaves))], " \\\\
      & ", paste("[", gates_ci_lower[1:(length(unique(leaves))-1)], ", ", gates_ci_upper[1:(length(unique(leaves))-1)], "] & ", sep = ""), paste("[", gates_ci_lower[length(unique(leaves))], ", ", gates_ci_upper[length(unique(leaves))], "]", sep = ""), " \\\\ \n\n", sep = "")
  cat("      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n")

  for (i in seq_len(length(table_names))) {
    cat("      \\texttt{", table_names[i], "} & ", paste(round(parms[[i]][1:(length(unique(leaves))-1), 1], 2), " & ", sep = ""), round(parms[[i]][length(unique(leaves)), 1], 2), " \\\\ \n",
        "      & ", paste("(", round(parms[[i]][1:(length(unique(leaves))-1), 2], 3), ")", " & ", sep = ""), paste("(", round(parms[[i]][length(unique(leaves)), 2], 3), ")", sep = ""), " \\\\ \n", sep = "")
  }

  cat("\n      \\addlinespace[3pt]
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
    \\end{tabular}
    \\end{adjustbox}
    \\caption{Average characteristics of units in each leaf, obtained by regressing each covariate on a set of dummies denoting leaf membership. Standard errors are estimated via the Eicker-Huber-White estimator and reported in parenthesis under each point estimate. For each leaf, point estimates and $95\\%$ confidence intervals for GATEs are displayed.}
    \\label{table:average.characteristics.leaves}
    \\end{table}
\\endgroup \n\n")

  no_variation_names <- names(unlist(lapply(parms, function(x) { if (any(x[, 2] == 0)) idx <- 1})))

  ## Warn the user for zero variation in leaves.
  if (length(no_variation_names) > 1) {
    names_string <- paste0(no_variation_names, collapse = ", ")
    warning(paste0("Variables '", names_string, "' have no variation in one or more leaves. Please correct the table by removing the associated standard errors."))
  } else if (length(no_variation_names) == 1) {
    warning(paste0("Variable '", no_variation_names, "' has no variation in one or more leaves. Please correct the table by removing the associated standard errors."))

  }
}
