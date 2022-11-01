#' Leaves Average Characteristics
#'
#' Compute the average characteristics of units in each leaf of an \code{aggTrees} object.
#'
#' @param object An \code{aggTrees} object.
#' @param cates Estimated CATEs.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' Prints LATEX code in the console.
#'
#' @details
#' The fact that a tree does not split on a particular covariate does not imply that such covariate is not related
#' to treatment effect heterogeneity. In fact, there are many ways in which groups with low or high treatment effect can be formed.
#' Looking at how the average characteristics of the units vary across the leaves of the tree is a more systematic way to assess
#' how treatment effects relate to the covariates.\cr
#'
#' To achieve this, \code{avg_characteristics_aggtree} regresses each covariate on a categorical variable denoting membership to
#' a particular leaf. This way, we get the average characteristics of units in each leaf, together with a standard error.\cr
#'
#' Standard errors are estimated via the Eicker-Huber-White estimator.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{honest_ols}}
#'
#' @export
avg_characteristics_aggtree <- function(object, cates, X) {
  ## Handling inputs and checks.
  if (!(inherits(object, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(object$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(object$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- object$tree

  ## Compute average characteristics.
  avg_char <- avg_characteristics_rpart(tree, cates, X)
}


#' Leaves Average Characteristics
#'
#' Compute the average characteristics of units in each leaf of a \code{\link[rpart]{rpart}} object.
#'
#' @param tree A \code{rpart} object.
#' @param y Outcome vector.
#' @param X Covariate matrix (no intercept).
#'
#' @return
#' Prints LATEX code in the console.
#'
#' @details
#' \code{avg_characteristics_aggtree} regresses each covariate on a categorical variable denoting membership to
#' a particular leaf. This way, we get the average characteristics of units in each leaf, together with a standard error.\cr
#'
#' Standard errors are estimated via the Eicker-Huber-White estimator.\cr
#'
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#'
#' @import estimatr
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{honest_ols}}, \code{\link{aggregation_tree}}
#'
#' @export
avg_characteristics_rpart <- function(tree, y, X) {
  ## Handling inputs and checks.
  if (!inherits(tree, "rpart")) stop("'tree' must be a rpart object.", call. = FALSE)

  ## Generate leaves indicators. Inspired by https://bookdown.org/halflearned/tutorial/hte1.html.
  tree_predictions <- predict(tree, data.frame(X))
  n_leaves <- length(unique(tree_predictions))
  leaves <- factor(tree_predictions, levels = sort(unique(tree_predictions)), labels = seq(n_leaves))

  ## Regress each leaf on the leaf indicator.
  regressions <- apply(X, MARGIN = 2, function(x) {estimatr::lm_robust(x ~ 0 + leaf, data = data.frame("x" = x, "leaf" = leaves), , se_type = "HC1")})

  ## Extract information.
  parms <- lapply(regressions, function(x) {coef(summary(x))[, c("Estimate", "Std. Error")]})
  gates <- round(sapply(sort(unique(leaves)), function(x) {mean(y[leaves == x])}), 2)

  ## Write table.
  table_names <- rename_latex(colnames(X))

  cat("      \\begingroup
        \\setlength{\\tabcolsep}{8pt}
        \\renewcommand{\\arraystretch}{1.1}
        \\begin{table}[H]
          \\centering
          \\begin{adjustbox}{width = 0.75\\textwidth}
          \\begin{tabular}{@{\\extracolsep{5pt}}l", rep("c", length(unique(leaves))), "}
          \\\\[-1.8ex]\\hline
          \\hline \\\\[-1.8ex]
          & ", c(paste("Leaf", 1:(length(unique(leaves))-1), " & ", sep = ""), paste("Leaf", length(unique(leaves)), sep = "")) ," \\\\
          \\addlinespace[2pt]
          \\hline \\\\[-1.8ex] \n\n", sep = "")

  cat("          GATEs & ", paste(gates[1:(length(unique(leaves))-1)], " & ", sep = ""), gates[length(unique(leaves))], "\\\\ \n", sep = "")


  for (i in seq_len(length(table_names))) {
  cat("          \\texttt{", table_names[i], "} & ", paste(round(parms[[i]][1:(length(unique(leaves))-1), 1], 2), " & ", sep = ""), round(parms[[i]][length(unique(leaves)), 1], 2), " \\\\ \n",
      "                        & ", paste("(", round(parms[[i]][1:(length(unique(leaves))-1), 2], 3), ")", " & ", sep = ""), paste("(", round(parms[[i]][length(unique(leaves)), 2], 3), ")", sep = ""), " \\\\ \n", sep = "")
  }

  cat("\n          \\addlinespace[3pt]
          \\\\[-1.8ex]\\hline
          \\hline \\\\[-1.8ex]
          \\end{tabular}
          \\end{adjustbox}
          \\caption{https://soundcloud.com/theplasticchairband}
          \\label{table:average.characterstics.leaves}
        \\end{table}
      \\endgroup")
}
