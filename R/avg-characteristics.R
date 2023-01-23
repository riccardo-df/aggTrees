#' Leaves Average Characteristics
#'
#' Compute the average characteristics of units in each leaf of an \code{\link[rpart]{rpart}} object.
#'
#' @param tree A \code{rpart} object.
#' @param X Covariate matrix (no intercept).
#' @param gates_point Estimated GATEs, one for each leaf.
#' @param gates_point Standard errors for estimated GATEs, one for each leaf.
#'
#' @return
#' Prints LATEX code in the console.
#'
#' @details
#' \code{avg_characteristics_rpart} regresses each covariate on a set of dummies denoting leaf membership.
#' This way, we get the average characteristics of units in each leaf, together with a standard error. Standard errors are
#' estimated via the Eicker-Huber-White estimator.\cr
#'
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#'
#' @import estimatr stats
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item V Chernozhukov, M Demirer, E Duflo, I Fernandez-Val (2018). Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments, with an application to immunization in India. arXiv preprint arXiv:1712.04802. \doi{10.48550/ARXIV.1712.04802}.
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso \code{\link{aggregation_tree}}, \code{\link{causal_ols_rpart}}, \code{\link{estimate_rpart}}
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

  ## Write table.
  table_names <- rename_latex(colnames(X))

  cat("\\begingroup
        \\setlength{\\tabcolsep}{8pt}
        \\renewcommand{\\arraystretch}{1.1}
        \\begin{table}[H]
          \\centering
          \\begin{adjustbox}{width = 0.75\\textwidth}
          \\begin{tabular}{@{\\extracolsep{5pt}}l ", rep("c ", length(unique(leaves))), "}
          \\\\[-1.8ex]\\hline
          \\hline \\\\[-1.8ex]
          & ", c(paste("\\textit{Leaf ", 1:(length(unique(leaves))-1), "} & ", sep = ""), paste("\\textit{Leaf ", length(unique(leaves)), "}", sep = "")) ," \\\\
          \\addlinespace[2pt]
          \\hline \\\\[-1.8ex] \n\n", sep = "")

  cat("          GATEs & ", paste(gates_point[1:(length(unique(leaves))-1)], " & ", sep = ""), gates_point[length(unique(leaves))], " \\\\
                & ", paste("(", gates_sd[1:(length(unique(leaves))-1)], ") & ", sep = ""), paste("(", gates_sd[length(unique(leaves))], ")", sep = ""), " \\\\ \n\n", sep = "")
  cat("          \\addlinespace[2pt]
          \\hline \\\\[-1.8ex] \n\n")

  for (i in seq_len(length(table_names))) {
  cat("          \\texttt{", table_names[i], "} & ", paste(round(parms[[i]][1:(length(unique(leaves))-1), 1], 2), " & ", sep = ""), round(parms[[i]][length(unique(leaves)), 1], 2), " \\\\ \n",
      "                        & ", paste("(", round(parms[[i]][1:(length(unique(leaves))-1), 2], 3), ")", " & ", sep = ""), paste("(", round(parms[[i]][length(unique(leaves)), 2], 3), ")", sep = ""), " \\\\ \n", sep = "")
  }

  cat("\n          \\addlinespace[3pt]
          \\\\[-1.8ex]\\hline
          \\hline \\\\[-1.8ex]
          \\end{tabular}
          \\end{adjustbox}
          \\caption{\\href{https://soundcloud.com/theplasticchairband}{Secret link.}}
          \\label{table:average.characterstics.leaves}
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
