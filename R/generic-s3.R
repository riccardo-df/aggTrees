#' Plot Method for aggTrees Objects
#'
#' Plots an \code{aggTrees} object.
#'
#' @param x An \code{aggTrees} object.
#' @param leaves Number of leaves of the desired tree. This can be used to plot subtrees.
#' @param type Plotting style.
#' @param palette Palette to color the nodes.
#' @param sequence If \code{TRUE}, the whole sequence of optimal partitions is displayed in a short animation.
#' @param ... Further arguments from \code{\link[rpart.plot]{prp}}.
#'
#' @details
#' Due to coding limitations, honest trees from \code{aggTrees} objects show honest estimates only in their leaves. Thus, if
#' the user wants to plot an honest tree, \code{type} 3 or 5 is recommended to avoid misleading plots.\cr
#'
#' \code{palette} can be either a vector of colors, or a function that takes as an argument the number of nodes of a tree and
#' returns a vector of colors. Please refer to \code{\link[colorspace]{choose_palette}} for high-quality palettes.\cr
#'
#' @import rpart.plot
#'
#' @author Riccardo Di Francesco
#'
#' @seealso
#' \code{\link{aggregation_tree}}, \code{\link{subtree_aggtree}}, \code{\link{recursive_partitioning_plot}}
#'
#' @export
plot.aggTrees <- function(x, leaves = get_leaves(x$tree),
                          type = 2, palette = c("BuRd"),
                          sequence = FALSE, ...) {
  ## Handling inputs and checks.
  if (!(inherits(x, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(inherits(x$tree, "rpart"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(x$honesty %in% c(TRUE, FALSE))) stop("You must provide a valid aggTrees object.", call. = FALSE)

  tree <- x$tree

  if (x$honesty == TRUE & !(type %in% c(3, 5))) warning("Only the leaf estimates are honest. Internal nodes show non-honest estimates. Consider using 'type' 3 or 5.")

  ## Plotting.
  if (sequence) {
    prefix <- c("ATE = ", rep("GATE = ", times = (length(tree$frame$n) - 1)))
    sizes <- tree$frame$n
    percentages <- round(tree$frame$n / tree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    full_nodes <- rownames(tree$frame) # Storing rownames of full tree (they serve as node labels).
    alpha_values <- rev(tree$cptable[, "CP"]) # Threshold values of cost-complexity parameter.

    for (alpha in alpha_values) {
      temp_nodes <- rownames(rpart::prune.rpart(tree, alpha)$frame) # Node labels of pruned tree.
      colors <- ifelse(full_nodes %in% temp_nodes, 1, "white") # Graying out collapsed nodes.

      grDevices::dev.hold()

      rpart.plot::prp(tree,
                      type = type,
                      extra = 0,
                      under = FALSE,
                      fallen.leaves = TRUE,
                      round = 0,
                      leaf.round = 0,
                      prefix = prefix,
                      suffix = suffix,
                      pal.thresh = tree$frame$yval[2],
                      box.palette = if (is.function(palette)) palette(nrow(tree$frame)) else palette,
                      branch = 0.3,
                      branch.col = colors,
                      split.col = colors,
                      col = colors,
                      ...)

      grDevices::dev.flush()
      Sys.sleep(1)
    }
  } else {
    subtree <- subtree_rpart(tree, leaves)
    prefix <- c("ATE = ", rep("GATE = ", times = (length(subtree$frame$n) - 1)))
    sizes <- subtree$frame$n
    percentages <- round(subtree$frame$n / subtree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    rpart.plot::prp(subtree,
                    type = type,
                    extra = 0,
                    under = FALSE,
                    fallen.leaves = TRUE,
                    round = 0,
                    leaf.round = 0,
                    prefix = prefix,
                    suffix = suffix,
                    pal.thresh = subtree$frame$yval[2],
                    box.palette = if (!is.function(palette)) palette else palette(nrow(subtree$frame)),
                    branch = 0.3,
                    ...)
  }
}


#' Summary Method for aggTrees Objects
#'
#' Summarizes an \code{aggTrees} object.
#'
#' @param object \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @author Riccardo Di Francesco
#'
#' @export
summary.aggTrees <- function(object, ...) {
  cat("Honesty: ",object$honesty, "\n", sep = "")
  summary(object$tree)
}


#' Print Method for aggTrees Objects
#'
#' Prints an \code{aggTrees} object.
#'
#' @param x \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @seealso \code{\link{aggregation_tree}}
#'
#' @author Riccardo Di Francesco
#'
#' @export
print.aggTrees <- function(x, ...) {
  cat("Honesty: ",x$honesty, "\n", sep = "")
  print(x$tree)
}


#' Print Method for \code{phased} objects
#'
#' Prints a \code{phased} object.
#'
#' @param x A \code{phased} object.
#' @param latex If TRUE, prints LATEX code for the table.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#'
#' @export
print.phased <- function(x, latex = FALSE, ...) {
  ## Handling input and check.
  if (!inherits(x, "phased")) stop("You need to pass a 'phased' object.")
  if (!(latex %in% c(TRUE, FALSE))) stop("latex must be either TRUE or FALSE.")

  measures <- data.frame(cbind(t(x$descriptive_stats), t(x$norm_diff), t(x$log_ratio_sd)))
  col_names <- c("Mean", "S.D.", "Mean", "S.D.", "Norm.Diff.", "Log S.D.")
  row_names <- rownames(measures)

  ## Arranging for printing.
  negative <- measures < 0 # This is useful later to handle negative values.

  # Converting to char, and automating desired widths.
  measures <- lapply(measures, sprintf, fmt = "%.3f")

  var_width <- max(sapply(x$var_names, nchar))
  metrics_width <- sapply(measures, function(x) max(nchar(x)))
  col_names_width <- sapply(col_names, nchar)
  max_col_width <- apply(rbind(metrics_width, col_names_width), 2, max)

  # Formatting according to desired widths.
  measures <- lapply(1:6, function(x) format(measures[[x]], width = max_col_width[x], justify = "centre"))
  col_names <- sapply(1:6, function(x) format(col_names[x], width = max_col_width[x], justify = "centre"))
  measures <- as.data.frame(do.call(cbind, measures))
  colnames(measures) <- col_names

  measures[negative] <- paste0(measures[negative], " ")
  measures <- format(measures, width = max(max_col_width), justify = "right")
  measures[, 5] <- paste0("   ", measures[, 5])
  rownames(measures) <- row_names

  # Final utils.
  spacing_base <- nchar(row_names[which.max(nchar(row_names))])
  line <- paste(rep("-", 120), collapse = "")
  line_thick <- paste(rep("=", 120), collapse = "")

  ## Printing.
  if (latex == FALSE) {
    cat(line_thick, "\n")
    cat(paste(rep(" ", spacing_base + 9)), "Treated               Controls              Overlap measures \n", sep = "")
    cat(paste(rep(" ", spacing_base + 7)), "(N_t = ", x$arm_sizes["treated"], ")      (N_c = ", x$arm_sizes["control"], ")      ------------------------- \n", sep = "")
    cat(line, "\n")
    print(measures, digits = 2, nsmall = 2, big.mark = ",")
    cat(line_thick, "\n")
  } else {
    table_names <- rename_latex(rownames(measures))

    cat("      \\begingroup
        \\setlength{\\tabcolsep}{8pt}
        \\renewcommand{\\arraystretch}{1.1}
        \\begin{table}[H]
          \\centering
          \\begin{adjustbox}{width = 0.75\\textwidth}
          \\begin{tabular}{@{\\extracolsep{5pt}}l c c c c c c}
          \\\\[-1.8ex]\\hline
          \\hline \\\\[-1.8ex]
          & \\multicolumn{2}{c}{Treated} & \\multicolumn{2}{c}{Controls} & \\multicolumn{2}{c}{Overlap measures} \\\\ \\cmidrule{6-7}
          & \\multicolumn{2}{c}{($n_t = ", x$arm_sizes["treated"], "$)} & \\multicolumn{2}{c}{($n_c =", x$arm_sizes["control"], "$)} & \\\\ \\cmidrule{2-5}
          & Mean & (S.D.) & Mean & (S.D.) & $\\hat{\\Delta}_j$ & $\\hat{\\Gamma}_j$ \\\\
          \\addlinespace[2pt]
          \\hline \\\\[-1.8ex] \n\n")

    for (i in seq_len(length(table_names))) {
      cat("            \\texttt{", table_names[i], "} & ", stringr::str_replace_all(string = measures[i, 1], pattern = " ", repl = ""),
          " & (", stringr::str_replace_all(string = measures[i, 2], pattern = " ", repl = ""),
          ") & ", stringr::str_replace_all(string = measures[i, 3], pattern = " ", repl = ""),
          " & (", stringr::str_replace_all(string = measures[i, 4], pattern = " ", repl = ""),
          ") & ", stringr::str_replace_all(string = measures[i, 5], pattern = " ", repl = ""),
          " & ", stringr::str_replace_all(string = measures[i, 6], pattern = " ", repl = ""), " \\\\ \n",
          sep = ""
      )
    }

    cat("\n\n          \\addlinespace[3pt]
          \\\\[-1.8ex]\\hline
          \\hline \\\\[-1.8ex]
          \\end{tabular}
          \\end{adjustbox}
          \\caption{https://soundcloud.com/theplasticchairband}
          \\label{table:descriptive.stats}
        \\end{table}
      \\endgroup")
  }
}
