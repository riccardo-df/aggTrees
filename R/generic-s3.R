#' Plot Method for aggTrees Objects
#'
#' Plots an \code{aggTrees} object.
#'
#' @param x An \code{aggTrees} object.
#' @param leaves Number of leaves of the desired tree. This can be used to plot subtrees.
#' @param sequence If \code{TRUE}, the whole sequence of optimal partitions is displayed in a short animation.
#' @param ... Further arguments from \code{\link[rpart.plot]{prp}}.
#'
#' @details
#' Nodes are colored using a diverging palette. Nodes with predictions smaller than the ATE (i.e., the root
#' prediction) are colored in blue shades, and nodes with predictions larger than the ATE are colored in red
#' shades. Moreover, predictions that are more distant in absolute value from the ATE get darker shades.
#' This way, we have an immediate understanding of the groups with extreme GATEs.
#'
#' @import rpart.plot
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @seealso
#' \code{\link{build_aggtree}}, \code{\link{analyze_aggtree}}
#'
#' @export
plot.aggTrees <- function(x, leaves = get_leaves(x$tree), sequence = FALSE, ...) {
  ## Handling inputs and checks.
  if (!(inherits(x, "aggTrees"))) stop("You must provide a valid aggTrees object.", call. = FALSE)
  if (!(sequence %in% c(TRUE, FALSE))) stop("Invalid 'sequence'. This must be either TRUE or FALSE.", call. = FALSE)

  tree <- x$tree

  ## Plotting.
  if (sequence) {
    ## Node labels.
    prefix <- c("ATE = ", rep("GATE = ", times = (length(tree$frame$n) - 1)))
    sizes <- tree$frame$n
    percentages <- round(tree$frame$n / tree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    ## Threshold values of cost-complexity parameter and node numbers.
    nodes <- as.numeric(rownames(tree$frame))
    alpha_values <- rev(tree$cptable[, "CP"])

    ## Settings for node colors.
    low_colors <- c("#002F70", "#6889D0", "#F6F6F6")
    high_colors <- c("#F6F6F6", "#CB6F70", "#5F1415")
    color_fun_low <- colorRampPalette(low_colors)
    color_fun_high <- colorRampPalette(high_colors)

    low_palette <- color_fun_low(sum(tree$frame$yval <= tree$frame$yval[1]))
    high_palette <- color_fun_high(sum(tree$frame$yval >= tree$frame$yval[1]))[-1]

    ## Plot sequence.
    for (alpha in alpha_values[1:length(alpha_values[-1])]) {
      ## Graying out pruned leaves.
      temp_nodes <- as.numeric(rownames(rpart::prune.rpart(tree, alpha)$frame))
      colors <- ifelse(nodes %in% temp_nodes, 1, "white")

      ## Get rid of pruned leaves, branches and split labels by using with colors. Also, order colors as in Details.
      box_colors <- c(low_palette, high_palette)
      ate <- tree$frame$yval[1]
      to_order_bottom <- sort(tree$frame$yval[-1][tree$frame$yval[-1] < ate])
      to_order_top <- sort(tree$frame$yval[-1][tree$frame$yval[-1] > ate])
      to_order_box_colors <- c(to_order_bottom, ate, to_order_top)
      box_colors <- box_colors[match(tree$frame$yval, to_order_box_colors)]
      box_colors[colors == "white"] <- "white"

      temp_parents <- temp_nodes %/% 2
      split_colors <- ifelse(nodes %in% temp_parents, "black", "white")

      ## Plot.
      grDevices::dev.hold()

      rpart.plot::prp(tree,
                      type = 2,
                      extra = 0,
                      under = FALSE,
                      fallen.leaves = TRUE,
                      round = 0,
                      leaf.round = 0,
                      prefix = prefix,
                      suffix = suffix,
                      box.col = box_colors,
                      branch = 0.3,
                      branch.col = colors,
                      split.col = split_colors,
                      col = colors,
                      ...)

      grDevices::dev.flush()
      Sys.sleep(1)
    }
  } else {
    ## Get subtree.
    subtree <- subtree(tree, leaves)

    ## Node labels.
    prefix <- c("ATE = ", rep("GATE = ", times = (length(subtree$frame$n) - 1)))
    sizes <- subtree$frame$n
    percentages <- round(subtree$frame$n / subtree$frame$n[1] * 100, 0)
    suffix <- paste("\n", "n = ", sizes, "\n (", percentages, "%)", sep = "")

    ## Settings for node colors.
    low_colors <- c("#002F70", "#6889D0", "#F6F6F6")
    high_colors <- c("#F6F6F6", "#CB6F70", "#5F1415")
    color_fun_low <- colorRampPalette(low_colors)
    color_fun_high <- colorRampPalette(high_colors)

    low_palette <- color_fun_low(sum(subtree$frame$yval <= subtree$frame$yval[1]))
    high_palette <- color_fun_high(sum(subtree$frame$yval >= subtree$frame$yval[1]))[-1]
    box_colors <- c(low_palette, high_palette)

    ate <- subtree$frame$yval[1]
    to_order_bottom <- sort(subtree$frame$yval[-1][subtree$frame$yval[-1] < ate])
    to_order_top <- sort(subtree$frame$yval[-1][subtree$frame$yval[-1] > ate])
    to_order_box_colors <- c(to_order_bottom, ate, to_order_top)
    box_colors <- box_colors[match(subtree$frame$yval, to_order_box_colors)]

    ## Plot.
    rpart.plot::prp(subtree,
                    type = 2,
                    extra = 0,
                    under = FALSE,
                    fallen.leaves = TRUE,
                    round = 0,
                    leaf.round = 0,
                    prefix = prefix,
                    suffix = suffix,
                    box.col = box_colors,
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
#' \code{\link{build_aggtree}}, \code{\link{analyze_aggtree}}
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
summary.aggTrees <- function(object, ...) {
  if (is.null(object$idx$honest_idx)) cat("Honest estimates:", FALSE, "\n") else cat("Honest estimates:", TRUE, "\n")
  summary(object$tree)
}


#' Print Method for aggTrees Objects
#'
#' Prints an \code{aggTrees} object.
#'
#' @param x \code{aggTrees} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' \code{\link{build_aggtree}}, \code{\link{analyze_aggtree}}
#'
#' @references
#' \itemize{
#'   \item R Di Francesco (2022). Aggregation Trees. CEIS Research Paper, 546. \doi{10.2139/ssrn.4304256}.
#' }
#'
#' @author Riccardo Di Francesco
#'
#' @export
print.aggTrees <- function(x, ...) {
  if (is.null(x$idx$honest_idx)) cat("Honest estimates:", FALSE, "\n") else cat("Honest estimates:", TRUE, "\n")
  print(x$tree)
}
