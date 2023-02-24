#' Descriptive Statistics by Treatment Arm
#'
#' Computes sample averages and standard deviations of the covariates across treatment arms.
#'
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' 4xp array, storing the desired statistics.
#'
#' @keywords internal
#'
#' @details
#' Sample means and standard deviations across treatment arms are a first, useful insight to assess covariate balance.
descriptive_arm <- function(X, D) {
  ## Computing measures.
  mean_treated <- round(apply(X[D == 1, ], MARGIN = 2, mean), 3)
  sd_treated <- round(apply(X[D == 1, ], MARGIN = 2, stats::sd), 3)

  mean_control <- round(apply(X[D == 0, ], MARGIN = 2, mean), 3)
  sd_control <- round(apply(X[D == 0, ], MARGIN = 2, stats::sd), 3)

  ## Output.
  out <- rbind("Mean T" = mean_treated, "S.D. T" = sd_treated, "Mean C" = mean_control, "S.D. C" = sd_control)
  return(out)
}


#' Normalized Differences
#'
#' Computes a measure of the difference between locations of the covariate distributions
#' across treatment arms.
#'
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' 1xp data frame storing the normalized difference of each covariate.
#'
#' @keywords internal
#'
#' @details
#' Normalized differences are computed as the difference in the means of each covariate across treatment arms, normalized
#' by the sum of the within-arm variances.
normalized_diff <- function(X, D) {
  ## Computing measure.
  results <- round(apply(X, MARGIN = 2, function(x) {
    (mean(x[D == 1]) - mean(x[D == 0])) / sqrt((stats::var(x[D == 1]) + stats::var(x[D == 0])) / 2)
  }), 3)

  ## Handling output.
  temp_mat <- matrix(results, nrow = 1)
  colnames(temp_mat) <- colnames(X)
  rownames(temp_mat) <- "Norm. Diff."

  ## Output.
  out <- data.frame(temp_mat)
  return(out)
}


#' Log Ratio of Standard Deviations
#'
#' Computes a measure of the difference in the dispersion of the
#' covariate distributions across treatment arms.
#'
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' 1xp data frame storing logarithm of the ratio of standard deviations of each covariate.
#'
#' @keywords internal
#'
#' @details
#' Log ratio of standard deviations are computed as the logarithm of the ratio of the within-arm standard deviations.
log_ratio_sd <- function(X, D) {
  # Computing measure.
  results <- round(apply(X, MARGIN = 2, function(x) log(stats::sd(x[D == 1])) - log(stats::sd(x[D == 0]))), 3)

  # Handling output.
  temp_mat <- matrix(results, nrow = 1)
  colnames(temp_mat) <- colnames(X)
  rownames(temp_mat) <- "Log. S.D."

  # Output.
  out <- data.frame(temp_mat)
  return(out)
}


#' Balance Measures
#'
#' Compute several balance measures to check whether the covariate distributions are balanced across
#' treatment arms.
#'
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' Prints LATEX code in the console.
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
#' ## Print table.
#' balance_measures(X, D)
#'
#' @details
#' For each covariate in \code{X}, \code{balance_measures} computes sample averages and standard deviations
#' for both treatment arms. Additionally, two balance measures are computed:
#' \describe{
#'   \item{\code{Norm. Diff.}}{Normalized differences, computed as the differences in the means of each covariate
#'   across treatment arms, normalized by the sum of the within-arm variances. They provide a measure of the
#'   discrepancy between locations of the covariate distributions across treatment arms.}
#'   \item{\code{Log S.D.}}{Log ratio of standard deviations are computed as the logarithm of the ratio of the
#'   within-arm standard deviations. They provide a measure of the
#'   discrepancy in the dispersion of the covariate distributions across treatment arms.}
#'   }
#'
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#'
#' @author Elena Dal Torrione, Riccardo Di Francesco
#'
#' @export
balance_measures <- function(X, D) {
  ## Handling inputs.
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)

  names <- colnames(X)
  n_t <- sum(D)
  n_c <- length(D) - n_t

  ## Computing measures.
  descriptive <- descriptive_arm(X, D)
  normalized_diff <- normalized_diff(X, D)
  log_ratio_sd <- log_ratio_sd(X, D)

  ## Write table.
  x <- list("descriptive_stats" = descriptive, "norm_diff" = normalized_diff, "log_ratio_sd" = log_ratio_sd,
    "var_names" = names, "arm_sizes" = c("treated" = n_t, "control" = n_c))

  measures <- data.frame(cbind(t(x$descriptive_stats), t(x$norm_diff), t(x$log_ratio_sd)))
  col_names <- c("Mean", "S.D.", "Mean", "S.D.", "Norm.Diff.", "Log S.D.")
  row_names <- rownames(measures)

  # Arranging for printing.
  negative <- measures < 0

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

  table_names <- rename_latex(rownames(measures))

  cat("\\begingroup
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
    cat("    \\texttt{", table_names[i], "} & ", stringr::str_replace_all(string = measures[i, 1], pattern = " ", repl = ""),
        " & (", stringr::str_replace_all(string = measures[i, 2], pattern = " ", repl = ""),
        ") & ", stringr::str_replace_all(string = measures[i, 3], pattern = " ", repl = ""),
        " & (", stringr::str_replace_all(string = measures[i, 4], pattern = " ", repl = ""),
        ") & ", stringr::str_replace_all(string = measures[i, 5], pattern = " ", repl = ""),
        " & ", stringr::str_replace_all(string = measures[i, 6], pattern = " ", repl = ""), " \\\\ \n",
        sep = ""
    )
  }

  cat("\n    \\addlinespace[3pt]
    \\\\[-1.8ex]\\hline
    \\hline \\\\[-1.8ex]
    \\end{tabular}
    \\end{adjustbox}
    \\caption{Balance between treatment and control groups. The last two columns report the estimated normalized differences ($\\hat{\\Delta}_j$) and logarithms of the ratio of standard deviations ($\\hat{\\Gamma}_j$).}
    \\label{table:descriptive.stats}
    \\end{table}
\\endgroup")
}
