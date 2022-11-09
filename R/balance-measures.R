#' Descriptive Statistics by Treatment Arm (Internal Use)
#'
#' Computes sample averages and standard deviations of the covariates across treatment arms.
#'
#' @param X Covariate matrix (no intercept).
#' @param D Treatment assignment vector.
#'
#' @return
#' 4xp array, storing the desired statistics.
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


#' Normalized Differences (Internal Use)
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


#' Log Ratio of Standard Deviations (Internal Use)
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
#' A \code{phased} object.
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
#' @author Elena Dal Torrione, Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item G. W. Imbens, D. B. Rubin (2015). Causal inference in statistics, social, and biomedical sciences. Cambridge University Press. \doi{10.1017/CBO9781139025751}.
#' }
#'
#' @seealso \code{\link{print.phased}}
#'
#' @export
balance_measures <- function(X, D) {
  ## Handling inputs.
  names <- colnames(X)
  n_t <- sum(D)
  n_c <- length(D) - n_t

  ## Computing measures.
  descriptive <- descriptive_arm(X, D)
  normalized_diff <- normalized_diff(X, D)
  log_ratio_sd <- log_ratio_sd(X, D)

  ## Handling output.
  out <- list(
    "descriptive_stats" = descriptive, "norm_diff" = normalized_diff, "log_ratio_sd" = log_ratio_sd,
    "var_names" = names, "arm_sizes" = c("treated" = n_t, "control" = n_c)
  )
  class(out) <- "phased"

  ## Output.
  return(out)
}


