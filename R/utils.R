#' Sample Splitting
#'
#' Splits the sample into training and honest subsamples.
#'
#' @param n Size of the sample to be split.
#' @param training_frac Fraction of units for the training sample.
#'
#' @return
#' A list storing the indexes for the two different subsamples.
#'
#' @author Riccardo Di Francesco
#'
#' @export
sample_split <- function(n, training_frac = 0.5) {
  ## Handling inputs and checks.
  if (n <= 0) stop("'n' cannot be equal or lower than zero.", call. = FALSE)
  if (training_frac <= 0 | training_frac > 1) stop("'training_frac' must lie in the interval (0, 1].", call. = FALSE)

  ## Split the sample.
  training_idx <- sample(1:n, floor(training_frac * n), replace = FALSE)
  if (training_frac < 1) honest_idx <- setdiff(1:n, training_idx) else if (training_frac == 1) honest_idx <- NULL

  ## Output.
  return(list("training_idx" = training_idx, "honest_idx" = honest_idx))
}


#' Doubly-Robust Scores
#'
#' Constructs doubly-robust scores via K-fold cross-fitting.
#'
#' @param Y Outcome vector.
#' @param D Treatment assignment vector.
#' @param X Covariate matrix (no intercept).
#' @param k Number of folds.
#'
#' @return
#' A vector of scores.
#'
#' @details
#' Honest regression forests are used to estimate the propensity score and the conditional mean function of the outcome.
#'
#' @importFrom caret createFolds
#' @import grf
#'
#' @author Riccardo Di Francesco
#'
#' @export
dr_scores <- function(Y, D, X, k = 5) {
  ## Handling inputs and checks.
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)
  if (k %% 1 != 0 | k < 2) stop("Invalid 'k'. This must be an integer greater than or equal to 2.", call. = FALSE)

  ## Generate folds.
  folds <- caret::createFolds(Y, k = k)

  ## Cross-fitting nuisances.
  nuisances_mat <- matrix(NA, nrow = length(Y), ncol = 3)
  colnames(nuisances_mat) <- c("pscore", "mu_0", "mu_1")

  for (fold in seq_len(length(folds))) {
    ## Leave one fold out and build forests on other folds.
    left_out_idx <- folds[[fold]]

    cond_mean_forest <- grf::regression_forest(data.frame("D" = D[-left_out_idx], X[-left_out_idx, ]), Y[-left_out_idx])
    pscore_forest <- grf::regression_forest(matrix(X[-left_out_idx, ], ncol = dim(X)[2]), D[-left_out_idx])

    ## Predict on left-out fold.
    nuisances_mat[left_out_idx, "mu_0"] <- predict(cond_mean_forest, data.frame("D" = rep(0, length(left_out_idx)), X[left_out_idx, ]))$predictions
    nuisances_mat[left_out_idx, "mu_1"] <- predict(cond_mean_forest, data.frame("D" = rep(1, length(left_out_idx)), X[left_out_idx, ]))$predictions
    nuisances_mat[left_out_idx, "pscore"] <- predict(pscore_forest, matrix(X[left_out_idx, ], ncol = dim(X)[2]))$predictions
  }

  ## Construct doubly-robust scores.
  scores <- nuisances_mat[, "mu_1"] - nuisances_mat[, "mu_0"] + (D * (Y - nuisances_mat[, "mu_1"])) / (nuisances_mat[, "pscore"]) -
    ((1 - D) * (Y - nuisances_mat[, "mu_0"])) / (1 - nuisances_mat[, "pscore"])

  ## Output.
  return(scores)
}


#' Covariate Matrix Expansion
#'
#' Expands the covariate matrix, adding interactions and polynomials. This is particularly useful for penalized regressions.
#'
#' @param X Covariate matrix (no intercept).
#' @param int_order Order of interactions to be added. Set equal to one if no interactions are desired.
#' @param poly_order Order of the polynomials to be added. Set equal to one if no polynomials are desired.
#' @param threshold Drop binary variables representing less than \code{threshold}\% of the population. Useful to speed up computation.
#'
#' @return
#' The expanded covariate matrix, as a data frame.
#'
#' @details
#' \code{expand_df} assumes that categorical variables are coded as \code{factors}. Also, no missing values are allowed.\cr
#'
#' \code{expand_df} uses \code{\link[stats]{model.matrix}} to expand factors to a set of dummy variables. Then, it identifies continuous covariates as those
#' not having 0 and 1 as unique values.\cr
#'
#' \code{expand_df} first introduces all the \code{int_order}-way interactions between the variables (using the expanded set of dummies), and then adds
#' \code{poly_order}-order polynomials for continuous covariates.
#'
#' @author Riccardo Di Francesco
#'
#' @export
expand_df <- function(X, int_order = 2, poly_order = 4, threshold = 0) {
  ## Handling inputs and checks.
  if (int_order < 1 | int_order > 4) stop("Wrong order of interactions! Must be either 1, 2, 3 or 4.")
  if (poly_order < 1) stop("Wrong order of polynomials! Must be greater than zero.")
  if (threshold < 0 | threshold > 1) stop("Wrong threshold! Must lie in the open interval (0, 1).")
  if (threshold == 1) stop("Wrong threshold! Cannot drop all the units.")

  X <- as.data.frame(X)
  X <- stats::model.matrix(~., data = data.frame(X))[, -1]

  idx_continuous <- !apply(X, MARGIN = 2, function(x) all(x %in% 0:1))
  X_continuous <- X[, idx_continuous]

  ## Adding int_order-way interactions.
  if (int_order == 1) {
    expanded_X <- X
  } else if (int_order == 2) {
    expanded_X <- stats::model.matrix(~ .^2, data = data.frame(X))[, -1]
  } else if (int_order == 3) {
    expanded_X <- stats::model.matrix(~ .^3, data = data.frame(X))[, -1]
  } else if (int_order == 4) {
    expanded_X <- stats::model.matrix(~ .^4, data = data.frame(X))[, -1]
  }

  ## Adding polynomials for continuous variables (works on original continuous covariates).
  if (poly_order > 1) {
    if (sum(idx_continuous) > 1) {
      for (i in seq_len(dim(X_continuous)[2])) {
        temp.poly <- stats::poly(X_continuous[, i], degree = poly_order, raw = TRUE)[, -1]
        expanded_X <- data.frame(expanded_X, temp.poly)

        for (j in 2:poly_order) {
          colnames(expanded_X)[(dim(expanded_X)[2]) - poly_order + j] <- paste(paste(colnames(X_continuous)[i], "..", sep = ""), j, sep = "")
        }
      }
    } else  if (sum(idx_continuous) == 1) {
      temp.poly <- stats::poly(X_continuous, degree = poly_order, raw = TRUE)[, -1]
      colnames(temp.poly) <- paste0(colnames(X)[idx_continuous], "..", 2:poly_order)
      expanded_X <- data.frame(expanded_X, temp.poly)
    }
  }

  ## Dropping binary variables with low variability.
  # Bit tricky: the following is an index which equals TRUE iff the column of expanded_X is binary and with low variability.
  temp_idx <- apply(expanded_X, MARGIN = 2, function(x) all(x %in% 0:1)) * (apply(expanded_X, MARGIN = 2, mean) < threshold)
  expanded_X <- expanded_X[, !temp_idx]

  ## Handling output.
  out <- data.frame(expanded_X)

  ## Output.
  return(out)
}


#' Renaming Variables for LATEX Usage
#'
#' Renames variables where the character "_" is used, which causes clashes in LATEX.
#'
#' @param names string vector.
#'
#' @return
#' The renamed string vector. Strings where "_" is not found are not modified by \code{rename_latex}.
#'
#' @keywords internal
rename_latex <- function(names) {
  ## Locating variables that need renaming.
  idx <- grepl("_", names, fixed = TRUE)

  if (sum(idx) == 0) return(names)

  ## Renaming variables.
  split_names <- stringr::str_split(string = names[idx], pattern = "_", simplify = TRUE)
  attach_names <- paste(split_names[, 1], split_names[, 2], sep = "\\_")

  ## Replacing.
  names[idx] <- attach_names

  ## Output.
  return(names)
}
