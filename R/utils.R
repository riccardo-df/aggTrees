#' Renaming Variables for LATEX Usage (Internal Use)
#'
#' Renames variables where the character "_" is used, which causes clashes in LATEX. Useful for the \code{phased} print method.
#'
#' @param names string vector.
#'
#' @return
#' The renamed string vector. Strings where "_" is not found are not modified by \code{rename_latex}.
rename_latex <- function(names) {
  ## Locating variables that need renaming.
  idx <- grepl("_", names, fixed = TRUE)

  ## Renaming variables.
  split_names <- stringr::str_split(string = names[idx], pattern = "_", simplify = TRUE)
  attach_names <- paste(split_names[, 1], split_names[, 2], sep = "\\_")

  ## Replacing.
  names[idx] <- attach_names

  ## Output.
  return(names)
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
#' It is assumed that categorical variables are coded as \code{factors}. Also, no missing values are allowed.
#'
#' @export
expand_df <- function(X, int_order = 2, poly_order = 4, threshold = 0) {
  ## Handling inputs and checks.
  if (int_order < 1 | int_order > 4) stop("Wrong order of interactions! Must be either 1, 2, 3 or 4.")
  if (poly_order < 1) stop("Wrong order of polynomials! Must be either greater than zero.")
  if (threshold < 0 | threshold > 1) stop("Wrong threshold! Must lie in the open interval (0, 1).")
  if (threshold == 1) stop("Wrong threshold! Cannot drop all the units.")

  X <- as.data.frame(X)
  X <- stats::model.matrix(~., data = data.frame(X))[, -1]

  X_continuous <- X[, !apply(X, MARGIN = 2, function(x) all(x %in% 0:1))]

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
    for (i in seq_len(dim(X_continuous)[2]))
    {
      temp.poly <- stats::poly(X_continuous[, i], degree = poly_order, raw = TRUE)[, -1]
      expanded_X <- data.frame(expanded_X, temp.poly)

      for (j in 2:poly_order)
      {
        colnames(expanded_X)[(dim(expanded_X)[2]) - poly_order + j] <- paste(paste(colnames(X_continuous)[i], "..", sep = ""), j, sep = "")
      }
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
