#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics glance
#' @export
generics::augment

#' Tidy method to extract results from a `mipo` object
#'
#' @param x An object of class \code{mipo}
#' @param conf.int Logical. Should confidence intervals be returned?
#' @param conf.level Confidence level for intervals. Defaults to .95
#' @param ... extra arguments (not used)
#' @export
#' @keywords internal
#' @return A dataframe withh these columns:
#' \itemize{
#'      \item term
#'      \item estimate
#'      \item ubar
#'      \item b
#'      \item t
#'      \item dfcom
#'      \item df
#'      \item riv
#'      \item lambda
#'      \item fmi
#'      \item p.value
#'      \item conf.low (if called with conf.int = TRUE)
#'      \item conf.high (if called with conf.int = TRUE)
#' }
tidy.mipo <- function(x, conf.int = FALSE, conf.level = .95, ...) {
  out <- summary(x,
    type = "all",
    conf.int = conf.int,
    conf.level = conf.level,
    ...
  )

  if ("term" %in% names(out)) out$term <- as.character(out$term)
  if ("contrast" %in% names(out)) out$contrast <- as.character(out$contrast)

  # needed for broom <= 0.5.6
  # rename variables if present
  idx <- grepl("%", names(out))
  names(out)[idx] <- c("conf.low", "conf.high")
  idx <- names(out) == "t"
  names(out)[idx] <- "statistic"

  # order columns
  cols_a <- c(
    "term", "estimate", "std.error", "statistic", "p.value",
    "conf.low", "conf.high"
  )
  cols_a <- base::intersect(cols_a, colnames(out))
  cols_b <- sort(base::setdiff(colnames(out), cols_a))
  out[, c(cols_a, cols_b)]
}

#' Glance method to extract information from a `mipo` object
#'
#' @param x An object with multiply-imputed models from `mice` (class: `mipo`)
#' @param ... extra arguments (not used)
#' @return a dataframe with one row and the following columns:
#' \itemize{
#'   \item nimp
#'   \item nobs
#' }
#' @note If x contains `lm` models, R2 and Adj.R2 are included in the output
#' @export
#' @keywords internal
#' @family tidiers
glance.mipo <- function(x, ...) {
  out <- data.frame(nimp = nrow(x$glanced))
  out$nobs <- tryCatch(x$glanced$nobs[1],
    error = function(e) NULL
  )

  # R2 in lm models
  out$r.squared <- tryCatch(pool.r.squared(x, adjusted = FALSE)[1],
    error = function(e) NULL
  )
  out$adj.r.squared <- tryCatch(pool.r.squared(x, adjusted = TRUE)[1],
    error = function(e) NULL
  )

  out
}

#' Augment method to extract fitted values from a `mira` object
#'
#' @param x A \code{mira} object with multiply-imputed models created from `mice`
#' @param newdata An optional [`base::data.frame()`] or [`tibble::tibble()`] in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param ... extra arguments (not used)
#' @return a dataframe with one row and the following columns:
#' \itemize{
#'   \item `.fitted` the pooled fitted values
#'   \item `.se.fit` the pooled standard error of the fitted values
#'   \item `.dffits` the degrees of freedom associated with the pooled fitted values
#' }
#' @export
#' @keywords internal
#' @family tidiers
augment.mira <- function(x, newdata = NULL, ...){
  if (!inherits(object, "mira")){
    stop("object must be of the class 'mira'")
  }
  predm <- lapply(
    getfit(object), predict, newdata = newdata, se.fit = TRUE, ...
  )
  # Obtain predictions Q and prediction variances U, for m imputations
  Q <- sapply(predm, `[[`, "fit")
  U <- sapply(predm, `[[`, "se.fit")^2
  dfcom <- fit$analyses[[1]]$df.residual
  # Pool predictions
  pred <- matrix(NA, nrow = nrow(Q), ncol = 3,
                 dimnames = list(NULL, c(".fitted", ".se.fit", ".dffits")))
  for(i in 1:nrow(Q)) {
    pi <- pool.scalar(Q[i, ], U[i, ], n = dfcom + 1)
    pred[i, 1] <- pi[["qbar"]]
    pred[i, 2] <- sqrt(pi[["t"]])
    pred[i, 3] <- pi[["df"]]
  }
  pred
}
