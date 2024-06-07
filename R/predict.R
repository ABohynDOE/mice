#' Predict method to extract fitted values from a `mira` object
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
predict.mira <- function(object, newdata = NULL, ...){
  if (!inherits(object, "mira")){
    stop("object must be of the class 'mira'")
  }
  predm <- lapply(
    getfit(object), predict, newdata = newdata, se.fit = TRUE, ...
  )
  # Obtain predictions Q and prediction variances U, for m imputations
  Q <- sapply(predm, `[[`, "fit")
  U <- sapply(predm, `[[`, "se.fit")^2
  dfcom <- object$analyses[[1]]$df.residual
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
