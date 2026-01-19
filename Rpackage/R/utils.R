#' Extract rows matching a given key from a list of data frames
#'
#' @param df_list list of data frames
#' @param match_key_g value of matching_key to extract
#' @return data frame
extract_data_g <- function(df_list, match_key_g) {
  do.call(rbind, Map(function(df) {
    df[df$matching_key == match_key_g, , drop = FALSE]
  }, df_list))
}

#' Drop index columns and transpose data frame
#'
#' @param df data frame
#' @return transposed matrix
dropindex_and_transpose <- function(df) {
  t(df[, setdiff(colnames(df), c("unique_index", "matching_key")),
       drop = FALSE])
}

#' Time a code block and print elapsed time
#'
#' A lightweight helper to measure and report the running time of an
#' expression using \code{message()}. Mainly intended for internal
#' progress/diagnostic messages in package functions.
#'
#' @param label A short character string used as the message label.
#' @param expr An expression to be evaluated.
#'
#' @return The result of evaluating \code{expr}.
#'
#' @keywords internal
.time_block <- function(label, expr) {
  t0 <- Sys.time()
  out <- force(expr)
  message(sprintf(
    "[%s] %.2f min",
    label,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
  out
}


