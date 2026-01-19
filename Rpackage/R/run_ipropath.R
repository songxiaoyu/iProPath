#' Run iProPath core algorithm
#'
#' @param merged_df Data frame containing columns required by \code{ipropath()}:
#'   p_value_YX, p_value_MX, p_value_YM, idx_Xgene, idx_Mediator, idx_Ygene, idx_Ysite,
#'   Est_param_YX, Est_param_YM, Est_param_MX.
#' @param pi1_MX Estimated non-null proportion for M~X.
#' @param pi1_YM Estimated non-null proportion for Y~M.
#' @param num_cores Integer number of CPU cores to use.
#' @param y_name Character string passed to \code{ipropath()}.
#' @return A list returned by \code{ipropath()}, with an additional element
#'   \code{runtime_secs} giving the elapsed time in seconds.
#' @export
run_ipropath <- function(
    merged_df,
    pi1_MX,
    pi1_YM,
    num_cores = 1L,
    y_name
) {
  num_cores <- as.integer(num_cores)
  if (is.na(num_cores) || num_cores <= 0L) stop("num_cores must be a positive integer.", call. = FALSE)

  needed <- c(
    "p_value_YX","p_value_MX","p_value_YM",
    "idx_Xgene","idx_Mediator","idx_Ygene","idx_Ysite",
    "Est_param_YX","Est_param_YM","Est_param_MX"
  )
  miss <- setdiff(needed, colnames(merged_df))
  if (length(miss) > 0) {
    stop("merged_df is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  res <- .time_block("ipropath",{iProPath::ipropath(
          p_value_YX = merged_df$p_value_YX,
          p_value_MX = merged_df$p_value_MX,
          p_value_YM = merged_df$p_value_YM,
          pi1_MX = pi1_MX,
          pi1_YM = pi1_YM,
          idx_Xgene = merged_df$idx_Xgene,
          idx_Mediator = merged_df$idx_Mediator,
          idx_Ygene = merged_df$idx_Ygene,
          idx_Ysite = merged_df$idx_Ysite,
          Est_param_YX = merged_df$Est_param_YX,
          Est_param_YM = merged_df$Est_param_YM,
          Est_param_MX = merged_df$Est_param_MX,
          numCores = num_cores,
          y_name = y_name
  )})
  res
}
