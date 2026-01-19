#' Prepare regression results for iProPath
#'
#' Runs trans regression and pairwise regression, then returns the merged
#' regression table used by iProPath as well as pairwise regression results.
#' @importFrom rlang .data
#' @param x_df_list List of X data frames (as used by your pipeline).
#' @param fixed_cov Fixed covariates data frame/matrix (as used by your pipeline).
#' @param mediator_list List of mediators (as used by your pipeline).
#' @param y_covdf_list List of Y covariates data frames (as used by your pipeline).
#' @param Y_list List of outcome data frames (as used by your pipeline).
#' @param x_name Character string specifying X type/name (e.g., \code{"mutation"}).
#' @param y_name Character string specifying Y type/name (e.g., \code{"glyco"}).
#' @param num_cores Integer number of CPU cores to use.
#' @param adjust_methods Character vector of p-value adjustment methods to compute
#'   for pairwise results.
#' @return A list with components:
#' \describe{
#'   \item{merged_df}{Data frame of merged trans-regression results.}
#'   \item{y_transreg}{Data frame of Y trans-regression results (pre-merge).}
#'   \item{m_reg}{Data frame of M regression results (pre-merge).}
#'   \item{pairwise_res}{Data frame of pairwise regression results with adjusted p-values.}
#' }
#' @export
prepare_ipropath_regression <- function(
    x_df_list,
    fixed_cov,
    mediator_list,
    y_covdf_list,
    Y_list,
    x_name = "mutation",
    y_name,
    num_cores = 1L,
    adjust_methods = c("BH", "bonferroni", "BY")
) {
  num_cores <- as.integer(num_cores)
  if (is.na(num_cores) || num_cores <= 0L) stop("num_cores must be a positive integer.", call. = FALSE)
  # Module 1: Trans regression
  m_reg <- .time_block("get_M_RegDat",{get_M_RegDat(x_df_list, fixed_cov, mediator_list, x_name, num_cores)})
  y_transreg <- .time_block("get_Y_transRegDat",{get_Y_transRegDat(x_df_list, fixed_cov, y_covdf_list, mediator_list, Y_list, x_name, y_name, num_cores)})
  merged_df <- dplyr::left_join(y_transreg, m_reg, by = c("idx_Xgene", "idx_Mediator"))
  # Module 2: Pairwise regression
  pairwise_res <-.time_block("get_pairwise_transRegDat", {get_pairwise_transRegDat(x_df_list, fixed_cov, y_covdf_list, Y_list, x_name, y_name, num_cores)})

  # Add adjusted p-values
  if (!is.null(adjust_methods) && length(adjust_methods) > 0) {
    if ("BH" %in% adjust_methods) {
      pairwise_res <- dplyr::mutate(pairwise_res, Pairwise_BH = stats::p.adjust(.data$Pairwise_p_value, method = "BH"))
    }
    if ("bonferroni" %in% adjust_methods) {
      pairwise_res <- dplyr::mutate(pairwise_res, Pairwise_FWER = stats::p.adjust(.data$Pairwise_p_value, method = "bonferroni"))
    }
    if ("BY" %in% adjust_methods) {
      pairwise_res <- dplyr::mutate(pairwise_res, Pairwise_BY = stats::p.adjust(.data$Pairwise_p_value, method = "BY"))
    }
  }

  list(
    merged_df = merged_df,
    y_transreg = y_transreg,
    m_reg = m_reg,
    pairwise_res = pairwise_res
  )
}
