#' Assemble final iProPath results table
#'
#' @importFrom rlang .data
#' @param pairwise_res Data frame of pairwise regression results.
#' @param y_transreg Data frame returned by \code{get_Y_transRegDat()}.
#' @param ipropath_res_p Data frame \code{res_p} from \code{ipropath()} output.
#' @param adjust_methods Character vector specifying multiple-testing corrections for iProPath p-values.
#' @return A data frame containing merged pairwise and iProPath results.
#' @export
assemble_ipropath_results <- function(
    pairwise_res,
    y_transreg,
    ipropath_res_p,
    adjust_methods = c("BH", "bonferroni", "BY")
) {
  pairwise_res1 <- pairwise_res %>%
    dplyr::left_join(y_transreg, by = c("idx_Xgene" = "idx_Xgene", "idx_Ysite" = "idx_Ysite")) %>%
    dplyr::rename(
      StdErr_X_in_YMXreg = .data$Std_Error_YX,
      StdErr_M_in_YMXreg = .data$Std_Error_YM
    )

  res_all <- pairwise_res1 %>%
    dplyr::left_join(ipropath_res_p, by = c("idx_Xgene" = "idx_Xgene", "idx_Ysite" = "idx_Ysite")) %>%
    dplyr::rename(
      X_gene_name = .data$idx_Xgene,
      Y_gene_name = .data$idx_Ygene.x,
      Y_site = .data$idx_Ysite,
      idx_Mediator = .data$idx_Mediator.x
    ) %>%
    dplyr::mutate(
      iProPath_p_value = dplyr::coalesce(.data$iProPath_p_value, .data$Pairwise_p_value),
      iProPath_Est_param = dplyr::coalesce(.data$iProPath_Est_param, .data$Pairwise_Est_param)
    )

  if (!is.null(adjust_methods) && length(adjust_methods) > 0) {
    if ("BH" %in% adjust_methods) {
      res_all <- dplyr::mutate(res_all, iProPath_BH = stats::p.adjust(.data$iProPath_p_value, method = "BH"))
    }
    if ("bonferroni" %in% adjust_methods) {
      res_all <- dplyr::mutate(res_all, iProPath_FWER = stats::p.adjust(.data$iProPath_p_value, method = "bonferroni"))
    }
    if ("BY" %in% adjust_methods) {
      res_all <- dplyr::mutate(res_all, iProPath_BY = stats::p.adjust(.data$iProPath_p_value, method = "BY"))
    }
  }

  res_all %>%
    dplyr::select(
      .data$X_gene_name, .data$Y_gene_name, .data$Y_site,
      .data$Pairwise_p_value, .data$Pairwise_Est_param, .data$Pairwise_Std_Error,
      dplyr::any_of(c("Pairwise_BH","Pairwise_FWER","Pairwise_BY")),
      .data$idx_Mediator, .data$Path_YX_p_value, .data$Path_YMX_p_value,
      .data$StdErr_X_in_YMXreg, .data$StdErr_M_in_YMXreg,
      .data$iProPath_p_value, .data$iProPath_Est_param,
      dplyr::any_of(c("iProPath_BH","iProPath_FWER","iProPath_BY"))
    )
}
