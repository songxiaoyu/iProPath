#' Write iProPath outputs (xlsx + text summary)
#'
#' @param res_all Final results data frame (from \code{assemble_ipropath_results()}).
#' @param real_data_pi_df Data frame of \eqn{\pi} estimates from real data.
#' @param ipropath_pi_df Data frame of \eqn{\pi} estimates from iProPath.
#' @param file_stub Character string used to build output filenames (without extension).
#' @param output_dir Directory to write output files.
#' @return NULL. Files are written to disk.
#' @export
write_ipropath_outputs <- function(
    res_all,
    real_data_pi_df,
    ipropath_pi_df,
    file_stub,
    output_dir = "Results"
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  openxlsx::write.xlsx(
    list(
      "iProPath" = res_all,
      "Real_data.Nonnull Estimation" = real_data_pi_df,
      "iProPath.Nonnull Estimation" = ipropath_pi_df
    ),
    file = file.path(output_dir, sprintf("%s.iProPath.xlsx", file_stub)),
    rowNames = FALSE
  )

  summarize_ipropath_results(
    res_all = res_all,
    real_data.pi_df = real_data_pi_df,
    ipropath.pi_df = ipropath_pi_df,
    file_name = file_stub,
    output_dir = output_dir
  )

  invisible(NULL)
}
