#' Summarize iProPath results and write to file
#'
#' @importFrom rlang .data
#' @param res_all ...
#' @param real_data.pi_df ...
#' @param ipropath.pi_df ...
#' @param file_name ...
#' @param output_dir ...
#' @return NULL. Results are written to a text file.
#' @export

summarize_ipropath_results <- function(res_all, real_data.pi_df=NULL, ipropath.pi_df, file_name, output_dir = "Results") {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # Define output file path
  output_file <- sprintf("%s/%s.txt", output_dir, file_name)
  # Redirect output to the file
  n0 <- sink.number()
  sink(output_file, append = FALSE)
  on.exit({
    while (sink.number() > n0) sink()
  }, add = TRUE)
  print(sprintf("Filename: %s",file_name))
  cat("\n")
  # Summary statistics
  cat("Total rows:", nrow(res_all),
      " Total genes:", length(unique(res_all$X_gene_name)),
      " Total sites:", length(unique(res_all$Y_site)), "\n")
  res_all=res_all %>% dplyr::mutate(Pairwise_BY=stats::p.adjust(.data$Pairwise_p_value,method="BY"))
  res_all=res_all %>% dplyr::mutate(iProPath_BY=stats::p.adjust(.data$iProPath_p_value,method="BY"))
  cat("\n")
  # Display pi estimation
  print("real_data.pi_df:")
  print(real_data.pi_df)
  print("ipropath.pi_df:")
  print(ipropath.pi_df)
  # significant x_gene counts
  cat("\n")
  cat("Number of Pairwise FWER <= 0.05:", res_all %>% dplyr::filter(.data$Pairwise_FWER <= 0.05) %>% nrow(), "\n")
  cat("Number of Pairwise BH <= 0.1:", res_all %>% dplyr::filter(.data$Pairwise_BH <= 0.1) %>% nrow(), "\n")
  cat("Number of Pairwise BY <= 0.1:", res_all %>% dplyr::filter(.data$Pairwise_BY <= 0.1) %>% nrow(), "\n")
  cat("Number of iProPath FWER <= 0.05:", res_all %>% dplyr::filter(.data$iProPath_FWER <= 0.05) %>% nrow(), "\n")
  cat("Number of iProPath BH <= 0.1:", res_all %>% dplyr::filter(.data$iProPath_BH <= 0.1) %>% nrow(), "\n")
  cat("Number of iProPath BY <= 0.1:", res_all %>% dplyr::filter(.data$iProPath_BY <= 0.1) %>% nrow(), "\n")
  cat("\n")
  # Top hits display
  print("Pairwise FWER <= 0.05 Hubs:")
  print(res_all %>% dplyr::filter(.data$Pairwise_FWER <= 0.05) %>% dplyr::count(.data$X_gene_name, sort = TRUE))
  print("iProPath FWER <= 0.05 Hubs:")
  print(res_all %>% dplyr::filter(.data$iProPath_FWER <= 0.05) %>% dplyr::count(.data$X_gene_name, sort = TRUE))

  print("Pairwise BH <= 0.1 Hubs:")
  print(res_all %>% dplyr::filter(.data$Pairwise_BH <= 0.1) %>% dplyr::count(.data$X_gene_name, sort = TRUE))
  print("iProPath BH <= 0.1 Hubs:")
  print(res_all %>% dplyr::filter(.data$iProPath_BH <= 0.1) %>% dplyr::count(.data$X_gene_name, sort = TRUE))

  print("Pairwise BY <= 0.1 Hubs:")
  print(res_all %>% dplyr::filter(.data$Pairwise_BY <= 0.1) %>% dplyr::count(.data$X_gene_name, sort = TRUE))
  print("iProPath BY <= 0.1 Hubs:")
  print(res_all %>% dplyr::filter(.data$iProPath_BY <= 0.1) %>% dplyr::count(.data$X_gene_name, sort = TRUE))
  return(NULL)
}
