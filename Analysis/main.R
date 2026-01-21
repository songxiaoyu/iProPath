library(iProPath)
data("sample_data")

mediator_list = sample_data$mediator_list
x_df_list = sample_data$x_df_list
Y_list = sample_data$Y_list

fixed_cov = NULL
y_covdf_list = NULL
# args <- list("melanoma_final_v1", "6", "trans_glyco", "mutation" )
# Y_list$trans_glyco = Y_list$trans_glyco[sample(100),]

args <- list("sample_data", "6", "ptm", "genetic_feature" )
file_name <- unlist(args[1])
num_cores <- as.numeric(args[2])
y_name <- unlist(args[3])
x_name <- unlist(args[4])



reg_res <- prepare_ipropath_regression(
  x_df_list, fixed_cov, mediator_list, y_covdf_list, Y_list,
  x_name = x_name,
  y_name = y_name,
  num_cores = num_cores,
  adjust_methods = c("BH", "bonferroni", "BY")
)

pi_res <- estimate_ipropath_pi(
  reg_res$merged_df$p_value_MX,
  reg_res$merged_df$p_value_YM,
  method = "qvalue.smoother",
  run_all = FALSE
)

ipropath_res <- run_ipropath(
  reg_res$merged_df,
  pi_res$pi1_MX,
  pi_res$pi1_YM,
  num_cores = num_cores,
  y_name = y_name
)

res_all <- assemble_ipropath_results(
  pairwise_res = reg_res$pairwise_res,
  y_transreg = reg_res$y_transreg,
  ipropath_res_p = ipropath_res$res_p
)

file_stub <- paste(file_name, x_name, y_name, sep = ".")
write_ipropath_outputs(
  res_all = res_all,
  real_data_pi_df = pi_res$pi_df,
  ipropath_pi_df = ipropath_res$pi_df,
  file_stub = file_stub,
  output_dir = "Results"
)
