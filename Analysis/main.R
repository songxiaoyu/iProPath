# Manually define arguments (fixed order)
args <- list(
  "sample_data",        # file_name
  6,                    # num_cores
  "ptm",                # y_name
  "genetic_feature"    # x_name
)
file_name <- args[[1]]
num_cores <- as.numeric(args[[2]])
y_name <- args[[3]]
x_name <- args[[4]]

library(iProPath)
data("sample_data")
data=get(file_name)

x_df_list=data$x_df_list
fixed_cov=data$fixed_cov
mediator_list=data$mediator_list
y_covdf_list=data$y_covdf_list
Y_list = data$Y_list
x_name = args[[4]]
y_name =  args[[3]]
num_cores = args[[2]]
file_name = args[[1]]

reg_res <- prepare_ipropath_regression(
  x_df_list,
  fixed_cov,
  mediator_list,
  y_covdf_list,
  Y_list,
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



