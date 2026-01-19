#' Build mediator regression data for a given X
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom foreach foreach %do%
#' @param x_df_list named list of data frames for X (must include matching_key)
#' @param fixed_cov data frame of fixed covariates (optional, depends on your pipeline)
#' @param mediator_list list of mediator data frames (each must include matching_key)
#' @param x_name character, which element in x_df_list to use
#' @param num_cores integer, number of cores for parallel backend
#' @return data frame (row-bind results across matching keys)
#' @export
get_M_RegDat = function(x_df_list, fixed_cov, mediator_list, x_name, num_cores){
  x_interested_df = x_df_list[[x_name]]
  # x_m_common_keys: at least having one mediator for mediation analysis
  x_m_common_keys = intersect(
    x_interested_df$matching_key,
    unique(unlist(lapply(mediator_list, function(df) df[["matching_key"]]))))
  x_interested_df = x_interested_df[x_interested_df$matching_key %in% x_m_common_keys,]
  cl <- parallel::makeCluster(num_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  doParallel::registerDoParallel(cl)

  ii = NULL
  start_time <- Sys.time()
  final_res = foreach::foreach(
    ii = seq_len(nrow(x_interested_df)),
    .combine = rbind,
    .packages = c("stats","dplyr","iProPath"),
    .errorhandling = "pass"
  ) %dopar% {
           match_key_g = x_interested_df$matching_key[ii]
           x_interested = dropindex_and_transpose(x_interested_df[ii,])

           x_covs_df_list =x_df_list[setdiff(names(x_df_list),x_name)]
           x_covs = extract_data_g(x_covs_df_list, match_key_g)
           x_covs_value = dropindex_and_transpose(x_covs)

           m = extract_data_g(mediator_list, match_key_g)
           m_value = dropindex_and_transpose(m)

           reg_p = c(); reg_beta = c(); param_sd= c()
           for (i in 1:ncol(m_value)) {
             m_i = m_value[,i]
             x_reg_df = data.frame(x_interested, x_covs_value, fixed_cov)
             colnames(x_reg_df)[1] = x_name

             reg = summary(stats::lm(m_i ~ ., data = x_reg_df))
             reg_beta = c(reg_beta, reg$coefficients[x_name, 1])
             param_sd = c(param_sd, reg$coefficients[x_name, 2])
             reg_p = c(reg_p, reg$coefficients[x_name, 4])
           }
           l = ncol(m_value) # number of Mediators in this gene dataframe.
           data.frame(
             idx_Xgene = rep(match_key_g, l),
             idx_Mediator = paste("Mediator",rownames(m),m$matching_key,sep="_"),
             p_value_MX = reg_p,
             Est_param_MX = reg_beta,
             Std_Error_MX = param_sd
           )
  }
  return(final_res)
}

