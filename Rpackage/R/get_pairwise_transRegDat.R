#' Build pairwise trans-regression data for Y ～ X
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom foreach foreach %do%
#' @param x_df_list named list of data frames for X
#' @param fixed_cov data frame of fixed covariates (optional, depends on your pipeline)
#' @param y_covdf_list data frame of covariates associated with Y (optional, depends on your pipeline)
#' @param Y_list named list of outcome (Y) data frames.
#' @param x_name character, which element in x_df_list to use
#' @param y_name character, name of the outcome in Y_list
#' @param num_cores integer, number of cores for parallel backend
#' @return data frame (row-bind results across matching keys)
#' @export
get_pairwise_transRegDat <- function(x_df_list, fixed_cov, y_covdf_list, Y_list,x_name,y_name,num_cores) {
  if (!dir.exists("Progress")) dir.create("Progress")
  sinkfile=paste0("Progress/",y_name,"_Pairwise_transreg_progress.txt")
  cat(sprintf(">> File '%s' has been created. Check this file for progress updates.\n", sinkfile), file=sinkfile, append=TRUE)

  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  Y=Y_list[[y_name]]
  x_interested_df = x_df_list[[x_name]]
  ii = NULL
  start_time <- Sys.time()
  transreg_results<-foreach::foreach(ii = seq_len(nrow(x_interested_df)), .combine = 'rbind',.packages = c("stats","dplyr","iProPath")) %dopar% {
                                match_key_g <- x_interested_df$matching_key[ii]
                                x_interested <- dropindex_and_transpose(x_interested_df[ii,])
                                cat(sprintf("Processing gene: %s %d/%d\n", match_key_g, ii, nrow(x_interested_df)), file=sinkfile, append=TRUE)

                                x_covs_df_list <-x_df_list[setdiff(names(x_df_list),x_name)]
                                x_covs =  extract_data_g(x_covs_df_list, match_key_g)
                                x_covs_value <-dropindex_and_transpose(x_covs)

                                y_covs = extract_data_g(y_covdf_list, match_key_g)
                                y_covs_value <- dropindex_and_transpose(y_covs)
                                colnames(y_covs_value)<-paste0("y_",colnames(y_covs_value))

                                x_reg_df <- data.frame(x_interested, x_covs_value, fixed_cov, y_covs_value)
                                colnames(x_reg_df)[1] <- x_name

                                # Excluding all cis-Y sites
                                cis_Y_idx <- which(Y$matching_key==match_key_g)
                                if(length(cis_Y_idx)!=0){trans_Y <- Y[-cis_Y_idx, ]}else{trans_Y <-Y}
                                N_trans_sites <- nrow(trans_Y)
                                result_g_list <- lapply(1:nrow(trans_Y),function(s) {
                                  if (s %% 1000 == 0) {
                                    running_time<-difftime(Sys.time(),start_time, units = "mins")
                                    cat(sprintf("Processing gene %d out of %d - Sites processed: %d / %d | Elapsed time: %.3f mins\n",
                                                ii, nrow(x_interested_df),s,nrow(trans_Y),running_time), file=sinkfile, append=TRUE)
                                  }
                                  y <- dropindex_and_transpose(trans_Y[s,])
                                  reg <- summary(stats::lm(y ~ ., data = x_reg_df))

                                  param_yx <- reg$coefficients[x_name, 1]
                                  param_sd_yx <- reg$coefficients[x_name, 2]
                                  p_yx <- reg$coefficients[x_name, 4]

                                  data.frame(
                                    idx_Xgene = match_key_g, idx_Ygene = trans_Y[s,"matching_key"], idx_Ysite = trans_Y[s,"unique_index"],
                                    Pairwise_p_value = p_yx, Pairwise_Est_param = param_yx, Pairwise_Std_Error = param_sd_yx)
                                })
                                result_g <- do.call(rbind, result_g_list)
                                return(result_g)
                              }
  return(transreg_results)
}

