#' Build outcome regression data for a given X and its cis-mediators M
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom foreach foreach %do%
#' @param x_df_list named list of data frames for X
#' @param fixed_cov data frame of fixed covariates (optional, depends on your pipeline)
#' @param y_covdf_list data frame of covariates associated with Y (optional, depends on your pipeline)
#' @param mediator_list list of mediator data frames (each must include matching_key)
#' @param Y_list named list of outcome (Y) data frames.
#' @param x_name character, which element in x_df_list to use
#' @param y_name character, name of the outcome in Y_list
#' @param num_cores integer, number of cores for parallel backend
#' @return data frame (row-bind results across matching keys)
#' @export
get_Y_transRegDat <- function(x_df_list, fixed_cov, y_covdf_list, mediator_list,Y_list,x_name,y_name,num_cores) {
  if (!dir.exists("Progress")) dir.create("Progress")
  sinkfile=paste0("Progress/",y_name,"_transreg_progress.txt")
  cat(sprintf(">> File '%s' has been created. Check this file for progress updates.\n",sinkfile), file=sinkfile, append=TRUE)

  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  Y=Y_list[[y_name]]
  x_interested_df = x_df_list[[x_name]]
  # x_m_common_keys: at least having one mediator for mediation analysis
  x_m_common_keys = intersect(
    x_interested_df$matching_key,
    unique(unlist(lapply(mediator_list, function(df) df[["matching_key"]])))
  )
  x_interested_df = x_interested_df[x_interested_df$matching_key %in% x_m_common_keys,]
  ii = NULL
  start_time <- Sys.time()
  transreg_results<-foreach::foreach(ii = seq_len(nrow(x_interested_df)), .combine = 'rbind',.packages = c("stats","dplyr","iProPath")) %dopar% {
                                match_key_g <- x_interested_df$matching_key[ii]
                                x_interested <- dropindex_and_transpose(x_interested_df[ii,])
                                cat(sprintf("Processing gene with mediator: %s %d/%d\n", match_key_g, ii, nrow(x_interested_df)), file=sinkfile, append=TRUE)

                                x_covs_df_list <-x_df_list[setdiff(names(x_df_list),x_name)]
                                x_covs = extract_data_g(x_covs_df_list, match_key_g)
                                x_covs_value <- dropindex_and_transpose(x_covs)

                                m = extract_data_g(mediator_list, match_key_g)
                                m_value <- dropindex_and_transpose(
                                  m,
                                  name_style = "mediator_row_key",
                                  mediator_tag = "Mediator",
                                  sep = "_"
                                )

                                y_covs <- extract_data_g(y_covdf_list, match_key_g)
                                y_covs_value <- dropindex_and_transpose(
                                  y_covs,
                                  prefix = "y_",
                                  name_style = "prefix_col"
                                )

                                cols <- list(
                                  x_interested = x_interested,
                                  x_covs_value = x_covs_value,
                                  fixed_cov    = fixed_cov,
                                  m_value = m_value,
                                  y_covs_value = y_covs_value
                                )
                                cols <- Filter(Negate(is.null), cols)

                                x_reg_df <- do.call(
                                  data.frame,
                                  cols
                                )
                                colnames(x_reg_df)[1] = x_name

                                # Excluding all cis-Y sites
                                cis_Y_idx <- which(Y$matching_key==match_key_g)
                                if(length(cis_Y_idx)!=0) {trans_Y <- Y[-cis_Y_idx, ]}else{trans_Y <-Y}
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

                                  reg_coef_names = rownames(reg$coefficients)
                                  m_name <- grep("^Mediator", colnames(x_reg_df), value = TRUE)
                                  m_name <- intersect(reg_coef_names, m_name)

                                  param_ym <- reg$coefficients[m_name, 1]
                                  param_sd_ym <- reg$coefficients[m_name, 2]
                                  p_ym <- reg$coefficients[m_name, 4]

                                  l <- length(m_name) # number of Mediators in this gene dataframe.
                                  data.frame(
                                    idx_Xgene = rep(match_key_g, l),
                                    idx_Mediator = m_name,
                                    idx_Ygene = rep(trans_Y[s,"matching_key"], l),
                                    idx_Ysite = rep(trans_Y[s,"unique_index"], l),

                                    p_value_YX = rep(p_yx, l),
                                    Est_param_YX = rep(param_yx, l),
                                    Std_Error_YX = rep(param_sd_yx, l),

                                    p_value_YM = p_ym,
                                    Est_param_YM = param_ym,
                                    Std_Error_YM = param_sd_ym
                                  )
                                })
                                result_g <- do.call(rbind, result_g_list)
                                return(result_g)
                              }
  return(transreg_results)
}
