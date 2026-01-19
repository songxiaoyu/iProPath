#' Integrative pathway-based mediation analysis (iProPath
#'
#' @param p_value_YX numeric vector. P-values from get_Y_transRegDat.R
#' @param p_value_MX Numeric vector. P-values from get_M_RegDat.R
#' @param p_value_YM Numeric vector. P-values from get_Y_transRegDat.R
#' @param pi1_MX estimated proportion of non-null signals for \eqn{M \sim X} p-values.
#' @param pi1_YM estimated proportion of non-null signals for \eqn{Y \sim M} p-values.
#' @param idx_Xgene unique index for X gene
#' @param idx_Mediator unique index for Mediator
#' @param idx_Ygene unique index for Y gene
#' @param idx_Ysite unique index for Y site
#' @param Est_param_YX Estimated regression coefficient of \eqn{X} in the linear model \eqn{Y \sim M + X}.
#' @param Est_param_YM Estimated regression coefficient of \eqn{M} in the linear model \eqn{Y \sim M + X}.
#' @param Est_param_MX Estimated regression coefficient of \eqn{X} in the linear model \eqn{M \sim X}.
#' @param numCores integer, number of cores for parallel backend
#' @param y_name character, name of the outcome in Y_list
#' @return A data frame of row-bound results across matching keys.
#' @export
ipropath <-function(p_value_YX, p_value_MX, p_value_YM, pi1_MX, pi1_YM,
                    idx_Xgene, idx_Mediator,idx_Ygene, idx_Ysite,
                    Est_param_YX, Est_param_YM,Est_param_MX,
                    numCores=1, y_name){
  p.M.uni = p_value_MX[which(!duplicated(paste(idx_Xgene,idx_Mediator)))]
  PathM_res <- iProPath::mdact_sample_parallel(p.M=p_value_MX, p.M.uni=p.M.uni, p.Y=p_value_YM, pi1_MX, pi1_YM,
                                     numCores=numCores,y_name=y_name)
  PathM_p=PathM_res$MDACT
  pi_df=PathM_res$pi_df # Estimated non-zero proportions

  # Handling extreme p-value for Cauchy combination.
  PathM_p[PathM_p == 0] <- 1e-16
  PathM_p[PathM_p >= 1] <- 1 - 1e-16
  info = data.frame(idx_Xgene=idx_Xgene, idx_Mediator=idx_Mediator,
                    idx_Ygene=idx_Ygene, idx_Ysite=idx_Ysite,
                    Est_param_YX=Est_param_YX, Est_param_YM=Est_param_YM, Est_param_MX=Est_param_MX,
                    p_value_YX=p_value_YX)
  Path_all<-cbind(info, PathM_p)

  res_p=Path_all %>%dplyr::group_by(idx_Xgene, idx_Ysite) %>% dplyr::summarise(idx_Ygene=unique(idx_Ygene),idx_Mediator=paste(idx_Mediator,collapse = ";"),
                                                                 Path_YX_p_value = unique(p_value_YX), Path_YMX_p_value = ACAT::ACAT(Pvals =PathM_p),
                                                                 iProPath_p_value = ACAT::ACAT(Pvals = c(unique(p_value_YX),PathM_p)),
                                                                 iProPath_Est_param = sum(c(unique(Est_param_YX),Est_param_YM*Est_param_MX)))
  res=list(res_p=res_p, pi_df=pi_df)
  return(res)
}
