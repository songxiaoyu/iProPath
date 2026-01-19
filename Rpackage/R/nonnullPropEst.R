#' Estimate non-null proportion
#'
#' @param x Numeric vector of test statistics or p-values used to estimate the non-null proportion.
#'
#' @param method Character string specifying the estimation method. Options include:
#' \describe{
#'   \item{"qvalue.smoother"}{Storey--Tibshirani q-value method with smoothing.}
#'   \item{"qvalue.bootstrap"}{Storey--Tibshirani q-value method with bootstrap tuning.}
#'   \item{"JC"}{Jin--Cai estimator based on the empirical characteristic function.}
#'   \item{"locfdr.central_matching"}{Local FDR estimator using central matching.}
#'   \item{"locfdr.mle"}{Local FDR estimator using maximum likelihood estimation.}
#' }
#'
#' @param run_all Logical; if \code{TRUE}, estimates are computed for all available methods.
#' If \code{FALSE}, only the method specified by \code{method} is used.
#' @return data frame
#' @export
nonnullPropEst<-function(x,method=c("qvalue.smoother","qvalue.bootstrap","JC",
                                    "locfdr.central_matching","locfdr.mle"),
                         run_all = TRUE){
  choices = c("qvalue.smoother","qvalue.bootstrap","JC","locfdr.central_matching","locfdr.mle")
  if (run_all) {
    method <- choices
  }
  method <- match.arg(method, choices = choices, several.ok = TRUE)
  results <- list()
  if ("qvalue.smoother" %in% method) {
    results$qvalue.smoother= iProPath::nonnullPropEst_q(x,pi0.method="smoother")
  }
  if ("qvalue.bootstrap" %in% method) {
    results$qvalue.bootstrap= iProPath::nonnullPropEst_q(x,pi0.method="bootstrap")
  }
  if ("JC" %in% method) {
    results$JC=iProPath::nonnullPropEst_JC(x,u=0,sigma=1, parallel=F)
  }
  if ("locfdr.central_matching" %in% method) {
    results$locfdr.central_matching= iProPath::nonnullPropEst_locfdr(x,est="central matching",u=0,sigma=1)
  }
  if ("locfdr.mle" %in% method) {
    results$locfdr.mle= iProPath::nonnullPropEst_locfdr(x,est="mle",u=0,sigma=1)
  }

  df_results <- data.frame(method = names(results), pi1 = unlist(results))

  return(df_results)
}

#' Estimate non-null proportion (q-value based)
#'
#' @param x Numeric vector of p-values.
#' @param pi0.method Character string specifying the \eqn{\pi_0} estimation method;
#'   either \code{"smoother"} or \code{"bootstrap"}.
#' @return Numeric scalar giving the estimated non-null proportion.
#' @export
nonnullPropEst_q <- function(x,pi0.method=c("smoother","bootstrap")){
  l2 = max(x)
  pi0_est=qvalue::pi0est(x, lambda = seq(0.05, l2, 0.05), pi0.method =pi0.method,
                         smooth.df = 3, smooth.log.pi0 = FALSE)
  pi1_est=1-pi0_est$pi0
  return(pi1_est)
}

#' Estimate non-null proportion (JC method)
#'
#' @param p Numeric vector of p-values.
#' @param u Numeric; tuning parameter (default 0).
#' @param sigma Numeric; noise scale parameter (default 1).
#' @param parallel Logical; whether to run computation in parallel.
#'
#' @return epsest Estimated non-null proportion.
#' @export
nonnullPropEst_JC <- function(p,u=0,sigma=1, parallel=F)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  x = stats::qnorm(p,lower.tail = F) # ??
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)

  epsest_fun=function(j) {
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    if (parallel==F) {
      co=mapply(function(f, f2)  mean(cos(t*xi[f]*z)), 1:101, 1) #f2?
    } else {
      co=parallel::mcmapply(function(f, f2)  mean(cos(t*xi[f]*z)), 1:101, 1)
    }
    epshat = 1 - sum(w*f*co)/sum(w)

    return(epshat)
  }
  if (parallel==F) {
    epsest=sapply(1:length(tt), function(f) epsest_fun(f))
  } else {
    epsest=parallel::mcmapply(function(f,f2) epsest_fun(f), 1:length(tt), 1)
  }
  return(epsest=max(epsest))
}

#' Estimate non-null proportion (local FDR)
#'
#' @param p Numeric vector of p-values.
#' @param est Character string specifying the locfdr estimation method;
#'   either \code{"central matching"} or \code{"mle"}.
#' @param u Numeric value specifying the mean of the null distribution.
#' @param sigma Numeric value specifying the standard deviation of the null distribution.
#' @return Numeric scalar giving the estimated non-null proportion.
#' @export
nonnullPropEst_locfdr <- function(p,est=c("central matching","mle"),u=0,sigma=1)
{
  x = stats::qnorm(p, lower.tail = FALSE)
  z = (x - u)/sigma
  if(est=="central matching"){
    w <- locfdr::locfdr(z, plot=1,nulltype=1)$fp0["cmest","p0"]
  }
  if(est=="mle"){
    w <- locfdr::locfdr(z, plot=1,nulltype=1)$fp0["mlest","p0"]
  }
  pi1= 1-w
  return(pi1)
}


