#' Estimate non-null proportions for iProPath
#'
#' Estimates \eqn{\pi_1^{MX}} and \eqn{\pi_1^{YM}} using \code{nonnullPropEst()},
#' then constructs a summary data frame via \code{compute_pi_df()}.
#'
#' @param p_value_MX Numeric vector of p-values for M~X associations.
#' @param p_value_YM Numeric vector of p-values for Y~M associations.
#' @param method Character string specifying the non-null proportion estimation method.
#' @param run_all Logical; passed to \code{nonnullPropEst()}.
#' @return A list with components:
#' \describe{
#'   \item{pi1_MX}{Estimated non-null proportion for M~X.}
#'   \item{pi1_YM}{Estimated non-null proportion for Y~M.}
#'   \item{pi1_MX_res}{Full result table returned by \code{nonnullPropEst()} for M~X.}
#'   \item{pi1_YM_res}{Full result table returned by \code{nonnullPropEst()} for Y~M.}
#'   \item{pi_df}{Data frame summary created by \code{compute_pi_df(pi1_MX, pi1_YM)}.}
#' }
#' @export
estimate_ipropath_pi <- function(
    p_value_MX,
    p_value_YM,
    method = "qvalue.smoother",
    run_all = FALSE
) {
  if (!is.numeric(p_value_MX)) stop("p_value_MX must be numeric.", call. = FALSE)
  if (!is.numeric(p_value_YM)) stop("p_value_YM must be numeric.", call. = FALSE)

  pi1_MX_res <- iProPath::nonnullPropEst(x = unique(p_value_MX), method = method, run_all = run_all)
  pi1_YM_res <- iProPath::nonnullPropEst(x = p_value_YM, method = method, run_all = run_all)

  # Expect nonnullPropEst returns an object indexable by [method, "pi1"]
  pi1_MX <- pi1_MX_res[method, "pi1"]
  pi1_YM <- pi1_YM_res[method, "pi1"]

  pi_df <- iProPath::compute_pi_df(pi1_MX, pi1_YM)

  list(
    pi1_MX = pi1_MX,
    pi1_YM = pi1_YM,
    pi1_MX_res = pi1_MX_res,
    pi1_YM_res = pi1_YM_res,
    pi_df = pi_df
  )
}
