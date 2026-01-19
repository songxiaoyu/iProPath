#' Compute pi dataframe
#'
#' @param pi1_MX estimated proportion of non-null signals for M ~ X p-values.
#' @param pi1_YM estimated proportion of non-null signals for Y ~ M p-values.
#' @return data frame
#' @export
compute_pi_df <- function(pi1_MX, pi1_YM) {

  # ---- input checks ----
  if (!is.numeric(pi1_MX) || length(pi1_MX) != 1 || is.na(pi1_MX)) {
    stop("pi1_MX must be a non-missing numeric scalar.")
  }
  if (!is.numeric(pi1_YM) || length(pi1_YM) != 1 || is.na(pi1_YM)) {
    stop("pi1_YM must be a non-missing numeric scalar.")
  }

  pi.m.1 <- pi1_MX
  pi.y.1 <- pi1_YM

  pi.01.est <- (1 - pi.m.1) * pi.y.1
  pi.10.est <- pi.m.1 * (1 - pi.y.1)
  pi.00.est <- (1 - pi.m.1) * (1 - pi.y.1)
  pi.11.est <- pi.m.1 * pi.y.1

  df <- data.frame(
    pi.m.1 = pi.m.1,
    pi.y.1 = pi.y.1,
    pi.01.est = pi.01.est,
    pi.10.est = pi.10.est,
    pi.00.est = pi.00.est,
    pi.11.est = pi.11.est
  )

  df <- df |>
    dplyr::rename(
      pi1.mx    = pi.m.1,
      pi1.ym    = pi.y.1,
      pi01.mx_ym = pi.01.est,
      pi10.mx_ym = pi.10.est,
      pi00.mx_ym = pi.00.est
    ) |>
    base::t() |>
    base::as.data.frame() |>
    dplyr::rename(est = .data$V1) |>
    dplyr::mutate(
      formula = c(
        "M-X signal prop",
        "Y-M signal prop",
        "(1-pi1.mx)pi1.ym",
        "pi1.mx(1-pi1.ym)",
        "(1-pi1.mx)(1-pi1.ym)",
        "pi1.mx * pi1.ym"
      )
    )

  return(df)
}
