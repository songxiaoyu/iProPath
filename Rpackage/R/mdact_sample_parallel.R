#' Parallel mDACT sampling-based test
#'
#' @importFrom rlang .data
#' @importFrom foreach foreach %dopar%
#' @importFrom foreach foreach %do%
#' @param p.M numeric vector of p-values from M ~ X regressions
#' @param p.M.uni unique numeric vector of p-values from M ~ X regressions
#' @param p.Y numeric vector of p-values from Y ~ M + X  regressions
#' @param pi1_MX estimated proportion of non-null signals for M ~ X p-values.
#' @param pi1_YM estimated proportion of non-null signals for Y ~ M p-values.
#' @param numCores integer, number of cores for parallel backend
#' @param y_name character, name of the outcome in Y_list
#' @return A list with the following components:
#' \describe{
#'   \item{MDACT}{Numeric vector of MDACT results.}
#'   \item{pi_df}{A data frame containing estimated \eqn{\pi} values.}
#' }
#' @export
mdact_sample_parallel <- function(p.M, p.M.uni, p.Y, pi1_MX, pi1_YM, numCores,y_name){
  pi.m.1=pi1_MX
  pi.y.1=pi1_YM
  if(pi.y.1 == 0){
    warning("Estimated non-null proportion in Y-M is zero, this likely indicates no true signal.")
    pi.y.1 = min(100/length(p.Y), 0.001)
  }
  pi.01.est=(1-pi.m.1)*pi.y.1
  pi.10.est=pi.m.1*(1-pi.y.1)
  pi.00.est=(1-pi.m.1)*(1-pi.y.1)
  # MDACT estimation
  sim.num <- length(p.M)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- max(1 - c, 0) # pi.11.est <- pi.m.1*pi.y.1
  pi_df <- data.frame(
        pi.m.1, pi.y.1, pi.01.est, pi.10.est, pi.00.est, pi.11.est
      ) %>% t() %>%
    as.data.frame() %>%
    dplyr::rename(est = .data$V1) %>%
    dplyr::mutate(
      name = c("pi1.mx", "pi1.ym", "pi01.mx_ym", "pi10.mx_ym", "pi00.mx_ym", "pi11.mx_ym"),
      formula = c(
        "M-X", "Y-M", "(1-pi1.mx)*pi1.ym", "pi1.mx*(1-pi1.ym)",
        "(1-pi1.mx)*(1-pi1.ym)", "pi1.mx*pi1.ym"
      )
    )

  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2

  F00 <- function(t){
    F_0int <- function(p.y){
      c0 <- pi.00.est * p.y^2 + pi.01.est * p.y
      c1 <- t- pi.10.est * p.y
      c2 <- (t - pi.00.est * p.y^2 - pi.10.est * p.y)/pi.01.est
      m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.y))))/(2*pi.00.est)
      s1 <- pmin(pmax(c2,0), 1)
      s2 <- pmin(pmax(m2,0), 1)
      s <- rep(0, length(p.y))
      s[c0>=c1] <- s1[c0>=c1]
      s[c0<c1] <- s2[c0<c1]
      s
    }
    x <- tryCatch(stats::integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps), error = function(e) e)
    if(!inherits(x, "error")){
      res <- stats::integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps)
    }else{
      res <- stats::integrate(F_0int, 0, 1)
    }
    res$value
  }
  if(length(p.Y)>=10000){
    n_sample = 10000
    p.Y.smp = p.Y[sample(1:length(p.Y),n_sample)]
  }else{
    p.Y.smp = p.Y
  }
  F01 <- function(t){
    c0 <- pi.00.est * p.Y.smp^2 + pi.01.est * p.Y.smp
    c1 <- t- pi.10.est * p.Y.smp
    c2 <- (t - pi.00.est * p.Y.smp^2 - pi.10.est * p.Y.smp)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y.smp))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y.smp)
  }

  F10 <- function(t){
    c0 <- pi.00.est * p.M.uni^2 + pi.10.est * p.M.uni
    c1 <- t- pi.01.est * p.M.uni
    c2 <- (t - pi.00.est * p.M.uni^2 - pi.01.est * p.M.uni)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M.uni))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)

    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M.uni)
  }
  correct <- function(t) { # equation 8. below FDR(t)_est
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) +
                        pi.01.est/(pi.01.est+pi.11.est)*F01(t) +
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) -
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
  }
  numCores <- numCores
  Ts_df=data.frame(idx=seq_len(length(Ts)),Ts)
  batch_size <- length(Ts)/(10*numCores)
  group_ids <- ceiling(seq_len(length(Ts)) / batch_size)
  batches <- split(Ts_df, group_ids)

  if (!dir.exists("Progress")) dir.create("Progress")

  sinkfile <- paste0("Progress/", y_name, "_mdact_progress.txt")
  n0 <- sink.number()
  sink(sinkfile, append = FALSE)
  on.exit({
    while (sink.number() > n0) sink()
  }, add = TRUE)

  message(sprintf(
    ">> File '%s' has been created. Check this file for progress updates.",
    sinkfile
  ))

  cat(
    sprintf(
      "======== MDACT started at %s ========\n
      Running in Parallel Mode Using %d cores\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), numCores
    ),
    file = sinkfile,
    append = TRUE
  )

  cl <- NULL
  cl <- parallel::makeCluster(numCores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, varlist = "batches", envir = environment())

  cat(
    sprintf(
      "Data Size: %d\nTotal Batches: %d\nBatch Size: %.2f\n",
      length(Ts),
      10 * numCores,
      batch_size
    ),
    file = sinkfile,
    append = TRUE
  )

  ii = NULL
  start_time <- Sys.time()
  MDACT_results <- foreach::foreach(ii = seq_along(batches), .combine = 'rbind',.packages = c("iProPath")) %dopar% {
      batch <- batches[[ii]]
      lapply(seq_along(batch$Ts), function(i) {
        ts_i=batch$Ts[[i]]
        if (i %% 1000 == 0) {
          running_time<-difftime(Sys.time(),start_time, units = "mins")
          cat(sprintf("Processing batch %d out of %d - Records processed: %d / %d | Elapsed time: %.3f mins\n",
                      ii, length(batches),i,nrow(batch),as.numeric(running_time)), file=sinkfile, append=TRUE)
        }
        x <- tryCatch(correct(ts_i), error = function(e) e)
        k <- 10
        while (inherits(x, "error") && k >= 0) {
          x <- tryCatch(correct(round(ts_i, k)), error = function(e) e)
          k <- k - 1
        }
        return(c(batch$idx[[i]], x))
      }) |> do.call(what = rbind)
  }
  res=list(MDACT=MDACT_results[,2],pi_df=pi_df)
  return(res)
}
