
#--------------------------------------------------------------------
# Objective: online changepoint detection for one ewe using a marked Poisson point process and Monte Carlo calibration
#
# Main functions:
#   sim_changepointPP             Simulate marked Poisson events
#   compute_G_one_day             Compute the daily detection statistic G(d)
#   bootstrap_G_null              Simulate the null distribution of G(d)
#   calibrate_min_p_null          Calibrate the sequential p-value threshold
#   compute_pvalue_one_day        Compute the p-value for one ewe on one day
#   detect_rupture_one_ewe_online Run th e complete online detection procedure
#--------------------------------------------------------------------

#' Simulate a marked Poisson process with one changepoint

#'
#' @param tau0 Numeric. Changepoint time.
#' @param lambda Numeric vector of length 2. `lambda[1]` is the event rate
#'   before `tau0`, and `lambda[2]` is the event rate after `tau0`.
#' @param mu Numeric vector of length 2. `mu[1]` is the exponential mark rate
#'   before `tau0`, and `mu[2]` is the exponential mark rate after `tau0`.
#' @param t_final Numeric. Final observation time.
#'
#' @return A data frame sorted by event time with columns:
#'   `time`, `marks`, and `state`. State 0 denotes events before the
#'   changepoint and state 1 denotes events after the changepoint.
#'

sim_changepointPP <- function(tau0, lambda, mu, t_final) {
  
  N_minus <- rpois(n = 1, lambda = lambda[1] * tau0)
  N_plus  <- rpois(n = 1, lambda = lambda[2] * (t_final - tau0))
  
  time_minus <- sort(runif(n = N_minus, min = 0, max = tau0))
  time_plus  <- sort(runif(n = N_plus, min = tau0, max = t_final))
  
  marks_minus <- rexp(N_minus, rate = mu[1])
  marks_plus  <- rexp(N_plus, rate = mu[2])
  
  out <- data.frame(
    time  = c(time_minus, time_plus),
    marks = c(marks_minus, marks_plus),
    state = as.factor(c(rep(0, N_minus), rep(1, N_plus)))
  )
  
  out <- out[order(out$time), ]
  rownames(out) <- NULL
  
  return(out)
}

#--------------------------------------------------------------------
#' Compute the online detection statistic for one ewe on one day
#
#' @param ewe_ind Data frame for one ewe containing at least `time` and `marks`.
#' @param d Integer or numeric. Current monitoring day.
#' @param lambda0 Numeric. Event rate under the pre-change model.
#' @param mu0 Numeric. Mark rate under the pre-change model.
#' @param lambda1 Numeric. Event rate under the post-change model.
#' @param mu1 Numeric. Mark rate under the post-change model.
#' @param L_min Numeric. Minimum allowed duration between a candidate
#'   changepoint and the current day.
#'
#' @return A list with:
#'   \itemize{
#'     \item `G`: the maximum  log-likelihood-ratio statistic;
#'     \item `S_hat`: the maximizing changepoint candidate.
#'   }

compute_G_one_day <- function(ewe_ind, d,lambda0, mu0,lambda1, mu1,L_min = 0) {
  
  ewe_ind <- ewe_ind[ewe_ind$time < d, ]
  ewe_ind <- ewe_ind[order(ewe_ind$time), ]
  
  times <- ewe_ind$time
  marks <- ewe_ind$marks
  
  log_ratio <- log((lambda1 * mu1) / (lambda0 * mu0))
  dlambda   <- lambda1 - lambda0
  dmu       <- mu1 - mu0
  
  S_grid <- c(0, times[times < d])
  S_grid <- S_grid[(d - S_grid) >= L_min]
  
  if (length(S_grid) == 0) {
    return(list(G = 0, S_hat = NA_real_))
  }
  
  L_vals <- numeric(length(S_grid))
  
  for (k in seq_along(S_grid)) {
    S <- S_grid[k]
    idx <- which(times > S & times < d)
    N_sd <- length(idx)
    V_sd <- if (N_sd == 0) 0 else sum(marks[idx])
    L_vals[k] <- N_sd * log_ratio - dlambda * (d - S) -dmu * V_sd
  }
  
  k_best <- which.max(L_vals)
  L_best <- L_vals[k_best]
  
  G <-  L_best
  S_hat <-  S_grid[k_best] 
  
  return(list(G = G,S_hat = S_hat ))
}

#' Simulate the null distribution of the daily statistic
#'
#' @param B Integer. Number of Monte Carlo trajectories.
#' @param K Integer. Number of monitoring days and final simulation time.
#' @param lambda0 Numeric. Event rate under the null model.
#' @param mu0 Numeric. Mark rate under the null model.
#' @param lambda1 Numeric. Alternative event rate used in the likelihood ratio.
#' @param mu1 Numeric. Alternative mark rate used in the likelihood ratio.
#' @param L_min Numeric. Minimum allowed interval length after a candidate changepoint.
#' @param seed Optional integer random seed.
#'
#' @return  `G(d)` null distribution of the daily statistic.
bootstrap_G_null <- function(B, K, lambda0, mu0, lambda1, mu1, L_min = 0, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  G_null <- matrix(0, nrow = B, ncol = K)
  colnames(G_null) <- paste0("day_", seq_len(K))
  
  for (b in seq_len(B)) {
    dat_b <- sim_changepointPP( tau0  = K,  lambda  = c(lambda0, lambda0), mu = c(mu0, mu0), t_final = K)
    for (d in seq_len(K)) {
      out_d <- compute_G_one_day(dat_b,d,lambda0, mu0,lambda1,mu1,L_min )
      G_null[b, d] <- out_d$G
    }
  }
  return(G_null)
}

#--------------------------------------------------------------------
#' Calibrate the sequential minimum-p-value threshold

#' @param G_null Numeric matrix returned by `bootstrap_G_null()`.
#' @param d_min Integer. First monitoring day included in calibration.
#' @param alpha Numeric in `[0,1]`, false-alarm level.
#'
#' @return A list with:
#'   \itemize{
#'     \item `P_null`: matrix of empirical day-specific p-values;
#'     \item `R_null`: minimum p-value for each null trajectory;
#'     \item `u_alpha`: calibrated online decision threshold. }
#'     
calibrate_min_p_null <- function(G_null, d_min = 1, alpha = 0.05) {
  
  B <- nrow(G_null)
  K <- ncol(G_null)
  P_null <- matrix(NA_real_, nrow = B, ncol = K)
  colnames(P_null) <- colnames(G_null)
  
  for (d in seq_len(K)) {
    G_day <- G_null[, d]
    P_null[, d] <- sapply(G_day, function(x) { (sum(G_day >= x)) / (B ) })
  }
  
  R_null <- apply( P_null[, d_min:K, drop = FALSE],  1, min,  na.rm = TRUE )
  
  # Online threshold: alarm at first day p_d <= u_alpha
  u_alpha <- as.numeric(quantile(R_null, probs = alpha, type = 1))
  
  return(list( P_null = P_null, R_null = R_null, u_alpha = u_alpha
  ))
}


#--------------------------------------------------------------------
#' Compute the empirical p-value for one ewe on one day
#'
#' @param ewe_ind Data frame for one ewe containing at least `time` and `marks`.
#' @param d Integer. Monitoring day.
#' @param G_null Numeric null-statistic matrix with one column per day.
#' @param lambda0 Numeric. Event rate under the pre-change model.
#' @param mu0 Numeric. Mark rate under the pre-change model.
#' @param lambda1 Numeric. Event rate under the post-change model.
#' @param mu1 Numeric. Mark rate under the post-change model.
#' @param L_min Numeric. Minimum allowed interval length after a candidate changepoint.
#'
#' @return A list containing the day, observed statistic, empirical p-value, and estimated changepoint.
#' 
compute_pvalue_one_day <- function(ewe_ind, d, G_null, lambda0, mu0,lambda1, mu1,L_min = 0) {
  
  B <- nrow(G_null)
  G_out <- compute_G_one_day( ewe_ind, d,lambda0,mu0, lambda1,mu1,L_min)
  
  p_d <- (sum(G_null[, d] >= G_out$G)) / B

  return(list(day = d,G_obs = G_out$G, p_value = p_d, S_hat = G_out$S_hat ))
}

#--------------------------------------------------------------------
# Final function
#' Detect a changepoint online for one ewe
#'
#' @param ewe_ind Data frame for one ewe containing at least `time` and `marks`.
#' @param K Integer. Final monitoring day.
#' @param lambda0 Numeric. Event rate under the pre-change model.
#' @param mu0 Numeric. Mark rate under the pre-change model.
#' @param lambda1 Numeric. Event rate under the post-change model.
#' @param mu1 Numeric. Mark rate under the post-change model.
#' @param B Integer. Number of Monte Carlo trajectories when `G_null` is not
#'   supplied.
#' @param alpha Numeric in `[0,1]`. Target sequential false-alarm level.
#' @param d_min Integer. First day on which an alarm is allowed.
#' @param L_min Numeric. Minimum allowed interval length after a candidate
#'   changepoint.
#' @param seed Optional integer random seed.
#' @param G_null Optional precomputed null-statistic matrix.
#' @param verbose Logical. Print the detection summary when `TRUE`.
#' @param return_null Logical. Include null-calibration objects in the output
#'   when `TRUE`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `alarm`: whether a changepoint alarm occurred;
#'     \item `tau`: first alarm day, or `Inf` if no alarm occurred;
#'     \item `S_hat`: estimated changepoint at the alarm day;
#'     \item `u_alpha`: calibrated sequential p-value threshold;
#'     \item `alpha`: requested false-alarm level;
#'     \item `path`: daily statistics, p-values, changepoint estimates, and
#'       alarm indicators.
#'   }

detect_rupture_one_ewe_online <- function(ewe_ind, K,lambda0, mu0,lambda1, mu1,B = 1000,alpha = 0.05,d_min = 1,
                                          L_min = 0,seed = NULL,G_null = NULL, verbose = TRUE,return_null = TRUE) {
  

  #  Bootstrap under H0 
  if (is.null(G_null)) {
    G_null <- bootstrap_G_null(B, K,lambda0,mu0, lambda1, mu1, L_min,seed = seed )}
  
  B0 <- nrow(G_null)
  
  #  Calibrate min-p threshold
  calib <- calibrate_min_p_null(G_null,d_min,alpha )
  
  u_alpha <- calib$u_alpha
  
  #  Online monitoring
  path <- data.frame(day = seq_len(K), G_obs = NA_real_, p_value = NA_real_, S_hat = NA_real_, alarm_here = FALSE )
  
  alarm <- FALSE
  tau <- Inf
  S_alarm <- NA_real_
  
  for (d in seq_len(K)) {
    G_out <- compute_G_one_day( ewe_ind, d,lambda0,mu0, lambda1, mu1, L_min)
    p_d <- (sum(G_null[, d] >= G_out$G)) / (B0 )
    path$G_obs[d] <- G_out$G
    path$p_value[d] <- p_d
    path$S_hat[d] <- G_out$S_hat
    if (d >= d_min && p_d <= u_alpha) {
      alarm <- TRUE
      tau <- d
      S_alarm <- G_out$S_hat
      path$alarm_here[d] <- TRUE
      break
    }
  }
  
  if (verbose) {
    if (alarm) {
      cat("ALARM: YES\n")
      cat("Alarm day tau =", tau, "\n")
      cat("Estimated rupture time S_hat(tau) =", S_alarm, "\n")
      cat("p-value at alarm =", path$p_value[tau], "\n")
      cat("G at alarm =", path$G_obs[tau], "\n")
    } else {
      cat("ALARM: NO before T =", K, "\n")
      cat("Minimum observed p-value =", min(path$p_value[d_min:K], na.rm = TRUE), "\n")
    }
  }
  
  out <- list(alarm,tau, S_alarm, u_alpha, alpha,path)
  
  if (return_null) {
    out$G_null <- G_null
    out$R_null <- calib$R_null
    out$P_null <- calib$P_null
  }
  
  return(out)
}

