source("Online.R")


#--------------------------------------------------------------------
# Objective: online changepoint detection for a herd of ewes using a marked Poisson point process with ram-specific parameters
#            and whole-herd family-wise error rate (FWER) control
#
# Main functions:
#   sim_changepointPP_ram               Simulate one ewe with ram effects
#   sim_herd_changepoint_ram            Simulate a complete herd
#   compute_G_one_day_ram               Compute G(d) for one ewe
#   bootstrap_G_null_herd_ram           Simulate herd-level null statistics
#   calibrate_min_p_null_herd           Calibrate the whole-herd FWER threshold
#   compute_pvalues_herd_one_day_ram    Compute daily p-values for active ewes
#   detect_rupture_herd_online_FWER_ram Run the complete online procedure
#--------------------------------------------------------------------'

#' Determine the ram labels
get_ram_levels <- function(lambda0, ram_levels = NULL) {
  if (!is.null(ram_levels)) return(ram_levels)
  if (!is.null(names(lambda0)) && all(names(lambda0) != "")) return(names(lambda0))
  return(seq_along(lambda0))
}
#--------------------------------------------------------------------'

#' Simulate one marked Poisson process with a changepoint and ram effect
#' @param tau0 Numeric. Changepoint time.
#' @param lambda0 Numeric vector. Ram-specific event rates before `tau0`.
#' @param lambda1 Numeric vector. Ram-specific event rates after `tau0`.
#' @param mu0 Numeric vector. Ram-specific exponential mark rates before `tau0`.
#' @param mu1 Numeric vector. Ram-specific exponential mark rates after `tau0`.
#' @param t_final Numeric. Final observation time.
#' @param ram_levels Optional vector of ram labels. When `NULL`, labels are obtained through `get_ram_levels(lambda0, ram_levels)`.
#'
#' @return A data frame contient:time; mark; ram associated ; factor: 0 before the changepoint and 1 after it.

sim_changepointPP_ram <- function(tau0, lambda0, lambda1, mu0, mu1,
                                  t_final, ram_levels = NULL) {
  
  ram_levels <- get_ram_levels(lambda0, ram_levels)
  R <- length(ram_levels)
  
  out_list <- list()
  counter <- 1
  
  for (rr in seq_len(R)) {
    
    ram_r <- ram_levels[rr]
    
    N_minus <- rpois(1, lambda0[rr] * tau0)
    N_plus  <- rpois(1, lambda1[rr] * (t_final - tau0))
    
    time_minus <- runif(N_minus, min = 0, max = tau0)
    time_plus  <- runif(N_plus, min = tau0, max = t_final)
    
    marks_minus <- rexp(N_minus, rate = mu0[rr])
    marks_plus  <- rexp(N_plus, rate = mu1[rr])
    
    dat_r <- data.frame(
      time  = c(time_minus, time_plus),
      marks = c(marks_minus, marks_plus),
      ram   = c(rep(ram_r, N_minus), rep(ram_r, N_plus)),
      state = as.factor(c(rep(0, N_minus), rep(1, N_plus)))
    )
    
    if (nrow(dat_r) > 0) {
      out_list[[counter]] <- dat_r
      counter <- counter + 1
    }
  }
  
  if (length(out_list) == 0) {
    return(data.frame(
      time = numeric(0),
      marks = numeric(0),
      ram = character(0),
      state = factor()
    ))
  }
  
  out <- do.call(rbind, out_list)
  out <- out[order(out$time), ]
  rownames(out) <- NULL
  
  return(out)
}

#--------------------------------------------------------------------
#' Simulate a herd with ewe-specific changepoints and ram effects
#'
#' @param tau_vec Numeric vector. One changepoint time per ewe.
#' @param lambda0 Numeric vector. Ram-specific event rates before a changepoint.
#' @param lambda1 Numeric vector. Ram-specific event rates after a changepoint.
#' @param mu0 Numeric vector. Ram-specific mark rates before a changepoint.
#' @param mu1 Numeric vector. Ram-specific mark rates after a changepoint.
#' @param t_final Numeric. Final observation time.
#' @param ram_levels Optional vector of ram labels.
#'
#' @return A data frame containing simulated events for all ewes. 
sim_herd_changepoint_ram <- function(tau_vec, lambda0, lambda1, mu0, mu1,t_final, ram_levels = NULL) {
  
  n_ewe <- length(tau_vec)
  out_list <- vector("list", n_ewe)
  
  for (i in seq_len(n_ewe)) {
    
    dat_i <- sim_changepointPP_ram(tau_vec[i],lambda0,lambda1,mu0, mu1,t_final, ram_levels )
    dat_i$id <- i
    out_list[[i]] <- dat_i
  }
  
  out <- do.call(rbind, out_list)
  out <- out[order(out$id, out$time), ]
  rownames(out) <- NULL
  
  return(out)
}

#--------------------------------------------------------------------
#' Compute the daily GLR statistic for one ewe with ram effects
#'
#' @param ewe_ind Data frame for one ewe containing at least `time`, `marks`,
#'   and `ram`.
#' @param d Integer or numeric. Current monitoring day.
#' @param lambda0 Numeric vector. Ram-specific pre-change event rates.
#' @param mu0 Numeric vector. Ram-specific pre-change mark rates.
#' @param lambda1 Numeric vector. Ram-specific post-change event rates.
#' @param mu1 Numeric vector. Ram-specific post-change mark rates.
#' @param L_min Numeric. Minimum allowed duration between a candidate
#'   changepoint and the current day.
#' @param ram_levels Optional vector defining ram labels and their order.
#'
#' @return A list with:
#'   \itemize{
#'     \item `G`: maximum log-likelihood-ratio statistic;
#'     \item `S_hat`: maximizing changepoint candidate`.
#'   }
#'
compute_G_one_day_ram <- function(ewe_ind, d,lambda0, mu0, lambda1, mu1, L_min = 0, ram_levels = NULL) {
  
  ram_levels <- get_ram_levels(lambda0, ram_levels)
  R <- length(ram_levels)
  
  ewe_ind <- ewe_ind[ewe_ind$time < d, ]
  ewe_ind <- ewe_ind[order(ewe_ind$time), ]
  
  times <- ewe_ind$time
  marks <- ewe_ind$marks
  rams  <- ewe_ind$ram
  
  log_ratio <- log((lambda1 * mu1) / (lambda0 * mu0))
  dlambda   <- lambda1 - lambda0
  dmu       <- mu1 - mu0
  
  S_grid <- c(0, times[times < d])
  S_grid <- unique(S_grid)
  S_grid <- S_grid[(d - S_grid) >= L_min]
  
  if (length(S_grid) == 0) {
    return(list(G = 0, S_hat = NA_real_))
  }
  
  L_vals <- numeric(length(S_grid))
  
  for (k in seq_along(S_grid)) {
    
    S <- S_grid[k]
    idx_post <- which(times > S & times < d)
    val <- 0
    
    for (rr in seq_len(R)) {
      
      ram_r <- ram_levels[rr]
      idx_r <- idx_post[as.character(rams[idx_post]) == as.character(ram_r)]
      
      N_r <- length(idx_r)
      V_r <- if (N_r == 0) 0 else sum(marks[idx_r])
      
      val <- val +N_r * log_ratio[rr] -dlambda[rr] * (d - S) -  dmu[rr] * V_r
    }
    
    L_vals[k] <- val
  }
  
  k_best <- which.max(L_vals)
  L_best <- L_vals[k_best]
  
  G <- L_best
  S_hat <-  S_grid[k_best] 
  
  return(list(  G = G, S_hat = S_hat))
}

#--------------------------------------------------------------------
#' Bootstrap the null GLR statistic for the whole herd
#'
#' @param B Integer. Number of Monte Carlo replications.
#' @param n_ewe Integer. Number of ewes in each simulated herd.
#' @param K Integer. Number of monitoring days and final simulation time.
#' @param lambda0 Numeric vector. Ram-specific null event rates.
#' @param mu0 Numeric vector. Ram-specific null mark rates.
#' @param lambda1 Numeric vector. Alternative event rates used in the GLR.
#' @param mu1 Numeric vector. Alternative mark rates used in the GLR.
#' @param L_min Numeric. Minimum allowed interval after a candidate changepoint.
#' @param seed Optional integer random seed.
#' @param ram_levels Optional vector defining ram labels and their order.
#'
#' @return A numeric array with dimensions `B x n_ewe x K`. Element
#'   `[b, i, d]` is the null statistic for replication `b`, ewe `i`, and day  `d`.
#'
bootstrap_G_null_herd_ram <- function(B, n_ewe, K, lambda0, mu0, lambda1, mu1, L_min = 0,  seed = NULL,ram_levels = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ram_levels <- get_ram_levels(lambda0, ram_levels)
  
  G_null <- array( 0,
    dim = c(B, n_ewe, K),
    dimnames = list(
      paste0("boot_", seq_len(B)),
      paste0("id_", seq_len(n_ewe)),
      paste0("day_", seq_len(K))
    )
  )
  
  for (b in seq_len(B)) {
    
    for (i in seq_len(n_ewe)) {

      ### Under H0: no changepoint, so lambda1=lambda0 and mu1=mu0
      dat_bi <- sim_changepointPP_ram( K,lambda0,lambda0,mu0,mu0,K,ram_levels = ram_levels  )
      
      for (d in seq_len(K)) {
        
        out_d <- compute_G_one_day_ram(dat_bi, d,lambda0,mu0,lambda1, mu1,L_min,ram_levels)
        G_null[b, i, d] <- out_d$G
      }
    }
  }
  
  return(G_null)
}


#--------------------------------------------------------------------
#' Calibrate the whole-herd online FWER threshold
#'
#' For bootstrap replication `b`, `M_null[b]` is the minimum p-value across all
#' ewes and monitored days from `d_min` through `K`. The lower `alpha` quantile
#' of `M_null` defines the common online threshold `u_alpha`.
#'
#' @param G_null Numeric array with dimensions `B x n_ewe x K`.
#' @param d_min Integer. First day included in sequential calibration.
#' @param alpha Numeric in `[0,1]`. Target whole-herd FWER level.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `P_null`: null p-value array with the same dimensions as `G_null`;
#'     \item `M_null`: minimum p-value for each bootstrap replication;
#'     \item `u_alpha`: calibrated whole-herd alarm threshold.
#'   }
calibrate_min_p_null_herd <- function(G_null, d_min = 1, alpha = 0.05) {
  
  B     <- dim(G_null)[1]
  n_ewe <- dim(G_null)[2]
  K     <- dim(G_null)[3]

  P_null <- array( NA_real_, dim = c(B, n_ewe, K),dimnames = dimnames(G_null))
  
  for (d in seq_len(K)) {
    
    G_day <- as.vector(G_null[, , d])
    Npool <- length(G_day)
    
    p_day <- sapply(G_day, function(x) {
      sum(G_day >= x) / Npool
    })
    
    P_null[, , d] <- array(p_day, dim = c(B, n_ewe))
  }
  
  M_null <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    M_null[b] <- min(P_null[b, , d_min:K, drop = FALSE], na.rm = TRUE)
    }
  cand <- sort(unique(M_null))
  u_alpha <- max(c(0, cand[ecdf(M_null)(cand) <= alpha]))
  
  return(list(P_null = P_null, M_null = M_null,  u_alpha = u_alpha))
}

#--------------------------------------------------------------------
#' Compute daily p-values for active ewes in the herd
#'
#' @param dat_herd Data frame containing herd observations with columns `id`,
#'   `time`, `marks`, and `ram`.
#' @param ids_active Vector of ewe identifiers still under monitoring.
#' @param d Integer. Current monitoring day.
#' @param G_null Numeric null-statistic array with dimensions `B x n_ewe x K`.
#' @param lambda0 Numeric vector. Ram-specific pre-change event rates.
#' @param mu0 Numeric vector. Ram-specific pre-change mark rates.
#' @param lambda1 Numeric vector. Ram-specific post-change event rates.
#' @param mu1 Numeric vector. Ram-specific post-change mark rates.
#' @param L_min Numeric. Minimum allowed interval after a candidate changepoint.
#' @param ram_levels Optional vector defining ram labels and their order.
#'
#' @return A data frame with one row per active ewe and columns `id`, `day`,`G_obs`, `p_value`, and `S_hat`.
#'
compute_pvalues_herd_one_day_ram <- function(dat_herd, ids_active, d,  G_null, lambda0, mu0, lambda1, mu1,
                                             L_min = 0,     ram_levels = NULL) {
  
  G_day_null <- as.vector(G_null[, , d])
  Npool <- length(G_day_null)
  
  out <- data.frame(id = ids_active,day = d, G_obs = NA_real_, p_value = NA_real_,S_hat = NA_real_ )
  
  for (a in seq_along(ids_active)) {
    
    id_i <- ids_active[a]
    dat_i <- dat_herd[dat_herd$id == id_i, ]
    
    G_out <- compute_G_one_day_ram(dat_i,d,lambda0, mu0,lambda1, mu1,L_min, ram_levels )
    
    p_i_d <- sum(G_day_null >= G_out$G) / Npool
    
    out$G_obs[a] <- G_out$G
    out$p_value[a] <- p_i_d
    out$S_hat[a] <- G_out$S_hat
  }
  
  return(out)
}


#--------------------------------------------------------------------
#' Run online herd detection with whole-herd FWER control and ram effects
#'
#' @param dat_herd Data frame containing herd observations with columns `id`,
#'   `time`, `marks`, and `ram`.
#' @param ids Vector of ewe identifiers to monitor.
#' @param K Integer. Final monitoring day.
#' @param lambda0 Numeric vector. Ram-specific pre-change event rates.
#' @param mu0 Numeric vector. Ram-specific pre-change mark rates.
#' @param lambda1 Numeric vector. Ram-specific post-change event rates.
#' @param mu1 Numeric vector. Ram-specific post-change mark rates.
#' @param B Integer. Number of Monte Carlo replications when `G_null` is absent.
#' @param alpha Numeric in `[0,1]`. Target whole-herd FWER level.
#' @param d_min Integer. First day on which an alarm is allowed.
#' @param L_min Numeric. Minimum allowed interval after a candidate changepoint.
#' @param seed Optional integer random seed.
#' @param G_null Optional precomputed null-statistic array.
#' @param verbose Logical. Print a summary when `TRUE`.
#' @param return_null Logical. Include null-calibration objects in the output.
#' @param ram_levels Optional vector defining ram labels and their order.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `result`: one-row-per-ewe alarm summary;
#'     \item `path`: daily monitoring results for ewes while active;
#'     \item `u_alpha`: calibrated whole-herd threshold;
#'     \item `alpha`: requested FWER level.
#'   }
#'   When `return_null = TRUE`, the output also contains `G_null`, `M_null`,
#'   and `P_null`.
#'
detect_rupture_herd_online_FWER_ram <- function(dat_herd,ids, K,lambda0, mu0, lambda1, mu1, 
                                           B = 1000, alpha = 0.05,d_min = 1,   L_min = 0,  seed = NULL,
                                          G_null = NULL,  verbose = TRUE, return_null = TRUE,ram_levels = NULL) {
  
  n_ewe <- length(ids)
  ram_levels <- get_ram_levels(lambda0, ram_levels)
  
  if (is.null(G_null)) {
    G_null <- bootstrap_G_null_herd_ram(B,n_ewe,K, lambda0,mu0,lambda1, mu1,L_min, seed, ram_levels)
  }
  
  calib <- calibrate_min_p_null_herd(G_null = G_null, d_min = d_min,alpha = alpha)
  u_alpha <- calib$u_alpha
  result <- data.frame(
    id = ids,
    alarm = FALSE,
    tau = Inf,
    S_hat = NA_real_,
    p_at_alarm = NA_real_,
    G_at_alarm = NA_real_
  )
  
  active_ids <- ids
  path_list <- list()
  counter <- 1
  
  for (d in seq_len(K)) {
    
    if (length(active_ids) == 0) break
    
    day_out <- compute_pvalues_herd_one_day_ram( dat_herd,active_ids,d, G_null, lambda0, 
                                                 mu0,lambda1, mu1, L_min,ram_levels = ram_levels )
    
    day_out$alarm_here <- FALSE
    
    if (d >= d_min) {
      
      alarm_ids <- day_out$id[day_out$p_value <= u_alpha]
      
      if (length(alarm_ids) > 0) {
        
        for (id_alarm in alarm_ids) {
          
          row_day <- day_out[day_out$id == id_alarm, ]
          idx_res <- which(result$id == id_alarm)
          
          result$alarm[idx_res] <- TRUE
          result$tau[idx_res] <- d
          result$S_hat[idx_res] <- row_day$S_hat
          result$p_at_alarm[idx_res] <- row_day$p_value
          result$G_at_alarm[idx_res] <- row_day$G_obs
          
          day_out$alarm_here[day_out$id == id_alarm] <- TRUE
        }
        
        active_ids <- setdiff(active_ids, alarm_ids)
      }
    }
    
    path_list[[counter]] <- day_out
    counter <- counter + 1
  }
  
  path <- do.call(rbind, path_list)
  rownames(path) <- NULL
  
  if (verbose) {
   
    if (sum(result$alarm) > 0) {
      print(result[result$alarm, ])
    }
  }
  
  out <- list( result = result, path = path,u_alpha = u_alpha, alpha = alpha )
  
  if (return_null) {
    out$G_null <- G_null
    out$M_null <- calib$M_null
    out$P_null <- calib$P_null
  }
  
  return(out)
}

