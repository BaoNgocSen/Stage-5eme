# E-process ONLINE DETECTION FOR ONE EWE

### Compute daily likelihood-ratio increment
### This function computes log L^d using data in the interval [d-1, d).
compute_logLR_one_day <- function(ewe_ind, d,lambda0, mu0,lambda1, mu1) {
  
  dat_d <- ewe_ind[ewe_ind$time >= (d - 1) & ewe_ind$time < d, ]
  N_d <- nrow(dat_d)
  V_d <- if (N_d == 0) 0 else sum(dat_d$marks)
  logLR_d <- N_d * log((lambda1 * mu1) / (lambda0 * mu0)) -(lambda1 - lambda0) - (mu1 - mu0) * V_d
  return(list(day = d,N_d = N_d,V_d = V_d,logLR_d = logLR_d))
}

### This function updates log Y_d = sum_{t=1}^d log L^t and alarms when Y_d >= log( 1 / alpha.)
detect_martingale_one_ewe <- function(ewe_ind,K,lambda0, mu0,lambda1, mu1,alpha = 0.05,verbose = TRUE) {
  
  log_threshold <- log(1 / alpha)
  path <- data.frame(
    day = seq_len(K),
    N_d = NA_integer_,
    V_d = NA_real_,
    logLR_d = NA_real_,
    logY = NA_real_,
    Y = NA_real_,
    alarm_here = FALSE
  )
  
  logY <- 0
  alarm <- FALSE
  tau <- Inf
  
  for (d in seq_len(K)) {
    
    out_d <- compute_logLR_one_day(ewe_ind = ewe_ind,d = d,lambda0 = lambda0,
      mu0 = mu0,lambda1 = lambda1,mu1 = mu1)
    
    logY <- logY + out_d$logLR_d
    
    path$N_d[d] <- out_d$N_d
    path$V_d[d] <- out_d$V_d
    path$logLR_d[d] <- out_d$logLR_d
    path$logY[d] <- logY
    path$Y[d] <- exp(logY)
    
    if (logY >= log_threshold) {
      alarm <- TRUE
      tau <- d
      path$alarm_here[d] <- TRUE
      break
    }
  }
  
  return(list(alarm = alarm,tau = tau,alpha = alpha,threshold = 1 / alpha,path = path))
}

### Bonferroni correction: each ewe is tested at level alpha / n_ewe
detect_martingale_multiple_ewes <- function(data_all,  K,lambda0, mu0, lambda1, mu1, alpha = 0.05,ewe_col = "ewe",
                                            n_ewe = NULL,verbose = TRUE) {
  
  # If n_ewe is not provided, compute it from the data
  if (is.null(n_ewe)) {
    n_ewe <- length(unique(data_all[[ewe_col]]))
  }
  
  alpha_each <- alpha / n_ewe
  ewe_ids <- sort(unique(data_all[[ewe_col]]))
  
  results <- list()
  summary_table <- data.frame(
    ewe = ewe_ids,
    alarm = FALSE,
    tau = Inf,
    max_Y = NA_real_,
    max_logY = NA_real_
  )
  
  for (idx in seq_along(ewe_ids)) {
    
    ewe_id <- ewe_ids[idx]
    ewe_ind <- data_all[data_all[[ewe_col]] == ewe_id, ]
    
    res_i <- detect_martingale_one_ewe(
      ewe_ind = ewe_ind,
      K = K,
      lambda0 = lambda0,
      mu0 = mu0,
      lambda1 = lambda1,
      mu1 = mu1,
      alpha = alpha_each,
      verbose = FALSE
    )
    
    results[[as.character(ewe_id)]] <- res_i
    
    summary_table$alarm[idx] <- res_i$alarm
    summary_table$tau[idx] <- res_i$tau
    summary_table$max_Y[idx] <- max(res_i$path$Y, na.rm = TRUE)
    summary_table$max_logY[idx] <- max(res_i$path$logY, na.rm = TRUE)
  }
  
  global_alarm <- any(summary_table$alarm)
  
  if (verbose) {
    print(summary_table)
    
    if (global_alarm) {
      cat("\nGLOBAL ALARM: YES\n")
      cat("At least one ewe/process raised an alarm.\n")
    } else {
      cat("\nGLOBAL ALARM: NO\n")
      cat("No ewe/process raised an alarm before K =", K, "\n")
    }
  }
  
  return(list(
    global_alarm = global_alarm,
    alpha = alpha,
    n_ewe = n_ewe,
    alpha_each = alpha_each,
    threshold_each = 1 / alpha_each,
    summary = summary_table,
    results = results
  ))
}