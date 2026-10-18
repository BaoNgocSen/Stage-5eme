#--------------------------------------------------------------------
# Detect break points in a marked Poisson point process with  lambda and mu known
#--------------------------------------------------------------------
# sim_changepointPP_m: simulate events for one ewe and multiple rams
# loglikelihood_m: compute the log-likelihood for one ewe 
# rebuildlambdamu_m: rebuild ram-specific lambda and mu from S 
# rebuildS_L_m: rebuild S by likelihood maximisation
# rebuildS_BH_m: compute Monte Carlo p-values and BH correction
# detectS_BH_m: finally function
#--------------------------------------------------------------------

#' Function that simulates the Marked Poisson Point Process with a break point
#'
#' @param tau0 Numeric . Time of the break point.
#' @param lambda Numeric matrix with one row per ram and two columns.
#'   Column 1: event rates before the break point and column 2: event rates after the break point.
#' @param mu Column 1: exponential mark rates before the break point, column 2:  rates after the break point.
#' @param t_final Numeric. Final observation time.
#' @return A tibble with columns `time`, `marks`, `state`, and `numAlpha`.


sim_changepointPP_m <- function(tau0, lambda, mu, t_final) {
 
  n_ram <- nrow(lambda)
  ram_labels <- rownames(lambda)
  
  list_rams <- lapply(seq_len(n_ram), function(r) {

    N_minus <- rpois(n = 1, lambda = lambda[r, 1] * tau0)
    N_plus  <- rpois(n = 1, lambda = lambda[r, 2] * (t_final - tau0))
    
    time_minus <- if (N_minus > 0) sort(runif(N_minus, 0, tau0)) else numeric(0)
    time_plus  <- if (N_plus  > 0) sort(runif(N_plus, tau0, t_final)) else numeric(0)
    
    marks_minus <- if (N_minus > 0) rexp(N_minus, rate = mu[r, 1]) else numeric(0)
    marks_plus  <- if (N_plus  > 0) rexp(N_plus,  rate = mu[r, 2]) else numeric(0)
    
    tibble::tibble(
      time     = c(time_minus, time_plus),
      marks    = c(marks_minus, marks_plus),
      state    = as.factor(c(rep(0, N_minus), rep(1, N_plus))),
      numAlpha = ram_labels[r]
    )
  })
  
  out <- dplyr::bind_rows(list_rams)
  out <- dplyr::arrange(out, .data$time)
  out
}

#-------------------------------------------------------------------
#' Compute the log-likelihood for one ewe
#'
#' @param lambda Numeric matrix with one row per ram and two columns of event
#'   rates: before and after the candidate break point.
#' @param mu Numeric matrix with one row per ram and two columns of exponential
#'   mark rates: before and after the candidate break point.
#' @param S Numeric scalar. Candidate break-point time.
#' @param time Numeric vector of event times for one ewe, ordered increasingly.
#' @param U Numeric vector of event marks in the same order as `time`.
#' @param ram Vector of ram identifiers in the same order as `time`.
#' @param t_final Numeric scalar. Final observation time.
#' @return A numeric scalar containing the log-likelihood.

loglikelihood_m <- function(lambda, mu, S, time, U, ram, t_final){
  # Standardise ram identifiers and determine the total number of events
  ram <- as.character(ram)
  ntot <- length(time)
  
  if (ntot == 0){
    return(sum((lambda[,2] - lambda[,1]) * S - t_final * lambda[,2]))
  }
  
  L0 <- sum((lambda[,2] - lambda[,1]) * S - t_final * lambda[,2])
  
  idx_before <- which(time < S)
  idx_equal  <- which(time == S)
  idx_after  <- which(time > S)
  
  # Case 1: events occurring exactly at S belong to the pre-break 
  idx_before_1 <- c(idx_before, idx_equal)
  idx_after_1  <- idx_after
  
  L_right <- L0
  
  if (length(idx_before_1) > 0){
    L_right <- L_right +
      sum(log(lambda[ram[idx_before_1], 1] * mu[ram[idx_before_1], 1])) -
      sum(mu[ram[idx_before_1], 1] * U[idx_before_1], na.rm = TRUE)
  }
  
  if (length(idx_after_1) > 0){
    L_right <- L_right +
      sum(log(lambda[ram[idx_after_1], 2] * mu[ram[idx_after_1], 2])) -
      sum(mu[ram[idx_after_1], 2] * U[idx_after_1], na.rm = TRUE)
  }
  
  # Without an event exactly at S, the allocation is unambiguous
  if (length(idx_equal) == 0){
    return(L_right)
  }
  
  # Case 2: events occurring exactly at S belong to the post-break regime
  idx_before_2 <- idx_before
  idx_after_2  <- c(idx_equal, idx_after)
  
  L_left <- L0

  if (length(idx_before_2) > 0){
    L_left <- L_left +
      sum(log(lambda[ram[idx_before_2], 1] * mu[ram[idx_before_2], 1])) -
      sum(mu[ram[idx_before_2], 1] * U[idx_before_2], na.rm = TRUE)
  }
  
  if (length(idx_after_2) > 0){
    L_left <- L_left +
      sum(log(lambda[ram[idx_after_2], 2] * mu[ram[idx_after_2], 2])) -
      sum(mu[ram[idx_after_2], 2] * U[idx_after_2], na.rm = TRUE)
  }
  
  return(max(L_right, L_left))
}

#-----------------------------------------------------------
#' Rebuild ram-specific lambda and mu parameters
#'
#' @param datain Input data frame containing  `id`, `time`, `marks` and `numAlpha`.
#' @param S Numeric vector of estimated break points.
#' @param t_final Numeric scalar. Final observation time.
#' @return A list with two elements: `lambda`, an `n_ram` by 2 matrix of event
#'   rates, and `mu`, an `n_ram` by 2 matrix of mark rates. Columns correspond
#'   to the periods before and after the break point.

rebuildlambdamu_m <- function(datain, S, t_final){
  ids <- sort(unique(datain$id))
  nb_indiv <- length(ids)
  ram_levels <- sort(unique(as.character(datain$numAlpha)))
  n_ram <- length(ram_levels)

  SommeS <- 0
  SommeSautAvant <- numeric(n_ram)
  SommeSautApres <- numeric(n_ram)
  SommeUAvant <- numeric(n_ram)
  SommeUApres <- numeric(n_ram)
  
  names(SommeSautAvant) <- ram_levels
  names(SommeSautApres) <- ram_levels
  names(SommeUAvant) <- ram_levels
  names(SommeUApres) <- ram_levels
  

  for (ii in seq_along(ids)){
    ind <- ids[ii]
    ewe_ind <- datain[datain$id == ind, ]
    SommeS <- SommeS + S[ii]
    
    if (nrow(ewe_ind) > 0){
      ewe_ind <- ewe_ind[order(ewe_ind$time), ]
      
      sauts <- ewe_ind$time
      U <- ewe_ind$marks
      ram <- as.character(ewe_ind$numAlpha)
      
  
      idx_avant <- which(sauts < S[ii])
      idx_apres <- which(sauts >= S[ii])

      if (length(idx_avant) > 0){
        ram_avant <- ram[idx_avant]
        U_avant <- U[idx_avant]
        
        tab_avant <- table(factor(ram_avant, levels = ram_levels))
        SommeSautAvant <- SommeSautAvant + as.numeric(tab_avant)
        
        for (j in seq_along(idx_avant)){
          SommeUAvant[ram_avant[j]] <- SommeUAvant[ram_avant[j]] + U_avant[j]
        }
      }
  
      if (length(idx_apres) > 0){
        ram_apres <- ram[idx_apres]
        U_apres <- U[idx_apres]
        
        tab_apres <- table(factor(ram_apres, levels = ram_levels))
        SommeSautApres <- SommeSautApres + as.numeric(tab_apres)
        
        for (j in seq_along(idx_apres)){
          SommeUApres[ram_apres[j]] <- SommeUApres[ram_apres[j]] + U_apres[j]
        }
      }
    }
  }
  
  denom_before_time <- SommeS
  denom_after_time <- nb_indiv * t_final - SommeS
  
  if (denom_before_time <= 0) denom_before_time <- 1
  if (denom_after_time <= 0) denom_after_time <- 1
  
  # Avoid 0/0
  for (r in seq_len(n_ram)){
    if (SommeSautAvant[r] == 0) SommeUAvant[r] <- 1
    if (SommeSautApres[r] == 0) SommeUApres[r] <- 1
  }
  
  lambdan <- matrix(0.01, nrow = n_ram, ncol = 2)
  mun <- matrix(0.01, nrow = n_ram, ncol = 2)
  
  rownames(lambdan) <- ram_levels
  rownames(mun) <- ram_levels
  colnames(lambdan) <- c("before", "after")
  colnames(mun) <- c("before", "after")

  for (r in seq_len(n_ram)){
    lambdan[r, 1] <- max(SommeSautAvant[r] / denom_before_time, 0.01)
    lambdan[r, 2] <- max(SommeSautApres[r] / denom_after_time, 0.01)
    
    mun[r, 1] <- max(SommeSautAvant[r] / SommeUAvant[r], 0.01)
    mun[r, 2] <- max(SommeSautApres[r] / SommeUApres[r], 0.01)
  }
  
  return(list(lambda = lambdan, mu = mun))
}

#-----------------------------------------------------------
#' Rebuild break-point times by likelihood maximisation
#'
#' @param lambdan Numeric matrix with one row per ram and two columns of
#'   current event-rate estimates.
#' @param mun Numeric matrix with one row per ram and two columns of current
#'   mark-rate estimates.
#' @param ewe Data frame containing `id`, `time`, `marks`, and  `numAlpha`.
#' @param t_final Numeric scalar. Final observation time.
#' @return  the  break-point vector and the  likelihood differences `l1 - l0`.

rebuildS_L_m <- function(lambdan, mun, ewe, t_final){
 
  S_reconstruit <- c()
  l0 <- c()
  l1 <- c()
  
  ids <- sort(unique(ewe$id))
  nb_indiv <- length(ids)
  
  for (ii in seq_along(ids)){
    ind <- ids[ii]
    ewe_ind <- ewe[ewe$id == ind, ]
    ntot <- nrow(ewe_ind)
    
    if (ntot == 0){
      S_reconstruit <- c(S_reconstruit, t_final)
      l0 <- c(l0, 1)
      l1 <- c(l1, 1)} 
    else {
      ewe_ind <- ewe_ind[order(ewe_ind$time), ]
      time_i <- ewe_ind$time
      U_i    <- ewe_ind$marks
      ram_i  <- as.character(ewe_ind$numAlpha)
      
      S <- time_i[1]
      likelihoodopt <- loglikelihood_m(lambda = lambdan, mu = mun,  S = S,
        time = time_i,U = U_i, ram = ram_i,t_final = t_final )
      Sopt <- S

      if (ntot > 1){
        for (k in 2:ntot){
          S <- time_i[k]
          
          likelihoodk <- loglikelihood_m( lambda = lambdan, mu = mun,S = S,
                                          time = time_i, U = U_i,ram = ram_i, t_final = t_final )
          
          if (likelihoodk > likelihoodopt){
            likelihoodopt <- likelihoodk
            Sopt <- S
          } 
        }  
      }
      
      S_reconstruit <- c(S_reconstruit, Sopt)
      
      l0_i <- loglikelihood_m( lambda = lambdan,mu = mun,S = t_final,
                          time = time_i, U = U_i, ram = ram_i,t_final = t_final)
  
      l0 <- c(l0, l0_i)
      l1 <- c(l1, likelihoodopt)
    }
  }
  return(list(S_reconstruit, l1 - l0))
}

#-----------------------------------------------------------------------
#' Rebuild break points using Monte Carlo p-values and BH correction
#'
#' @param datain Input data frame containing `id`, `time`, `marks`, and `numAlpha`.
#' @param alpha Numeric scalar. Target threshold
#' @param nb_MC Positive integer. Number of Monte Carlo samples
#' @param lambdan Numeric matrix with one row per ram and two event-rate columns.
#' @param mun Numeric matrix with one row per ram and two mark-rate columns.
#' @param t_final Numeric scalar. Final observation time.
#' @return A named list containing `S`, `pvalue`, `pvalueclasse`, and `critBH`.

rebuildS_BH_m <- function(datain, alpha, nb_MC, lambdan, mun, t_final){
  # Parametric Monte Carlo simulation under H0: tau0 equals t_final
  list_individuals <- lapply(seq_len(nb_MC), function(i){
    tau0 <- t_final
    ewe_MC <- sim_changepointPP_m(
      tau0    = tau0,
      lambda  = lambdan,
      mu      = mun,
      t_final = t_final
    )
    
    if (nrow(ewe_MC) > 0) {
      ewe_MC$id <- i
      ewe_MC$tau <- tau0
    } else {
      ewe_MC <- tibble::tibble(
        time     = numeric(0), 
        marks    = numeric(0),
        state    = factor(levels = c(0, 1)),
        numAlpha = if (is.null(rownames(lambdan))) character(0) else character(0),
        id       = integer(0),
        tau      = numeric(0)
      )
    }
    ewe_MC
  })
  
  ewe_MC <- do.call(rbind, list_individuals)
  
  Sl_data <- rebuildS_L_m(lambdan, mun, datain, t_final)
  Sl_MC   <- rebuildS_L_m(lambdan, mun, ewe_MC, t_final)
  
  l_data <- Sl_data[[2]]  
  l_MC   <- Sl_MC[[2]]    
  S      <- Sl_data[[1]]  
  
  ids_data <- sort(unique(datain$id))
  nb_indiv <- length(ids_data)
  
  pvalue <- numeric(nb_indiv)
  
  for (ii in seq_len(nb_indiv)){
    pvalue[ii] <- sum(l_MC > l_data[ii]) / nb_MC
  }
  
  pvalueclasse <- sort(pvalue, index.return = TRUE)
  
  critBH <- numeric(nb_indiv)
  BH <- 0
  
  for (ind in seq_len(nb_indiv)){
    if (BH == 0 && (pvalueclasse$x[ind] > ind * alpha / nb_indiv)){
      BH <- 1
    }
    if (BH == 1){
      S[pvalueclasse$ix[ind]] <- t_final
    }
    critBH[ind] <- BH
  }
  
  resu <- list(S, pvalue, pvalueclasse, critBH)
  names(resu) <- c("S", "pvalue", "pvalueclasse", "critBH")
  return(resu)
}

#-----------------------------------------------------------
#' Detect break points with ram-specific parameters and BH correction
#'
#' @param datain Input data frame containing `id`, `time`, `marks`, and  `numAlpha`.
#' @param lambda0 Numeric matrix with one row per ram and two columns containing
#'   the initial event-rate estimates before and after the break point.
#' @param mu0 Numeric matrix with one row per ram and two columns containing the
#'   initial mark-rate estimates before and after the break point.
#' @param nprec Positive integer. Number of parameter-reconstruction iterations.
#' @param alpha Numeric scalar. Threshold used by the Benjamini-Hochberg rule.
#' @param nb_MC Positive integer. Number of Monte Carlo samples per iteration.
#' @param t_final Numeric scalar. Final observation time.
#' @return A named list with three elements: `outBH`, containing break points,
#'   p-values and BH results; `lambda_n`, containing the final event-rate matrix;
#'   and `mu_n`, containing the final mark-rate matrix.

detectS_BH_m <- function(datain, lambda0, mu0, nprec, alpha, nb_MC, t_final){
  lambdan <- lambda0
  mun <- mu0
  
  for (i in seq_len(nprec)){
    outBH <- rebuildS_BH_m(datain, alpha, nb_MC, lambdan, mun, t_final)
    names(outBH) <- c("S", "p.value", "p.value.classe", "BH")
    
    S <- outBH[[1]]
    
    Param <- rebuildlambdamu_m(datain, S, t_final)
    lambdan <- Param$lambda
    mun <- Param$mu
  }
  
  cat("Estimated lambda:\n")
  print(round(lambdan, 6))
  cat("\nEstimated mu:\n")
  print(round(mun, 6))
  cat("\nEstimated S:\n")
  print(round(outBH$S, 6))
  cat("\nP-values:\n")
  print(round(outBH$p.value, 6))
  cat("\n")

  output <- list(outBH, lambdan, mun)
  names(output) <- c("outBH", "lambda_n", "mu_n")
  return(output)
}
