source("Online.R")
#--------------------------------------------------------------------
# Online changepoint detection using a marked Poisson point 
# process with ram-specific parameters with known parameters          

# sim_changepointPP_ram  Simulate one ewe with ram effects
# sim_herd_changepoint_ram            Simulate a complete herd
# compute_G_one_day_ram               Compute G(d) for one ewe
# bootstrap_G_null_herd_ram           Simulate herd-level null statistics
# calibrate_min_p_null_herd           Calibrate the whole-herd FWER threshold
# compute_pvalues_herd_one_day_ram    Compute daily p-values for active ewes
# detect_rupture_herd_online_FWER_ram Run the complete online procedure
#--------------------------------------------------------------------'

# Simuler un processus de Poisson avec un point de rupture et un effet  bélier.
# Simulation les processus de Poisson avec :
# tau0 : 1 rupture
# lambda : liste des intensités
# mu : liste des paramètres des marques
# t_final : temps final d'observation

sim_changepointPP_ram <- function(tau0, lambda, mu, t_final) {
  
  n_ram <- nrow(lambda)
  ram_labels <- rownames(lambda)
  out <- lapply(seq_len(n_ram), function(r) {
    
    N_minus <- rpois(n = 1, lambda = lambda[r, 1] * tau0)
    N_plus  <- rpois(n = 1, lambda = lambda[r, 2] * (t_final - tau0))
    
    time_minus <- if (N_minus > 0) sort(runif(N_minus, 0, tau0)) else numeric(0)
    time_plus  <- if (N_plus  > 0) sort(runif(N_plus, tau0, t_final)) else numeric(0)
    
    marks_minus <- if (N_minus > 0) rexp(N_minus, rate = mu[r, 1]) else numeric(0)
    marks_plus  <- if (N_plus  > 0) rexp(N_plus,  rate = mu[r, 2]) else numeric(0)
    
    data.frame(
      time     = c(time_minus, time_plus),
      marks    = c(marks_minus, marks_plus),
      state    = as.factor(c(rep(0, N_minus), rep(1, N_plus))),
      numAlpha = rep(ram_labels[r], N_minus + N_plus)
    )
  })
  
  out<-bind_rows(out)
}

#--------------------------------------------------------------------
#'Simuler un troupeau avec des ruptures à chaque brebis et un effet du bélier.
# tau_vec : vecteur des temps de rupture, un par brebis.
# lambda0, lambda1 : intensités avant/après la rupture.
# mu0, mu1 : paramètres des marques avant/après la rupture.
# t_final : temps final d'observation.
# return Un data frame contenant les événements simulés pour toutes les brebis.

sim_herd_changepoint_ram <- function(tau_vec, lambda,mu,t_final) {
  
  n_ewe <- length(tau_vec)
  out_list <- vector("list", n_ewe)
  
  for (i in seq_len(n_ewe)) {
    
    dat_i <- sim_changepointPP_ram(tau_vec[i],lambda,mu,t_final)
    dat_i$ID <- i
    out_list[[i]] <- dat_i
  }
  
  out <- do.call(rbind, out_list)
  out <- out[order(out$ID, out$time), ]
  return(out)
}

#--------------------------------------------------------------------
#Calculer la statistique journalière pour une brebis avec effet du bélier
# ewe_ind : données contenant au moins `time`, `marks` et numAlpha.
# d : jour de surveillance actuel.
# lambda0, lambda1 : intensités avant/après la rupture.
# mu0, mu1 : paramètres des marques avant/après la rupture.
# return Une liste contenant :
#'   `G` : statistique maximale du rapport de vraisemblance ;
#'   `S_hat` : point de rupture candidat maximisant la statistique.
#--------------------------------------------------------------------
compute_G_one_day_ram <- function(ewe_ind, d,lambda,mu) {
  
  ram_levels <- rownames(lambda)
  R <- length(ram_levels)
  ewe_ind <- ewe_ind[ewe_ind$time < d, ]
  ewe_ind <- ewe_ind[order(ewe_ind$time), ]
  times <- ewe_ind$time
  marks <- ewe_ind$marks
  rams  <- ewe_ind$numAlpha
  log_ratio <- log((lambda[, 2]  * mu[, 2]) / (lambda[, 1] * mu[, 1]))
  dlambda   <- lambda[, 2] - lambda[, 1]
  dmu       <- mu[, 2] - mu[, 1]
  S_grid <- times[times < d]
  S_grid <- unique(S_grid)
  S_grid <- times[times < d]
  
  if (length(S_grid) == 0) {
    return(list( G = 0,S_hat = NA_real_))
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
      val <- val + N_r * log_ratio[rr] -dlambda[rr] * (d - S) -  dmu[rr] * V_r
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
#Simuler par bootstrap la statistique GLR nulle pour tout le troupeau
# B : nombre de simulations de Monte Carlo.
# n_ewe : nombre de brebis dans chaque troupeau simulé.
# K : nombre de jours de surveillance et temps final de simulation.
# lambda0, mu0 : paramètres nuls des événements et des marques.
# lambda1, mu1 : paramètres alternatifs utilisés dans la GLR.
# seed : graine aléatoire 
# return: Un tableau numérique de dimension `B x n_ewe x K`, où
# `[b, i, d]` correspond à la statistique nulle pour la simulation `b`,
# la brebis `i` et le jour `d`.
bootstrap_G_null_herd_ram <- function(B, n_ewe, K, lambda, mu,
                                     seed = 123) {
  
  set.seed(seed)
  ram_levels <-rownames(lambda)
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
      ### Under H0: no changepoin S=K
      dat_bi <- sim_changepointPP_ram(K,lambda,mu,K)
      for (d in seq_len(K)) {
        out_d <- compute_G_one_day_ram(dat_bi, d,lambda, mu)
        G_null[b, i, d] <- out_d$G
      }
    }
  }
  
  return(G_null)
}

#--------------------------------------------------------------------
# Calibrer le seuil FWER en ligne pour l'ensemble du troupeau

# G_null : tableau numérique de dimension `B x n_ewe x K`.
# d_min : premier jour pris en compte pour la calibration séquentielle.
# alpha : niveau de FWER souhaité pour l'ensemble du troupeau.
# return Une liste contenant :
#   `P_null` : tableau des p-values sous le modèle nul, de même dimension
#              que `G_null` ;
#   `M_null` : p-value minimale pour chaque simulation bootstrap ;
#  `u_alpha` : seuil d'alarme calibré pour l'ensemble du troupeau.
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
#Calculer les p-values journalières pour les brebis encore surveillées
# dat_herd : données du troupeau contenant `ID`, `time`, `marks` et `numAlpha`.
# ids_active : identifiants des brebis encore surveillées.
# d : jour de surveillance actuel.
# G_null : tableau des statistiques sous le modèle nul,
#          de dimension `B x n_ewe x K`.
# lambda, mu0 : paramètres avant la rupture.
# lambda1, mu1 : paramètres après la rupture.
# return Un data frame avec une ligne par brebis surveillée et les colonnes
# `ID`, `day`, `G_obs`, `p_value` et `S_hat`.
compute_pvalues_herd_one_day_ram <- function(dat_herd, ids_active, d,  G_null,
                                             lambda, mu) {
      
  G_day_null <- as.vector(G_null[, , d])
  Npool <- length(G_day_null)
  
  out <- data.frame(ID = ids_active,day = d, G_obs = NA_real_, 
                    p_value = NA_real_,S_hat = NA_real_ )
  
  for (a in seq_along(ids_active)) {
    
    id_i <- ids_active[a]
    dat_i <- dat_herd[dat_herd$ID == id_i, ]
    G_out <- compute_G_one_day_ram(dat_i,d,lambda, mu)
    p_i_d <- sum(G_day_null >= G_out$G) / Npool
    out$G_obs[a] <- G_out$G
    out$p_value[a] <- p_i_d
    out$S_hat[a] <- G_out$S_hat
  }
  
  return(out)
}

#--------------------------------------------------------------------
# Exécuter la détection en ligne pour du troupeau avec effet bélier
# dat_herd : données contenant `ID`, `time`, `marks` et numAlpha.
# K : dernier jour de surveillance.
# lambda0, mu0 : paramètres avant la rupture.
# lambda1, mu1 : paramètres après la rupture.
# B : nombre de simulations de Monte Carlo si `G_null` n'est pas fourni.
# alpha : niveau de FWER souhaité pour l'ensemble du troupeau.
# d_min : premier jour où une alarme peut être déclenchée.
# seed : graine aléatoire 
# verbose : affiche un résumé lorsque `TRUE`.
# ram_levels : labels des béliers et leur ordre (optionnel).
# @return Une liste contenant :
#   `result` : résumé des alarmes, une ligne par brebis ;
#   `u_alpha` : seuil calibré pour l'ensemble du troupeau ;
# `G_null`, `M_null` et `P_null`.
detect_rupture_herd_online_FWER_ram <- function(dat_herd, K, lambda, mu,
                                                B = 1000, alpha = 0.05,
                                                d_min = 1, seed = 123,
                                                verbose = TRUE) {
  
  ids <- sort(unique(dat_herd$ID))
  n_ewe <- length(ids)
  set.seed(seed)
  ram_levels <- rownames(lambda)
  G_null <- bootstrap_G_null_herd_ram( B, n_ewe, K, lambda, mu, seed )
  calib <- calibrate_min_p_null_herd(G_null, d_min, alpha )
  u_alpha <- calib$u_alpha
  result <- data.frame(
    ID = ids,
    alarm = FALSE,
    tau = K,
    S_hat = K,
    p_at_alarm = NA_real_,
    G_at_alarm = NA_real_
  )
  
  active_ids <- ids
  for (d in seq_len(K)) {
    if (length(active_ids) == 0) break
    day_out <- compute_pvalues_herd_one_day_ram(dat_herd,active_ids,d,G_null,
                                  lambda,mu)
    if (d >= d_min) {
      alarm_ids <- day_out$ID[ day_out$p_value <= u_alpha ]
      if (length(alarm_ids) > 0) {
        for (id_alarm in alarm_ids) {
          row_day <- day_out[day_out$ID == id_alarm,]
          idx_res <- which(result$ID == id_alarm)
          result$alarm[idx_res] <- TRUE
          result$tau[idx_res] <- d
          result$S_hat[idx_res] <- row_day$S_hat
          result$p_at_alarm[idx_res] <- row_day$p_value
          result$G_at_alarm[idx_res] <- row_day$G_obs
        }
        
        active_ids <- setdiff(active_ids, alarm_ids )
      }
    }
  }
  
  if (verbose) {
    if (sum(result$alarm) > 0) {
      print(result[result$alarm, ])
    }
  }
  
  out <- list(
    result = result,
    u_alpha = u_alpha,
    G_null = G_null,
    M_null = calib$M_null,
    P_null = calib$P_null
  )
  
  return(out)
}

