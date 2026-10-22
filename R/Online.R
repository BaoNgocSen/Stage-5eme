
#--------------------------------------------------------------------
# Objectif : détecter en ligne la rupture pour une brebis et un bélier (mark=1)

# sim_changepointPP : simuler des événements d'un processus de Poisson marqué
# compute_G_one_day : calculer la statistique de détection G(d) pour une journée
# bootstrap_G_null : simuler la distribution nulle de G(d)
# calibrate_min_p_null : calibrer le seuil des p-values séquentielles
# compute_pvalue_one_day : calculer la p-value pour une brebis sur une journée
# detect_rupture_one_ewe_online : exécuter toute la procédure 
#--------------------------------------------------------------------

# Simuler un processus de Poisson marqué avec une rupture
# tau0 : la rupture.
# lambda (lambda[1],lambda[2]) : intensités avant et après tau0.
# mu (mu[1], mu[2]) : paramètres des lois exponentielles  avant et après tau0.
# t_final : temps final d'observation.

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
  
  return(out)
}

#--------------------------------------------------------------------
# Calculer la statistique en ligne pour une brebis et une journée

# ewe_ind : Data frame contenant les colonnes `time` et `marks`.
# d : Jour de surveillance actuel.
# lambda0/lambda1 : intensités avant et après le point de rupture.
# mu0/mu1 : paramètres des lois des marques avant et après le point de rupture.
# Return (G, S_hat) : statistique de test et la rupture estimé.

compute_G_one_day <- function(ewe_ind, d,lambda0, mu0,lambda1, mu1) {
  
  ewe_ind <- ewe_ind[ewe_ind$time < d, ]
  ewe_ind <- ewe_ind[order(ewe_ind$time), ]
  
  times <- ewe_ind$time
  marks <- ewe_ind$marks
  
  log_ratio <- log((lambda1 * mu1) / (lambda0 * mu0))
  dlambda   <- lambda1 - lambda0
  dmu       <- mu1 - mu0
  
  S_grid <- c(0, times[times < d])
  
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

#--------------------------------------------------------------------
# Simuler la distribution nulle de la statistique journalière

# B : Nombre de simulations de Monte Carlo.
# K : Nombre de jours de surveillance.
# lambda0/lambda1 : Intensités avant et après la rupture.
# mu0/mu1 : Paramètres des lois des marques avant et après la rupture.
# seed : Graine aléatoire facultative (entier).
#  @return Distribution nulle de la statistique journalière G(d).
bootstrap_G_null <- function(B, K, lambda0, mu0, lambda1, mu1, seed = 123) {
  
  G_null <- matrix(0, nrow = B, ncol = K)
  colnames(G_null) <- paste0("day_", seq_len(K))
  
  for (b in seq_len(B)) {
    dat_b <- sim_changepointPP( tau0  = K,  lambda  = c(lambda0, lambda0), mu = c(mu0, mu0), t_final = K)
    for (d in seq_len(K)) {
      out_d <- compute_G_one_day(dat_b,d,lambda0, mu0,lambda1,mu1)
      G_null[b, d] <- out_d$G
    }
  }
  return(G_null)
}

#--------------------------------------------------------------------
# Calibrer le seuil de p-value minimale séquentielle

# G_null : matrice numérique retournée par `bootstrap_G_null()`.
# d_min : premier jour de surveillance.
# alpha : niveau de fausse alarme.
#return Une liste contenant :
#         P_null : matrice des p-values ;
#          R_null : p-value minimale pour chaque trajectoire sous le modèle nul ;
#          u_alpha : seuil calibré.
    
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
  
  return(list( P_null = P_null, R_null = R_null, u_alpha = u_alpha))
}


#--------------------------------------------------------------------
# Calculer la p-value empirique pour une brebis un jour donné

# ewe_ind : données correspondant à une seule brebis.
# d : jour de surveillance.
# G_null : matrice numérique de la statistique nulle,
#          avec une colonne par jour.
# lambda0/lambda1 : intensités dans le modèle avant/après le point de rupture.
# @return Une liste contenant : le jour, la statistique observée,
#'        la p-value empirique et le point de rupture estimé.

compute_pvalue_one_day <- function(ewe_ind, d, G_null, lambda0,mu0,lambda1,mu1) {
  B <- nrow(G_null)
  G_out <- compute_G_one_day( ewe_ind, d,lambda0,mu0, lambda1,mu1)
  
  p_d <- (sum(G_null[, d] >= G_out$G)) / B

  return(list(day = d,G_obs = G_out$G, p_value = p_d, S_hat = G_out$S_hat ))
}


#--------------------------------------------------------------------
# Fonction finale : détecter en ligne un point de rupture pour une brebis

# ewe_ind : données correspondant à une seule brebis.
# B : nombre de simulations de Monte Carlo.
# K : nombre de jours de surveillance.
# lambda0/lambda1 : intensités avant et après le point de rupture.
# d_min : premier jour de surveillance.
# alpha : niveau de fausse alarme.
# seed : graine aléatoire facultative (entier).
# verbose : affiche un résumé de la détection lorsque `TRUE`.

# return Une liste contenant :
#         `alarm` : indique si une alarme de rupture a été déclenchée ;
#         `tau` : premier jour de déclenchement de l'alarme,
#                ou `Inf` si aucune alarme n'a été déclenchée ;
#         `S_hat` : point de rupture estimé au jour de l'alarme ;
#         `u_alpha` : seuil calibré des p-values séquentielles ;
#         `alpha` : niveau de fausse alarme choisi ;
#         `path` : statistiques journalières, p-values,
#'                estimations des points de rupture et indicateurs d'alarme.  

detect_rupture_one_ewe_online <- function(ewe_ind, K,lambda0,mu0,lambda1, mu1,B = 1000,alpha = 0.05,
                                          d_min = 1,seed = 123, verbose = TRUE) {
  
  G_null <- bootstrap_G_null(B, K,lambda0,mu0, lambda1, mu1,seed = seed )
  B0 <- nrow(G_null)
  
  #  Calibrate min-p threshold
  calib <- calibrate_min_p_null(G_null,d_min,alpha )
  u_alpha <- calib$u_alpha
  
  #  Online monitoring
  path <- data.frame(day = seq_len(K), G_obs = NA_real_, p_value = NA_real_,
                     S_hat = NA_real_, alarm_here = FALSE )
  
  alarm <- FALSE
  tau <- Inf
  S_alarm <- NA_real_
  
  for (d in seq_len(K)) {
    G_out <- compute_G_one_day( ewe_ind, d,lambda0,mu0, lambda1, mu1)
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
    }
  }
  
  out <- list(alarm,tau, S_alarm, u_alpha, alpha,path)
  out$G_null <- G_null
  out$R_null <- calib$R_null
  out$P_null <- calib$P_null

  
  return(out)
}

