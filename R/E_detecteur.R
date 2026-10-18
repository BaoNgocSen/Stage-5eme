#Simulation
sim_changepointPP <- function(tau0, lambda, t_final) {
  
  N_minus <- rpois(n = 1, lambda = lambda[1] * tau0)
  N_plus  <- rpois(n = 1, lambda = lambda[2] * (t_final - tau0))
  
  time_minus <- sort(runif(n = N_minus, min = 0, max = tau0))
  time_plus  <- sort(runif(n = N_plus, min = tau0, max = t_final))

  
  out <- data.frame(
    time  = c(time_minus, time_plus),
    state = as.factor(c(rep(0, N_minus), rep(1, N_plus)))
  )
  return(out)
}
# Estimate lambda
intensity <- function(data,d0,en,W){
  sta<-max(d0,en-W)
  num_event<-nrow(data[(data$time<=en)&(data$time>=sta),])
  return(num_event/(en-sta))
}
# Calcule baseline increment
Ln<- function(data,n,lamb0,d0,W){
  if((n==d0)|| (n==(d0+1)))
    return(1)
  if (nrow(data[(data$time<=n)&(data$time>=(n-1)),])==0)
    return(0)
  else {
    lamb1<-intensity(data,d0,n-1,W)
    num<-nrow(data[(data$time<=n)&(data$time>=(n-1)),])
    L<- exp(lamb0-lamb1)*(lamb1/lamb0)^num
    return(L)
  }
}
# Calcule SR
SR<-function(data,n,lamb0,d0,W,SR_pre=0){
  if (n==d0){
    #cat("L",L, "\n")
    return(0)
  }
  else{
    L<-Ln(data,n,lamb0,d0,W)
    #cat("L",L, "\n")
    return(L*(SR_pre+1))
  }
}
# Alarme
Alarme_SR<-function(data,lamb0,d0,W,alpha,n,Verbose=TRUE){
  E_detec_SR<-0
  new_data<-NULL
  for (i in ((d0+1):n)){
    new_data <- rbind(new_data,data[(data$time>(i-1))&(data$time<=i),])
    SR_pre<- E_detec_SR
    #cat("Time excute:",i, "Value E_detec_SR",E_detec_SR , " ")
    E_detec_SR<-SR(new_data,i,lamb0,d0,W,SR_pre)
    if (E_detec_SR>=1/alpha){
      if (Vesbose==TRUE){
      cat("Alarm at", i, "\n")}
      return(i)
    }
  }
  if (Verbose==TRUE){
  print("No changepoint detection")}
  return(n)
}

