model{
  for (i in 1:model_obs){
    
    for (num_cutpt in 1:3){
      Z[i,num_cutpt] <- ( alpha_SD*SD[i] 
                          + alpha_N*Nitrogen[i]
                          + alpha_P*Phosphorus[i]
                          + alpha_DIN*DIN[i]
                          + alpha_DIP*DIP[i]
                          + alpha_SubR[Subregion[i]]
                          - C[num_cutpt])
      / s
      
      Q[i,num_cutpt] <- 1/(1+exp(-Z[i,num_cutpt]))
    }
    
    P[i,1] <- max(min(1 - Q[i,1],1),0) + 0.0001
    P[i,2] <- Q[i,1] - Q[i,2] + 0.0001
    P[i,3] <- Q[i,2] - Q[i,3] + 0.0001
    P[i,4] <- max(min(Q[i,3],1),0) + 0.0001
    
    # Farnaz - what is this line 23 doing?
    TS[i] ~ dcat(P[i,])
  }
  
  # PRIORS
  # Proportional Odds Logistic Regression (POLR)
  # POLR Coefficients
  # If using posterior distribution from 2010 as priors, these uninformed priors needs to be changed
  # to posteriors from 2010

  #alpha_SD ~ dnorm(0,.0001)
  #alpha_N ~ dnorm(0,.0001)
  #alpha_P ~ dnorm(0,.0001)
  #alpha_DIN ~ dnorm(0,.0001)
  #alpha_DIP ~ dnorm(0,.0001)

  alpha_SD ~ dnorm(mn_alpha_SD,tau_alpha_SD)
  alpha_N ~ dnorm(mn_alpha_N,tau_alpha_N)
  alpha_P ~ dnorm(mn_alpha_P,tau_alpha_P)
  alpha_DIN ~ dnorm(mn_alpha_DIN,tau_alpha_DIN)
  alpha_DIP ~ dnorm(mn_alpha_DIP,tau_alpha_DIP)

  for (j in 1:num_subregions){
    #alpha_SubR[j] ~ dnorm(0,.0001)
    alpha_SubR[j] ~ dnorm(mn_alpha_SubR[j],tau_alpha_SubR[j])
  }
  
  s ~ dlnorm(mu.log.s, tau.log.s)
  mu.log.s ~ dnorm(mn_s,0.0001)
  #mu.log.s ~ dnorm(0,0.0001)
  tau.log.s <- pow(sigma.log.s,-2)
  #tau.log.s <- tau_log_s
  sigma.log.s ~ dunif(0,1000)
  
  for (i in 1:3) {
    #cutpt_raw[i] ~ dnorm(0, .0001)
    cutpt_raw[i] ~ dnorm(mn_cut_pts[i],tau_cut_pts[i])
  }
  C <- sort(cutpt_raw)
  
}