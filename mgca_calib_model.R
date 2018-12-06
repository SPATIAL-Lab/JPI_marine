model {
  
  for(i in 1:length(t_m)){
    mgca[i] ~ dnorm(mgca_m[i], mgca_pre)
    mgca_m[i] = (a[1] + a[2] * t[i]) * mgca_sw[i] ^ a[3]
    
    t[i] ~ dnorm(t_m[i], 1 / t_sd[i]^2)
    mgca_sw[i] ~ dnorm(mgca_sw_m[i], 1 / mgca_sw_sd[i]^2) I (0.1,)
  }
  
  mgca_pre ~ dgamma(3, 1 / 30)
  
  #Using very loose constraints on MgCa calib parameters
  a[1] ~ dunif(0.9, 1.1)
  a[2] ~ dunif(-0.05, 0.35)
  a[3] ~ dunif(-0.02, 0.02)   
}