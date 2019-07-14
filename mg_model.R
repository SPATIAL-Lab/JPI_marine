model {
  
  #Data model for seawater MgCa observations
  
  for(i in 1:length(MgCa_sw)){
    MgCa_sw[i] ~ dnorm(MgCa_sw_m[MgCa_sw.age.ind[i]], 1 / MgCa_sw.sd[i] ^ 2)
    
  }
  
  #Process model for MgCa_sw timeseries
  
  for(i in 2:nmgca.ages){
    MgCa_sw_m[i] = MgCa_sw_m[i-1] * ((MgCa_sw_m.eps[i]) + 1)
    
    MgCa_sw_m.eps[i] ~ dnorm(MgCa_sw_m.eps[i-1] * MgCa_sw_m.eps.ac ^ mgca.tau[i], 
                             MgCa_sw_m.eps.num / 
                               (1 - MgCa_sw_m.eps.ac ^ (2 * mgca.tau[i])))
    
    mgca.tau[i] = mgca.ages[i-1] - mgca.ages[i]   
  }
  
  MgCa_sw_m.eps.num = (1 - MgCa_sw_m.eps.ac ^ 2) * MgCa_sw_m.pre
  
  MgCa_sw_m.eps[1] ~ dnorm(0, MgCa_sw_m.pre)
  MgCa_sw_m[1] ~ dunif(MgCa_sw_m.init.min, MgCa_sw_m.init.max)
  MgCa_sw_m.init.min = 1
  MgCa_sw_m.init.max = 3
  
  #Priors on MgCa_sw model parameters  
  
  MgCa_sw_m.eps.ac ~ dunif(0.9, 1)
  
  MgCa_sw_m.pre ~ dgamma(MgCa_sw_m.pre.shp, MgCa_sw_m.pre.rate)
  MgCa_sw_m.pre.shp = 100
  MgCa_sw_m.pre.rate = 0.01
  
}