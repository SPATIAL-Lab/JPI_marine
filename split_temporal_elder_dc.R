model {
  
  #####
  ##JPI model for Mg/Ca and d18O proxy data, site U1385, using down-core calibration
  #####
  
  #Data model for MgCa observations
  
  for(i in 1:length(MgCa)){
    MgCa[i] ~ dnorm(MgCa.m[i], MgCa_calib.pre)
    
    MgCa.m[i] = (a[1] + a[2] * BWT[MgCa.age.ind[i]]) * MgCa_sw[i] ^ a[3]
    
    MgCa_sw[i] ~ dnorm(MgCa_sw.m, 1 / 0.03 ^ 2)
    
  }
  
  MgCa_sw.m ~ dnorm(MgCa_sw.neo[1], 1 / MgCa_sw.neo[2] ^ 2)
  
  #Data model for MgCa_calib observations
  
  for(i in 1:length(MgCa_calib)){
    MgCa_calib[i] ~ dnorm(MgCa_calib.m[i], MgCa_calib.pre)
    
    MgCa_calib.m[i] = (a[1] + a[2] * MgCa_calib.bwt[i]) * MgCa_calib.sw[i] ^ a[3]
    
    MgCa_calib.bwt[i] ~ dnorm(MgCa_calib.bwt.m[i], 1 / MgCa_calib.bwt.sd[i] ^ 2)
    
    MgCa_calib.sw[i] ~ dnorm(5.2, 1 / 0.03 ^ 2)
    
  }
  
  ##Downcore Mg/Ca calibration
  #Data model for MgCa
  for(i in 1:length(MgCa_HOL)){
    MgCa_HOL[i] ~ dnorm(MgCa_HOL.m[i], MgCa_calib.pre)
    MgCa_HOL.m[i] = (a[1] + a[2] * BWT_HOL) * MgCa_sw_HOL ^ a[3]
  }
  
  for(i in 1:length(MgCa_LGM)){
    MgCa_LGM[i] ~ dnorm(MgCa_LGM.m[i], MgCa_calib.pre)
    MgCa_LGM.m[i] = (a[1] + a[2] * BWT_LGM) * MgCa_sw_LGM ^ a[3]
  }
  
  #Data model for d18O
  for(i in 1:length(d18O_HOL)){
    d18O_HOL[i] ~ dnorm(d18O_HOL.m[i], d18O_calib.pre)
    d18O_HOL.m[i] = d18O_sw_HOL + b[1] + b[2] * BWT_HOL + b[3] * BWT_HOL ^ 2
  }
  
  for(i in 1:length(d18O_LGM)){
    d18O_LGM[i] ~ dnorm(d18O_LGM.m[i], d18O_calib.pre)
    d18O_LGM.m[i] = d18O_sw_LGM + b[1] + b[2] * BWT_LGM + b[3] * BWT_LGM ^ 2
  }
  
  #Seawater d18O for both states
  d18O_sw_LGM = d18O_sw_HOL + D_d18O_sw_LGM
  D_d18O_sw_LGM ~ dnorm(1.1, 1 / 0.1 ^ 2) #Estimate from Adkins et al 2002
  d18O_sw_HOL ~ dnorm(-0.1, 1 / 0.2 ^ 2) #Estimated from Elderfield et al 2010, fig 6
  
  #Seawater Mg/Ca for both states is ~modern
  MgCa_sw_LGM ~ dnorm(5.2, 1 / 0.02 ^ 2)
  MgCa_sw_HOL ~ dnorm(5.2, 1 / 0.02 ^ 2)
  
  #Bottom water temperatures for both states
  BWT_LGM = BWT_HOL + D_BWT_LGM
  D_BWT_LGM ~ dunif(-5, 0)
  BWT_HOL ~ dnorm(3.7, 1 / 0.4 ^ 2) #Estimated from Elderfield et al 2010, fig 6
  
  #Priors on MgCa_calib data model parameters
  
  MgCa_calib.pre ~ dgamma(MgCa_calib.pre.shp, MgCa_calib.pre.rate)
  MgCa_calib.pre.shp = 2
  MgCa_calib.pre.rate = 1/30
  
  a[1] ~ dnorm(a.1.m, 1 / a.1.var)
  a[2] ~ dnorm(a.2.m, 1 / a.2.var)
  a[3] ~ dnorm(a.3.m, 1 / a.3.var)
  
  a.1.m = 0.985
  a.1.var = 0.1 ^ 2
  a.2.m = 0.08
  a.2.var = 0.05 ^ 2
  a.3.m = -0.001 
  a.3.var = 0.02 ^ 2
  
  #Data model for d18O observations
  
  for(i in 1:length(d18O)){
    d18O[i] ~ dnorm(d18O.m[i], d18O_calib.pre)
    
    d18O.m[i] = d18O_sw[d18O.age.ind[i]] + b[1] + b[2] * BWT[d18O.age.ind[i]] + b[3] * BWT[d18O.age.ind[i]] ^ 2
  }
  
  #Data model for d18O_calib observations
  
  for(i in 1:length(d18O_calib)){
    d18O_calib[i] ~ dnorm(d18O_calib.m[i], d18O_calib.pre)
    
    d18O_calib.m[i] = b[1] + b[2] * d18O_calib.bwt[i] + b[3] * d18O_calib.bwt[i] ^ 2
    
    d18O_calib.bwt[i] ~ dnorm(d18O_calib.bwt.m[i], 1 / d18O_calib.bwt.sd[i])
  }
  
  # Priors on d18O data model parameters
  
  d18O_calib.pre ~ dgamma(d18O_calib.pre.shp, d18O_calib.pre.rate)
  d18O_calib.pre.shp = 3
  d18O_calib.pre.rate = 1/30
  
  b[1] ~ dnorm(b.1.m, 1 / b.1.var)
  b[2] ~ dnorm(b.2.m, 1 / b.2.var)
  b[3] ~ dnorm(b.3.m, 1 / b.3.var)
  
  b.1.m = 4.05
  b.1.var = 0.06 ^ 2
  b.2.m = -0.215
  b.2.var = 0.02 ^ 2
  b.3.m = -0.001
  b.3.var = 0.001 ^ 2
  
  #Process model for BWT and d18O timeseries
  
  for(i in 2:nages){
    d18O_sw[i] = d18O_sw[i-1] + d18O_sw.eps[i]
    BWT[i] = BWT[i-1] + BWT.eps[i]
    
    d18O_sw.eps[i] ~ dnorm(d18O_sw.eps[i - 1] * d18O_sw.eps.ac ^ tau[i], 
                           d18O_sw.eps.num / (1 - d18O_sw.eps.ac ^ (2*tau[i])))
    BWT.eps[i] ~ dnorm(BWT.eps[i - 1] * BWT.eps.ac ^ tau[i], 
                       BWT.eps.num / (1 - BWT.eps.ac ^ (2 * tau[i])))
    
    tau[i] = (ages[i - 1] - ages[i])
    
  }
  
  #Calculate some values that are constant for all time steps
  
  d18O_sw.eps.num = (1 - d18O_sw.eps.ac ^ 2) * d18O_sw.pre
  BWT.eps.num = (1 - BWT.eps.ac ^ 2) * BWT.pre
  
  #Priors on BWT and d18O timeseries model parameters
  
  d18O_sw.eps[1] ~ dnorm(0, d18O_sw.pre)
  BWT.eps[1] ~ dnorm(0, BWT.pre) 
  d18O_sw[1] ~ dunif(d18O_sw.init.min, d18O_sw.init.max)
  BWT[1] ~ dunif(BWT.init.min, BWT.init.max)
  
  d18O_sw.init.min = -0.5
  d18O_sw.init.max = 0.5
  
  BWT.init.min = -1
  BWT.init.max = 4
  
  d18O_sw.eps.ac ~ dunif(0, 0.8)
  BWT.eps.ac ~ dunif(0, 0.8)
  
  d18O_sw.pre ~ dgamma(d18O_sw.pre.shp, d18O_sw.pre.rate)
  d18O_sw.pre.shp = 20
  d18O_sw.pre.rate = 0.1
  
  BWT.pre ~ dgamma(BWT.pre.shp, BWT.pre.rate)
  BWT.pre.shp = 20
  BWT.pre.rate = 1
  
}

