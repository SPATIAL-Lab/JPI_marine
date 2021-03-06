model {
  
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
  
  #Precision based on Uvigerina coretop variance
  MgCa_calib.pre ~ dgamma(MgCa_calib.pre.shp, MgCa_calib.pre.rate)
  MgCa_calib.pre.shp = 2
  MgCa_calib.pre.rate = 1/30
  
  #Using very loose constraints on MgCa calib parameters
  a[1] ~ dunif(0.9, 1.1)
  a[2] ~ dunif(-0.05, 0.35)
  a[3] ~ dunif(-0.02, 0.02) 
  
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
  
}