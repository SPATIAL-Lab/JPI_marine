model {

  #####
  ##JPI model for Mg/Ca and d18O proxy data, site 806
  #####
  
  #Data model for MgCa observations

  for(i in 1:length(MgCa)){
    MgCa[i] ~ dnorm(MgCa.m[i], MgCa_calib.pre)

    MgCa.m[i] = (a[1] + a[2] * BWT[MgCa.age.ind[i,1]]) * MgCa.sw[i] ^ a[3]

    MgCa.sw[i] ~ dnorm(MgCa_sw_m[MgCa.age.ind[i,2]], 1 / 0.03 ^ 2)

  }

  #Data model for MgCa_calib observations

  for(i in 1:length(MgCa_calib)){
    MgCa_calib[i] ~ dnorm(MgCa_calib.m[i], MgCa_calib.pre)

    MgCa_calib.m[i] = (a[1] + a[2] * MgCa_calib.bwt[i]) * MgCa_calib.sw[i] ^ a[3]

    MgCa_calib.bwt[i] ~ dnorm(MgCa_calib.bwt.m[i], 1 / MgCa_calib.bwt.sd[i] ^ 2)

    MgCa_calib.sw[i] ~ dnorm(MgCa_sw_m[MgCa.age.ind[i,2]], 1 / 0.03 ^ 2)

  }

  #Priors on MgCa_calib data model parameters

  MgCa_calib.pre ~ dgamma(MgCa_calib.pre.shp, MgCa_calib.pre.rate)
  MgCa_calib.pre.shp = 2
  MgCa_calib.pre.rate = 1/30

  a[1] ~ dnorm(a.1.m, 1 / a.1.var)
  a[2] ~ dnorm(a.2.m, 1 / a.2.var)
  a[3] ~ dnorm(a.3.m, 1 / a.3.var)

  a.1.m = 1.5
  a.1.var = 0.1 ^ 2
  a.2.m = 0.1
  a.2.var = 0.01 ^ 2
  a.3.m = -0.02 
  a.3.var = 0.03 ^ 2

  #Data model for d18O observations

  for(i in 1:length(d18O)){
    d18O[i] ~ dnorm(d18O.m[i], d18O.pre[i])

    d18O.pre[i] = ifelse(d18O.age.ind[i] < 345, d18O_calib.pre, d18O_calib.pre.2)    
    d18O.m[i] = d18O_sw[d18O.age.ind[i]] + b[1] + b[2] * BWT[d18O.age.ind[i]] + b[3] * BWT[d18O.age.ind[i]] ^ 2
  }

  d18O_calib.pre.2 ~ dgamma(d18O_calib.pre.2.shp, d18O_calib.pre.2.rate)
  d18O_calib.pre.2.shp = 6
  d18O_calib.pre.2.rate = 1
  
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

  b.1.m = 3.32
  b.1.var = 0.02 ^ 2
  b.2.m = -0.237
  b.2.var = 0.01 ^ 2
  b.3.m = 0.001
  b.3.var = 0.0005 ^ 2

  #Process model for BWT and d18O timeseries

  for(i in 2:nages){
    d18O_sw[i] = d18O_sw[i-1] + d18O_sw.eps[i] * tau[i]
    BWT[i] = BWT[i-1] + BWT.eps[i] * tau[i]
    
    d18O_sw.eps[i] ~ dnorm(exp(-(1 - d18O_sw.eps.ac) * tau[i]) * d18O_sw.eps[i - 1], 
                          1 / ((1 / d18O_sw.pre) / (2 * (1 - d18O_sw.eps.ac)) *
                             (1 - exp(-2 * (1 - d18O_sw.eps.ac) * tau[i]))))
    BWT.eps[i] ~ dnorm(exp(-(1 - BWT.eps.ac) * tau[i]) * BWT.eps[i - 1], 
                       1 / ((1 / BWT.pre) / (2 * (1 - BWT.eps.ac)) *
                         (1 - exp(-2 * (1 - BWT.eps.ac) * tau[i]))))
    
    tau[i] = ages[i - 1] - ages[i]

  }

  #Priors on BWT and d18O timeseries model parameters

  d18O_sw.eps[1] ~ dnorm(0, d18O_sw.pre)
  BWT.eps[1] ~ dnorm(0, BWT.pre) 
  d18O_sw[1] ~ dunif(d18O_sw.init.min, d18O_sw.init.max)
  BWT[1] ~ dunif(BWT.init.min, BWT.init.max)

  d18O_sw.init.min = -1
  d18O_sw.init.max = 1

  BWT.init.min = 3
  BWT.init.max = 8

  d18O_sw.eps.ac ~ dunif(0, 0.4)
  BWT.eps.ac ~ dunif(0, 0.4)

  d18O_sw.pre ~ dgamma(d18O_sw.pre.shp, d18O_sw.pre.rate)
  d18O_sw.pre.shp = 50
  d18O_sw.pre.rate = 0.5
  
  BWT.pre ~ dgamma(BWT.pre.shp, BWT.pre.rate)
  BWT.pre.shp = 25
  BWT.pre.rate = 0.5

  #Data model for seawater MgCa observations

  for(i in 1:length(MgCa_sw)){
    MgCa_sw[i] ~ dnorm(MgCa_sw_m[MgCa_sw.age.ind[i]], 1 / MgCa_sw.sd[i] ^ 2)
  
  }

  #Process model for MgCa_sw timeseries

  for(i in 2:nmgca.ages){
    MgCa_sw_m[i] = MgCa_sw_m[i-1] * ((MgCa_sw_m.eps[i] * mgca.tau[i]) + 1)
    
    MgCa_sw_m.eps[i] ~ dnorm(exp(-(1 - MgCa_sw_m.eps.ac) * mgca.tau[i]) * MgCa_sw_m.eps[i-1], 
                             1 / ((1 / MgCa_sw_m.pre) / (2 * (1 - MgCa_sw_m.eps.ac)) * 
                               (1 - exp(-2 * (1 - MgCa_sw_m.eps.ac) * mgca.tau[i]))))
    
    mgca.tau[i] = mgca.ages[i-1] - mgca.ages[i]   
  }
  
  #MgCa_sw_m.eps[1] ~ dunif(-MgCa_sw_m.eps.hr, MgCa_sw_m.eps.hr)
  MgCa_sw_m.eps[1] ~ dnorm(0, MgCa_sw_m.pre)
  MgCa_sw_m[1] ~ dunif(MgCa_sw_m.init.min, MgCa_sw_m.init.max)
  MgCa_sw_m.init.min = 1
  MgCa_sw_m.init.max = 3

  #Priors on MgCa_sw model parameters  
  
  MgCa_sw_m.eps.ac ~ dunif(0.90, 1)
  
  #MgCa_sw_m.eps.hr = 0.01
  
  MgCa_sw_m.pre ~ dgamma(MgCa_sw_m.pre.shp, MgCa_sw_m.pre.rate)
  MgCa_sw_m.pre.shp = 100
  MgCa_sw_m.pre.rate = 0.01

}

