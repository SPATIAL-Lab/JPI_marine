model {

  #Data model for MgCa observations

  for(i in 1:length(MgCa)){
    MgCa[i] ~ dnorm(MgCa.m[i], MgCa_calib.pre)

    MgCa.m[i] = (lc[1] + lc[2] * BWT[MgCa.age.ind[i,1]]) * MgCa.sw[i] ^ lc[3]
    #MgCa.m[i] = ec[1] * MgCa.sw[i] ^ ec[2] * exp(ec[3] * BWT[MgCa.age.ind[i,1]])

    MgCa.sw[i] ~ dnorm(MgCa_sw_m[MgCa.age.ind[i,2]], 1 / 0.03 ^ 2)

  }

  #Data model for MgCa_calib observations

  for(i in 1:length(MgCa_calib)){
    MgCa_calib[i] ~ dnorm(MgCa_calib.m[i], MgCa_calib.pre)

    MgCa_calib.m[i] = (lc[1] + lc[2] * MgCa_calib.bwt[i]) * MgCa_calib.sw[i] ^ lc[3]
    #MgCa_calib.m[i] = ec[1] * MgCa_calib.sw[i] ^ ec[2] * exp(ec[3] * MgCa_calib.bwt[i])

    MgCa_calib.bwt[i] ~ dnorm(MgCa_calib.bwt.m[i], 1 / MgCa_calib.bwt.sd[i] ^ 2)

    MgCa_calib.sw[i] ~ dnorm(MgCa_sw_m[MgCa.age.ind[i,2]], 1 / 0.03 ^ 2)

  }

  #Priors on MgCa_calib data model parameters

  MgCa_calib.pre ~ dgamma(MgCa_calib.pre.shp, MgCa_calib.pre.rate)
  MgCa_calib.pre.shp = 2
  MgCa_calib.pre.rate = 1/30

  #ec[1] ~ dnorm(ec.1.m, 1 / ec.1.var)
  #ec[2] ~ dnorm(ec.2.m, 1 / ec.2.var)
  #ec[3] ~ dnorm(ec.3.m, 1 / ec.3.var)

  lc[1] ~ dnorm(lc.1.m, 1 / lc.1.var)
  lc[2] ~ dnorm(lc.2.m, 1 / lc.2.var)
  lc[3] ~ dnorm(lc.3.m, 1 / lc.3.var)

  ec.1.m = 0.7
  ec.1.var = 0.04 ^ 2
  ec.2.m = 0.4
  ec.2.var = 0.03 ^ 2
  ec.3.m = 0.1
  ec.3.var = 0.01 ^ 2

  lc.1.m = 1.5
  lc.1.var = 0.28 ^ 2
  lc.2.m = 0.1
  lc.2.var = 0.03 ^ 2
  lc.3.m = -0.02 
  lc.3.var = 0.047 ^ 2

  #Data model for d18O observations

  for(i in 1:length(d18O)){
    d18O[i] ~ dnorm(d18O.m[i], d18O.pre[i])

    d18O.pre[i] = ifelse(d18O.age.ind[i] < 345, d18O_calib.pre, 1 / (1 / d18O_calib.pre + (d18O_calib.sd.off ^ 2)))    
    d18O.m[i] = d18O_sw[d18O.age.ind[i]] + a[1] + a[2] * BWT[d18O.age.ind[i]] + a[3] * BWT[d18O.age.ind[i]] ^ 2
  }

  d18O_calib.sd.off ~ dnorm(0.25, 1 / 0.05 ^ 2)
  
  #Data model for d18O_calib observations

  for(i in 1:length(d18O_calib)){
    d18O_calib[i] ~ dnorm(d18O_calib.m[i], d18O_calib.pre)

    d18O_calib.m[i] = a[1] + a[2] * d18O_calib.bwt[i] + a[3] * d18O_calib.bwt[i] ^ 2

    d18O_calib.bwt[i] ~ dnorm(d18O_calib.bwt.m[i], 1 / d18O_calib.bwt.sd[i])
  }

  # Priors on d18O data model parameters

  d18O_calib.pre ~ dgamma(d18O_calib.pre.shp, d18O_calib.pre.rate)
  d18O_calib.pre.shp = 1
  d18O_calib.pre.rate = 1/100

  a[1] ~ dnorm(a.1.m, 1 / a.1.var)
  a[2] ~ dnorm(a.2.m, 1 / a.2.var)
  a[3] ~ dnorm(a.3.m, 1 / a.3.var)

  a.1.m = 3.32
  a.1.var = 0.02 ^ 2
  a.2.m = -0.237
  a.2.var = 0.01 ^ 2
  a.3.m = 0.0005
  a.3.var = 0.0005 ^ 2

  #System model for BWT and d18O timeseries

  for(i in 2:nages){
    d18O_sw[i] = d18O_sw[i-1] + d18O_sw.eps[i]
    BWT[i] = BWT[i-1] + BWT.eps[i]
    
    d18O_sw.eps[i] ~ dnorm(d18O_sw.eps[i - 1] * d18O_sw.eps.ac, d18O_sw.pre)
    BWT.eps[i] ~ dnorm(BWT.eps[i - 1] * BWT.eps.ac, BWT.pre)

  }

  #Priors on BWT and d18O timeseries model parameters

  d18O_sw.eps[1] ~ dnorm(0, d18O_sw.pre)
  BWT.eps[1] ~ dnorm(0, BWT.pre) 
  d18O_sw[1] = d18O_sw.init
  BWT[1] = BWT.init

  d18O_sw.init ~ dunif(d18O_sw.init.min, d18O_sw.init.max)
  d18O_sw.init.min = -1
  d18O_sw.init.max = 1

  BWT.init ~ dunif(BWT.init.min, BWT.init.max)
  BWT.init.min = 3
  BWT.init.max = 8

  d18O_sw.eps.ac ~ dunif(0, 0.4)
  BWT.eps.ac ~ dunif(0, 0.4)

  d18O_sw.pre ~ dgamma(d18O_sw.pre.shp, d18O_sw.pre.rate)
  d18O_sw.pre.shp = 10
  d18O_sw.pre.rate = 1/5
  
  BWT.pre ~ dgamma(BWT.pre.shp, BWT.pre.rate)
  BWT.pre.shp = 20
  BWT.pre.rate = 2

  #Data model for seawater MgCa observations

  for(i in 1:length(MgCa_sw)){
    MgCa_sw[i] ~ dnorm(MgCa_sw_m[MgCa_sw.age.ind[i]], 1 / MgCa_sw.sd[i] ^ 2)
  
  }

  #System model for MgCa_sw timeseries

  for(i in 2:nmgca.ages){
    MgCa_sw_m[i] = MgCa_sw_m[i-1] * (MgCa_sw_m.eps[i] + 1)
    
    #MgCa_sw_m.eps[i] ~ dunif(MgCa_sw_m.eps.co[i] - MgCa_sw_m.eps.hr, MgCa_sw_m.eps.co[i] + MgCa_sw_m.eps.hr)
    MgCa_sw_m.eps[i] ~ dnorm(MgCa_sw_m.eps[i - 1] * MgCa_sw_m.eps.ac, MgCa_sw_m.pre)
    MgCa_sw_m.eps.co[i] = MgCa_sw_m.eps[i - 1] * MgCa_sw_m.eps.ac
  }
  
  #MgCa_sw_m.eps[1] ~ dunif(-MgCa_sw_m.eps.hr, MgCa_sw_m.eps.hr)
  MgCa_sw_m.eps[1] ~ dnorm(0, MgCa_sw_m.pre)
  MgCa_sw_m[1] ~ dunif(MgCa_sw_m.init.min, MgCa_sw_m.init.max)
  MgCa_sw_m.init.min = 1
  MgCa_sw_m.init.max = 2

  #Priors on MgCa_sw model parameters  
  
  MgCa_sw_m.eps.ac ~ dunif(0.90, 1)
  
  #MgCa_sw_m.eps.hr = 0.01
  
  MgCa_sw_m.pre ~ dgamma(MgCa_sw_m.pre.shp, MgCa_sw_m.pre.rate)
  MgCa_sw_m.pre.shp = 1000
  MgCa_sw_m.pre.rate = 0.01

}

