model {

  #Data model for MgCa observations

  for(i in 1:length(MgCa.b)){
    MgCa.b[i] ~ dnorm(MgCa.b.m[i], MgCa_calib.pre)

    MgCa.b.m[i] = (a[1] + a[2] * BWT.b[MgCa.age.ind.b[i]]) * MgCa.sw.b[i] ^ a[3]

    MgCa.sw.b[i] ~ dnorm(MgCa.sw.m, 1 / 0.03 ^ 2)
    
  }
  
  for(i in 1:length(MgCa.e)){
    MgCa.e[i] ~ dnorm(MgCa.e.m[i], MgCa_calib.pre)
    
    MgCa.e.m[i] = (a[1] + a[2] * BWT.e[MgCa.age.ind.e[i]]) * MgCa.sw.e[i] ^ a[3]
    
    MgCa.sw.e[i] ~ dnorm(MgCa.sw.m, 1 / 0.03 ^ 2)
    
  }

  MgCa.sw.m ~ dnorm(MgCa_sw.neo[1], 1 / MgCa_sw.neo[2] ^ 2)
  
  #Data model for MgCa_calib observations

  for(i in 1:length(MgCa_calib)){
    MgCa_calib[i] ~ dnorm(MgCa_calib.m[i], MgCa_calib.pre)

    MgCa_calib.m[i] = (a[1] + a[2] * MgCa_calib.bwt[i]) * MgCa_calib.sw[i] ^ a[3]

    MgCa_calib.bwt[i] ~ dnorm(MgCa_calib.bwt.m[i], 1 / MgCa_calib.bwt.sd[i] ^ 2)

    MgCa_calib.sw[i] ~ dnorm(5.2, 1 / 0.03 ^ 2)

  }

  #Priors on MgCa_calib data model parameters

  MgCa_calib.pre ~ dgamma(MgCa_calib.pre.shp, MgCa_calib.pre.rate)
  MgCa_calib.pre.shp = 2
  MgCa_calib.pre.rate = 1/30

  a[1] ~ dnorm(a.1.m, 1 / a.1.var)
  a[2] ~ dnorm(a.2.m, 1 / a.2.var)
  a[3] ~ dnorm(a.3.m, 1 / a.3.var)

  a.1.m = 1.02
  a.1.var = 0.1 ^ 2
  a.2.m = 0.07
  a.2.var = 0.01 ^ 2
  a.3.m = -0.02 
  a.3.var = 0.03 ^ 2

  #Data model for d18O observations

  for(i in 1:length(d18O.b)){
    d18O.b[i] ~ dnorm(d18O.b.m[i], d18O_calib.c.pre)

    d18O.b.m[i] = d18O_sw.b[d18O.age.ind.b[i]] + b.c[1] + b.c[2] * BWT.b[d18O.age.ind.b[i]] + b.c[3] * BWT.b[d18O.age.ind.b[i]] ^ 2
  }

  for(i in 1:length(d18O.e)){
    d18O.e[i] ~ dnorm(d18O.e.m[i], d18O_calib.u.pre)
    
    d18O.e.m[i] = d18O_sw.e[d18O.age.ind.e[i]] + b.u[1] + b.u[2] * BWT.e[d18O.age.ind.e[i]] + b.u[3] * BWT.e[d18O.age.ind.e[i]] ^ 2
  }
  
  #Data model for d18O_calib observations - Uvigerina

  for(i in 1:length(d18O_calib.u)){
    d18O_calib.u[i] ~ dnorm(d18O_calib.u.m[i], d18O_calib.u.pre)

    d18O_calib.u.m[i] = b.u[1] + b.u[2] * d18O_calib.u.bwt[i] + b.u[3] * d18O_calib.u.bwt[i] ^ 2

    d18O_calib.u.bwt[i] ~ dnorm(d18O_calib.u.bwt.m[i], 1 / d18O_calib.u.bwt.sd[i])
  }

  # Priors on d18O data model parameters - Uvigerina

  d18O_calib.u.pre ~ dgamma(d18O_calib.u.pre.shp, d18O_calib.u.pre.rate)
  d18O_calib.u.pre.shp = 3
  d18O_calib.u.pre.rate = 1/30

  b.u[1] ~ dnorm(b.u.1.m, 1 / b.u.1.var)
  b.u[2] ~ dnorm(b.u.2.m, 1 / b.u.2.var)
  b.u[3] ~ dnorm(b.u.3.m, 1 / b.u.3.var)

  b.u.1.m = 4.05
  b.u.1.var = 0.06 ^ 2
  b.u.2.m = -0.215
  b.u.2.var = 0.02 ^ 2
  b.u.3.m = -0.001
  b.u.3.var = 0.001 ^ 2

  #Data model for d18O_calib observations - Cib
  
  for(i in 1:length(d18O_calib.c)){
    d18O_calib.c[i] ~ dnorm(d18O_calib.c.m[i], d18O_calib.c.pre)
    
    d18O_calib.c.m[i] = b.c[1] + b.c[2] * d18O_calib.c.bwt[i] + b.c[3] * d18O_calib.c.bwt[i] ^ 2
    
    d18O_calib.c.bwt[i] ~ dnorm(d18O_calib.c.bwt.m[i], 1 / d18O_calib.c.bwt.sd[i])
  }
  
  # Priors on d18O data model parameters - Cib
  
  d18O_calib.c.pre ~ dgamma(d18O_calib.c.pre.shp, d18O_calib.c.pre.rate)
  d18O_calib.c.pre.shp = 3
  d18O_calib.c.pre.rate = 1/30
  
  b.c[1] ~ dnorm(b.c.1.m, 1 / b.c.1.var)
  b.c[2] ~ dnorm(b.c.2.m, 1 / b.c.2.var)
  b.c[3] ~ dnorm(b.c.3.m, 1 / b.c.3.var)
  
  b.c.1.m = 3.32
  b.c.1.var = 0.02 ^ 2
  b.c.2.m = -0.237
  b.c.2.var = 0.01 ^ 2
  b.c.3.m = 0.001
  b.c.3.var = 0.0005 ^ 2

  #Process model for BWT and d18O timeseries

  for(i in 2:nages){
    d18O_sw.b[i] = d18O_sw.b[i-1] + d18O_sw.b.eps[i]
    BWT.b[i] = BWT.b[i-1] + BWT.b.eps[i]
    
    d18O_sw.e[i] = d18O_sw.e[i-1] + d18O_sw.e.eps[i]
    BWT.e[i] = BWT.e[i-1] + BWT.e.eps[i]
    
    d18O_sw.b.eps[i] ~ dnorm(exp(-(1 - d18O_sw.b.eps.ac) * tau[i]) * d18O_sw.b.eps[i - 1], 
                           1 / ((1 / d18O_sw.b.pre) / (2 * (1 - d18O_sw.b.eps.ac)) *
                                  (1 - exp(-2 * (1 - d18O_sw.b.eps.ac) * tau[i]))))
    BWT.b.eps[i] ~ dnorm(exp(-(1 - BWT.b.eps.ac) * tau[i]) * BWT.b.eps[i - 1], 
                       1 / ((1 / BWT.b.pre) / (2 * (1 - BWT.b.eps.ac)) *
                              (1 - exp(-2 * (1 - BWT.b.eps.ac) * tau[i]))))

    
    d18O_sw.e.eps[i] ~ dnorm(exp(-(1 - d18O_sw.e.eps.ac) * tau[i]) * d18O_sw.e.eps[i - 1], 
                             1 / ((1 / d18O_sw.e.pre) / (2 * (1 - d18O_sw.e.eps.ac)) *
                                    (1 - exp(-2 * (1 - d18O_sw.e.eps.ac) * tau[i]))))
    BWT.e.eps[i] ~ dnorm(exp(-(1 - BWT.e.eps.ac) * tau[i]) * BWT.e.eps[i - 1], 
                         1 / ((1 / BWT.e.pre) / (2 * (1 - BWT.e.eps.ac)) *
                                (1 - exp(-2 * (1 - BWT.e.eps.ac) * tau[i]))))    

    tau[i] = ages[i - 1] - ages[i]
    
  }

  #Priors on BWT and d18O timeseries model parameters

  d18O_sw.b.eps[1] ~ dnorm(0, d18O_sw.b.pre)
  BWT.b.eps[1] ~ dnorm(0, BWT.b.pre) 
  d18O_sw.b[1] ~ dunif(0, 1.2)
  BWT.b[1] ~ dunif(2, 6)
  
  d18O_sw.e.eps[1] ~ dnorm(0, d18O_sw.e.pre)
  BWT.e.eps[1] ~ dnorm(0, BWT.e.pre) 
  d18O_sw.e[1] ~ dunif(-0.5, 0.5)
  BWT.e[1] ~ dunif(-1, 3)

  d18O_sw.b.eps.ac ~ dunif(0, 0.8)
  BWT.b.eps.ac ~ dunif(0, 0.8)

  d18O_sw.e.eps.ac ~ dunif(0, 0.8)
  BWT.e.eps.ac ~ dunif(0, 0.8)

  d18O_sw.b.pre ~ dgamma(d18O_sw.pre.shp, d18O_sw.pre.rate)
  d18O_sw.e.pre ~ dgamma(d18O_sw.pre.shp, d18O_sw.pre.rate)
  d18O_sw.pre.shp = 10
  d18O_sw.pre.rate = 1/5
  
  BWT.b.pre ~ dgamma(BWT.pre.shp, BWT.pre.rate)
  BWT.e.pre ~ dgamma(BWT.pre.shp, BWT.pre.rate)
  BWT.pre.shp = 20
  BWT.pre.rate = 2

}

