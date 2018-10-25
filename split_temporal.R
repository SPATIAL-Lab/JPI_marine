split_AR = "model {

  for(i in 1:length(MgCa)){
    MgCa[i] ~ dnorm(MgCa.m[i], 1 / MgCa.var[i])

    MgCa.var[i] = (MgCa.m[i] * 0.01) ^ 2

    MgCa.m[i] = MgCa.m.e[i]
    
    MgCa.m.e[i] = e.1 * MgCa_sw[i] ^ e.2 * exp(e.3 * BWT[MgCa.age.ind[i]])
    
    MgCa_sw[i] ~ dunif(MgCa_sw.m[i] - MgCa_sw.r, MgCa_sw.m[i] + MgCa_sw.r)
    MgCa_sw.m[i] = 5.2 - 0.238 * Age[i] + 0.00661 * Age[i] ^ 2 - 6.66e-5 * Age[i] ^ 3

  }

  MgCa_sw.r = 0.5

  e.1 ~ dnorm(e.1.m, 1 / e.1.var)
  e.2 ~ dnorm(e.2.m, 1 / e.2.var)
  e.3 ~ dnorm(e.3.m, 1 / e.3.var)
  e.1.m = 0.66
  e.1.var = 0.04 ^ 2
  e.2.m = 0.27
  e.2.var = 0.03 ^ 2
  e.3.m = 0.114
  e.3.var = 0.01 ^ 2

  for(i in 1:length(d18O)){
    d18O[i] ~ dnorm(d18O.m[i], 1 / d18O.var)
    
    d18O.m[i] = d18O_sw[d18O.age.ind[i]] + a.1[i] + a.2[i] * BWT[d18O.age.ind[i]] + a.3[i] * BWT[d18O.age.ind[i]] ^ 2
    a.1[i] ~ dnorm(a.1.m, 1 / a.1.var)
    a.2[i] ~ dnorm(a.2.m, 1 / a.2.var)
    a.3[i] ~ dnorm(a.3.m, 1 / a.3.var)
  }

  a.1.m = 3.31
  a.1.var = 0.02 ^ 2
  a.2.m = -0.245
  a.2.var = 0.005 ^ 2
  a.3.m = 0.0011
  a.3.var = 0.0002 ^ 2

  d18O.var = 0.1 ^ 2

  for(i in 2:nages){
    d18O_sw[i] = d18O_sw[i-1] + d18O_sw.eps[i]
    BWT[i] = BWT[i-1] + BWT.eps[i]
    
    d18O_sw.eps[i] ~ dnorm(d18O_sw.eps[i - 1] * d18O_sw.eps.ac, 1 / d18O_sw.var)
    BWT.eps[i] ~ dnorm(BWT.eps[i - 1] * BWT.eps.ac, 1 / BWT.var)

  }

  d18O_sw.eps[1] ~ dnorm(0, 1 / d18O_sw.var)
  BWT.eps[1] ~ dnorm(0, 1 / BWT.var) 
  d18O_sw[1] = d18O_sw.init
  BWT[1] = BWT.init

  d18O_sw.init ~ dunif(d18O_sw.init.min, d18O_sw.init.max)
  d18O_sw.init.min = -0.5
  d18O_sw.init.max = 1

  BWT.init ~ dunif(BWT.init.min, BWT.init.max)
  BWT.init.min = 5
  BWT.init.max = 9

  d18O_sw.eps.ac ~ dunif(0, 1)
  BWT.eps.ac ~ dunif(0, 1)

  d18O_sw.var ~ dgamma(d18O_sw.var.k, 1 / d18O_sw.var.theta)
  d18O_sw.var.k = d18O_sw.var.m / d18O_sw.var.theta
  d18O_sw.var.theta = d18O_sw.var.var / d18O_sw.var.m
  d18O_sw.var.m = 0.05
  d18O_sw.var.var = 0.1 ^ 2
  
  BWT.var ~ dgamma(BWT.var.k, 1 / BWT.var.theta)
  BWT.var.k = BWT.var.m / BWT.var.theta
  BWT.var.theta = BWT.var.var / BWT.var.m
  BWT.var.m = 0.1
  BWT.var.var = 0.2 ^ 2
}
"