lear_AR = "model {

  for(i in 1:nobs){
    MgCa[i] ~ dnorm(MgCa.m[i], 1 / MgCa.var[i])
    d18O[i] ~ dnorm(d18O.m[i], 1 / d18O.var)

    MgCa.var[i] = (MgCa.m[i] * 0.01) ^ 2

    MgCa.m[i] = MgCa.m.e[i]
    
    MgCa.m.l[i] = l.1[i] + l.2[i] * BWT[age.ind[i]] * MgCa_sw[i] ^ l.3[i]
    l.1[i] ~ dnorm(l.1.m, 1 / l.1.var)
    l.2[i] ~ dnorm(l.2.m, 1 / l.2.var)
    l.3[i] ~ dnorm(l.3.m, 1 / l.3.var)
    
    MgCa.m.e[i] = e.1[i] * MgCa_sw[i] ^ e.2[i] * exp(e.3[i] * BWT[age.ind[i]])
    e.1[i] ~ dnorm(e.1.m, 1 / e.1.var)
    e.2[i] ~ dnorm(e.2.m, 1 / e.2.var)
    e.3[i] ~ dnorm(e.3.m, 1 / e.3.var)
    
    MgCa_sw[i] ~ dunif(MgCa_sw.m[i] - MgCa_sw.r, MgCa_sw.m[i] + MgCa_sw.r)
    MgCa_sw.m[i] = 5.2 - 0.238 * Age[i] + 0.00661 * Age[i] ^ 2 - 6.66e-5 * Age[i] ^ 3
    
    d18O.m[i] = d18O_sw[age.ind[i]] + a.1[i] + a.2[i] * BWT[age.ind[i]] + a.3[i] * BWT[age.ind[i]] ^ 2
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

  MgCa_sw.r = 0.5

  e.1.m = 0.66
  e.1.var = 0.04 ^ 2
  e.2.m = 0.27
  e.2.var = 0.03 ^ 2
  e.3.m = 0.114
  e.3.var = 0.01 ^ 2

  l.1.m = 1.21
  l.1.var = 0.02 ^ 2
  l.2.m = 0.12
  l.2.var = 0.002 ^ 2
  l.3.m = -0.003
  l.3.var = 0.01 ^ 2

  l.prob = 1/3

  d18O.var = 0.04 ^ 2


  for(i in 2:nages){
    d18O_sw[i] = d18O_sw.init + d18O_sw.ac * (d18O_sw[i-1] - d18O_sw.init) + d18O_sw.eps[i]
    BWT[i] = BWT.init + BWT.ac * (BWT[i-1] - BWT.init) + BWT.eps[i]

    d18O_sw.eps[i] ~ dnorm(d18O_sw.eps.ac * d18O_sw.eps[i-1], 1 / d18O_sw.var)
    BWT.eps[i] ~ dnorm(BWT.eps.ac * BWT.eps[i-1], 1 / BWT.var)
  }

  d18O_sw.eps[1] ~ dnorm(0, 1 / d18O_sw.var)
  BWT.eps[1] ~ dnorm(0, 1 / BWT.var)
  d18O_sw[1] = d18O_sw.init
  BWT[1] = BWT.init

  d18O_sw.eps.ac ~ dbeta(d18O_sw.eps.ac.alpha, d18O_sw.eps.ac.beta)
  d18O_sw.eps.ac.alpha = d18O_sw.eps.ac.m * d18O_sw.eps.ac.size
  d18O_sw.eps.ac.beta = (1 - d18O_sw.eps.ac.m) * d18O_sw.eps.ac.size
  d18O_sw.eps.ac.size = d18O_sw.eps.ac.m * (1 - d18O_sw.eps.ac.m) / d18O_sw.eps.ac.var - 1
  d18O_sw.eps.ac.var = 0.2 ^ 2
  d18O_sw.eps.ac.m = 0.5

  d18O_sw.ac ~ dbeta(d18O_sw.ac.alpha, d18O_sw.ac.beta)
  d18O_sw.ac.alpha = d18O_sw.ac.m * d18O_sw.ac.size
  d18O_sw.ac.beta = (1 - d18O_sw.ac.m) * d18O_sw.ac.size
  d18O_sw.ac.size = d18O_sw.ac.m * (1 - d18O_sw.ac.m) / d18O_sw.ac.var - 1
  d18O_sw.ac.var = 0.04 ^ 2
  d18O_sw.ac.m = 0.88

  BWT.eps.ac ~ dbeta(BWT.eps.ac.alpha, BWT.eps.ac.beta)
  BWT.eps.ac.alpha = BWT.eps.ac.m * BWT.eps.ac.size
  BWT.eps.ac.beta = (1 - BWT.eps.ac.m) * BWT.eps.ac.size
  BWT.eps.ac.size = BWT.eps.ac.m * (1 - BWT.eps.ac.m) / BWT.eps.ac.var - 1
  BWT.eps.ac.var = 0.2 ^ 2
  BWT.eps.ac.m = 0.5

  BWT.ac ~ dbeta(BWT.ac.alpha, BWT.ac.beta)
  BWT.ac.alpha = BWT.ac.m * BWT.ac.size
  BWT.ac.beta = (1 - BWT.ac.m) * BWT.ac.size
  BWT.ac.size = BWT.ac.m * (1 - BWT.ac.m) / BWT.ac.var - 1
  BWT.ac.var = 0.04 ^ 2
  BWT.ac.m = 0.88

  d18O_sw.init ~ dunif(d18O_sw.init.min, d18O_sw.init.max)
  d18O_sw.init.min = -2
  d18O_sw.init.max = 5

  BWT.init ~ dunif(BWT.init.min, BWT.init.max)
  BWT.init.min = 0
  BWT.init.max = 10

  d18O_sw.var ~ dnorm(d18O_sw.var.m, 1 / d18O_sw.var.var) I (0, )
  d18O_sw.var.m = 0.5 ^ 2
  d18O_sw.var.var = 0.1 ^ 4
  
  BWT.var ~ dnorm(BWT.var.m, 1 / BWT.var.var) I (0, )
  BWT.var.m = 2 ^ 2
  BWT.var.var = 0.4 ^ 4
}
"