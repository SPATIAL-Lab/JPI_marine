lear = "model {

  MgCa ~ dnorm(MgCa.m, 1 / MgCa.var)
  d18O ~ dnorm(d18O.m, 1 / d18O.var)

  MgCa.var = (MgCa.m * 0.01) ^ 2
  d18O.var = 0.04 ^ 2

  MgCa.m = ifelse(le == 1, MgCa.m.l, MgCa.m.e) 
  le ~ dbern(l.prob)
  l.prob = 1/3
  
  MgCa.m.l = l.1 + l.2 * BWT * MgCa.sw ^ l.3
  l.1 ~ dnorm(l.1.m, 1 / l.1.var)
  l.2 ~ dnorm(l.2.m, 1 / l.2.var)
  l.3 ~ dnorm(l.3.m, 1 / l.3.var)
  l.1.m = 1.21
  l.1.var = 0.02 ^ 2
  l.2.m = 0.12
  l.2.var = 0.002 ^ 2
  l.3.m = -0.003
  l.3.var = 0.01 ^ 2
  
  MgCa.m.e = e.1 * MgCa.sw ^ e.2 * exp(e.3 * BWT)
  e.1 ~ dnorm(e.1.m, 1 / e.1.var)
  e.2 ~ dnorm(e.2.m, 1 / e.2.var)
  e.3 ~ dnorm(e.3.m, 1 / e.3.var)
  e.1.m = 0.66
  e.1.var = 0.04 ^ 2
  e.2.m = 0.27
  e.2.var = 0.03 ^ 2
  e.3.m = 0.114
  e.3.var = 0.01 ^ 2
  
  MgCa.sw ~ dunif(MgCa.sw.m - MgCa.sw.r, MgCa.sw.m + MgCa.sw.r)
  MgCa.sw.m = 5.2 - 0.238 * Age + 0.00661 * Age ^ 2 - 6.66e-5 * Age ^ 3
  MgCa.sw.r = 0.5
  
  d18O.m = d18O.sw + a.1 + a.2 * BWT + a.3 * BWT ^ 2
  a.1 ~ dnorm(a.1.m, 1 / a.1.var)
  a.2 ~ dnorm(a.2.m, 1 / a.2.var)
  a.3 ~ dnorm(a.3.m, 1 / a.3.var)
  a.1.m = 3.31
  a.1.var = 0.02 ^ 2
  a.2.m = -0.245
  a.2.var = 0.005 ^ 2
  a.3.m = 0.0011
  a.3.var = 0.0002 ^ 2
  
  BWT ~ dunif(BWT.m - BWT.r, BWT.m + BWT.r)
  BWT.m = 7
  BWT.r = 10
  
  d18O.sw ~ dunif(d18O.sw.m - d18O.sw.r, d18O.sw.m + d18O.sw.r)
  d18O.sw.m = 1
  d18O.sw.r = 3
}
"