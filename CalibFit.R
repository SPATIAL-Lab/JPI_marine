mgca_mod = " model {


  for(i in 1:length(t_m)){
    mgca[i] ~ dnorm(mgca_m[i], mgca_var)
    mgca_m[i] = (l[1] + l[2] * t[i]) * mgca_sw[i] ^ l[3]

    t[i] ~ dnorm(t_m[i], 1 / t_sd[i]^2)
    mgca_sw[i] ~ dnorm(mgca_sw_m[i], 1 / mgca_sw_sd[i]^2) I (0.1,)
  }

  mgca_var ~ dgamma(mgca_var.k, 1 / mgca_var.theta)
  mgca_var.k = mgca_var.m / mgca_var.theta
  mgca_var.theta = mgca_var.var / mgca_var.m
  mgca_var.m = 1 ^ 2
  mgca_var.var = 0.01

  l[1] ~ dnorm(l1.m, 1 / l1.var)
  l[2] ~ dnorm(l2.m, 1 / l2.var)
  l[3] ~ dnorm(l3.m, 1 / l3.var)

  l1.m = 1.4
  l1.var = 0.5 ^ 2
  l2.m = 0.11
  l2.var = 0.1 ^ 2
  l3.m = -0.019 
  l3.var = 0.05 ^ 2

}
"

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

parameters = c("l", "mgca_var")
set.seed(1395)

d = read.csv("mgca_calib.csv")
d$MgCa_sw = rep(5.2, nrow(d))
d$MgCa_sw_sd = rep(0.03, nrow(d))
for(i in 1:nrow(d)){ 
  d$MgCa_sw[i] = ifelse(d$Age[i] == 0, 5.2, 1.5)
  d$MgCa_sw_sd[i] = ifelse(d$Age[i] == 0, 0.03, 0.3)
}

rdat = list(t_m = d$BWT, t_sd = d$BWT_sd, mgca_sw_m = d$MgCa_sw, mgca_sw_sd = d$MgCa_sw_sd, mgca = d$MgCa)

rmod <- jags(model.file = textConnection(mgca_mod), parameters.to.save = parameters, 
                  data = rdat, inits = NULL, 
                  n.chains=3, n.iter = 10000, n.burnin = 500, n.thin = 5)

rmod
rmod.mcmc = as.mcmc(rmod)
plot(rmod.mcmc)
