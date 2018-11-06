#####
#Independent Bayesian regressions to generate priors for proxy model parameters
#####

#d18O proxy calibration
o_mod = " model {

  for(i in 1:length(t_m)){
    d18O[i] ~ dnorm(d18O_m[i], d18O_pre)
    d18O_m[i] = a[1] + a[2] * t[i] + a[3] * t[i] ^ 2

    t[i] ~ dnorm(t_m[i], 1 / t_sd[i]^2)

  }

  d18O_pre ~ dgamma(3, 1 / 30)

  a[1] ~ dnorm(a1.m, 1 / a1.var)
  a[2] ~ dnorm(a2.m, 1 / a2.var)
  a[3] ~ dnorm(a3.m, 1 / a3.var)

  a1.m = 3.3
  a1.var = 0.5 ^ 2
  a2.m = -0.25
  a2.var = 0.1 ^ 2
  a3.m = 0 
  a3.var = 0.05 ^ 2

}
"

#Mg/Ca proxy calibration
mgca_mod = " model {

for(i in 1:length(t_m)){
mgca[i] ~ dnorm(mgca_m[i], mgca_pre)
mgca_m[i] = (l[1] + l[2] * t[i]) * mgca_sw[i] ^ l[3]

t[i] ~ dnorm(t_m[i], 1 / t_sd[i]^2)
mgca_sw[i] ~ dnorm(mgca_sw_m[i], 1 / mgca_sw_sd[i]^2) I (0.1,)
}

mgca_pre ~ dgamma(3, 1 / 30)

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

#####
#libraries
library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

##Mg/Ca

#Get the data
d = read.csv("mgca_calib.csv")

#Assign seawater Mg/Ca estimates and uncertainty 
d$MgCa_sw = rep(5.2, nrow(d))
d$MgCa_sw_sd = rep(0.03, nrow(d))
for(i in 1:nrow(d)){ 
  d$MgCa_sw[i] = ifelse(d$Age[i] == 0, 5.2, 1.5)
  d$MgCa_sw_sd[i] = ifelse(d$Age[i] == 0, 0.03, 0.3)
}

#Setup
parameters = c("l", "mgca_pre")
rdat = list(t_m = d$BWT, t_sd = d$BWT_sd, mgca_sw_m = d$MgCa_sw, mgca_sw_sd = d$MgCa_sw_sd, mgca = d$MgCa)

#Run it
set.seed(proc.time()[3])
rmod <- jags(model.file = textConnection(mgca_mod), parameters.to.save = parameters, 
                  data = rdat, inits = NULL, 
                  n.chains=3, n.iter = 10000, n.burnin = 500, n.thin = 5)

#Results
rmod
rmod.mcmc = as.mcmc(rmod)
plot(rmod.mcmc)

##d18O

#Get the data
d = read.csv("C_comp.csv")

#Setup
parameters = c("a", "d18O_pre")
rdat = list(t_m = d$Temperature_C, t_sd = rep(0.2, nrow(d)), d18O = d$C.SW_d18O)

#Run it
set.seed(proc.time()[3])
rmod <- jags(model.file = textConnection(o_mod), parameters.to.save = parameters, 
             data = rdat, inits = NULL, 
             n.chains=3, n.iter = 10000, n.burnin = 500, n.thin = 5)

#Results
rmod
rmod.mcmc = as.mcmc(rmod)
plot(rmod.mcmc)