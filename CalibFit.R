#####
#Independent Bayesian regressions to generate priors for proxy model parameters
#####

#d18O proxy calibration
o_mod = " model {

  for(i in 1:length(t_m)){
    d18O[i] ~ dnorm(d18O_m[i], d18O_pre)
    d18O_m[i] = b[1] + b[2] * t[i] + b[3] * t[i] ^ 2

    t[i] ~ dnorm(t_m[i], 1 / t_sd[i]^2)

  }

  d18O_pre ~ dgamma(3, 1 / 30)

  b[1] ~ dnorm(b1.m, 1 / b1.var)
  b[2] ~ dnorm(b2.m, 1 / b2.var)
  b[3] ~ dnorm(b3.m, 1 / b3.var)

  b1.m = 3.3
  b1.var = 0.5 ^ 2
  b2.m = -0.25
  b2.var = 0.1 ^ 2
  b3.m = 0 
  b3.var = 0.05 ^ 2

}
"

#Mg/Ca proxy calibration
mgca_mod = " model {

for(i in 1:length(t_m)){
mgca[i] ~ dnorm(mgca_m[i], mgca_pre)
mgca_m[i] = (a[1] + a[2] * t[i]) * mgca_sw[i] ^ a[3]

t[i] ~ dnorm(t_m[i], 1 / t_sd[i]^2)
mgca_sw[i] ~ dnorm(mgca_sw_m[i], 1 / mgca_sw_sd[i]^2) I (0.1,)
}

mgca_pre ~ dgamma(3, 1 / 30)

a[1] ~ dnorm(a1.m, 1 / a1.var)
a[2] ~ dnorm(a2.m, 1 / a2.var)
a[3] ~ dnorm(a3.m, 1 / a3.var)

a1.m = 1.4
a1.var = 0.5 ^ 2
a2.m = 0.11
a2.var = 0.1 ^ 2
a3.m = -0.019 
a3.var = 0.05 ^ 2

}
"

#####
#libraries
library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

setwd("C:/Users/gjbowen/Dropbox/Hypomirror/JPI_marine/")

##Mg/Ca

#Get the data for O. spp
d = read.csv("O_mgca_calib.csv")

#Assign seawater Mg/Ca estimates and uncertainty 
d$MgCa_sw = rep(5.2, nrow(d))
d$MgCa_sw_sd = rep(0.03, nrow(d))
for(i in 1:nrow(d)){ 
  d$MgCa_sw[i] = ifelse(d$Age[i] == 0, 5.2, 1.5)
  d$MgCa_sw_sd[i] = ifelse(d$Age[i] == 0, 0.03, 0.3)
}

#Setup
parameters = c("a", "mgca_pre")
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

#Get the data now for Uvigerina
d = read.csv("U_mgca_calib.csv")

#Assign seawater Mg/Ca estimates and uncertainty 
d$MgCa_sw = rep(5.2, nrow(d))
d$MgCa_sw_sd = rep(0.03, nrow(d))

#Setup
parameters = c("a", "mgca_pre")
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

#Get the data for Cib
d = read.csv("C_d18O_calib.csv")
d = d[is.na(d$Ignore),]

#Setup
parameters = c("b", "d18O_pre")
rdat = list(t_m = d$BWT, t_sd = d$BWT_sd, d18O = d$d18O_f.sw)

#Run it
set.seed(proc.time()[3])
rmod <- jags(model.file = textConnection(o_mod), parameters.to.save = parameters, 
             data = rdat, inits = NULL, 
             n.chains=3, n.iter = 10000, n.burnin = 500, n.thin = 5)

#Results
rmod
rmod.mcmc = as.mcmc(rmod)
plot(rmod.mcmc)

#Get the data for Uvigerina
d = read.csv("U_d18O_calib.csv")
d = d[is.na(d$Ignore),]

#Setup
parameters = c("b", "d18O_pre")
rdat = list(t_m = d$BWT, t_sd = d$BWT_sd, d18O = d$d18O_f.sw)

#Run it
set.seed(proc.time()[3])
rmod <- jags(model.file = textConnection(o_mod), parameters.to.save = parameters, 
             data = rdat, inits = NULL, 
             n.chains=3, n.iter = 10000, n.burnin = 500, n.thin = 5)

#Results
rmod
rmod.mcmc = as.mcmc(rmod)
plot(rmod.mcmc)

