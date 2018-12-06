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

mgca_uvi = "model {

#Data model for downcore MgCa_calib constraint
LGM[1] ~ dnorm(D_MgCa_LGM.m, MgCa_calib.pre / 2)
LGM[2] ~ dnorm(D_d18O_LGM.m, d18O_calib.pre / 2)

D_MgCa_LGM.m = (a[1] + a[2] * BWT_LGM) * MgCa_sw_LGM ^ a[3] - (a[1] + a[2] * BWT_HOL) * MgCa_sw_HOL ^ a[3]
D_d18O_LGM.m = D_d18O_sw_LGM + (b[1] + b[2] * BWT_LGM + b[3] * BWT_LGM ^ 2) - (b[1] + b[2] * BWT_HOL + b[3] * BWT_HOL ^ 2)

D_d18O_sw_LGM ~ dnorm(1.1, 1 / 0.1 ^ 2) #Estimate from Adkins et al 2002
#hold seawater composition fixed to limit assocaited variance
MgCa_sw_LGM = 5.2
MgCa_sw_HOL = 5.2

BWT_LGM = BWT_HOL + D_BWT_LGM
D_BWT_LGM ~ dunif(-5, 0)
BWT_HOL ~ dnorm(3.7, 1 / 0.2 ^ 2) #Estimated from Elderfield et al 2010, fig 6

#Priors on MgCa_calib data model parameters

#Precision based on Uvigerina coretop variance
MgCa_calib.pre ~ dgamma(MgCa_calib.pre.shp, MgCa_calib.pre.rate)
MgCa_calib.pre.shp = 2
MgCa_calib.pre.rate = 1/30

#Using very loose constraints on MgCa calib parameters, holding a[3] fixed
a[1] ~ dunif(0.6, 1.4)
a[2] ~ dunif(-0.05, 0.35)
a[3] = -0.02 

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

b.1.m = 4.25
b.1.var = 0.06 ^ 2
b.2.m = -0.215
b.2.var = 0.02 ^ 2
b.3.m = -0.001
b.3.var = 0.001 ^ 2

}
"


#####
#libraries
library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)

setwd("C:/Users/gjbowen/Dropbox/Hypomirror/JPI_marine/code/")
setwd("C:/Users/u0133977/Dropbox/Hypomirror/JPI_marine/code/")

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
rmod <- jags(model.file = "mgca_calib_model.R", parameters.to.save = parameters, 
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

##Now downcore for MgCa of Uvigerina, testing viability of this method to give
##statistically robust constraints on MgCa calibration

d = read.csv("U_d18O_calib.csv")
d = d[is.na(d$Ignore),]

#Setup
parameters = c("a", "MgCa_calib.pre", "b", "d18O_calib.pre"," D_d18O_sw_LGM", 
               "D_BWT_LGM", "BWT_HOL")
rdat = list(d18O_calib.bwt.m = d$BWT, d18O_calib.bwt.sd = d$BWT_sd, d18O_calib = d$d18O_f.sw, LGM = c(-0.24, 1.7))

#Run it
set.seed(proc.time()[3])
rmod <- jags(model.file = "mgca_dc_calib_model.R", parameters.to.save = parameters, 
             data = rdat, inits = NULL, 
             n.chains=3, n.iter = 50000, n.burnin = 500, n.thin = 10)

source("helpers.R")

#This plot shows the lack of constraint on slope of MgCa calib
plotd(runif(100000, -0.05, 0.35), ylim=c(0,60), 
      xlab = expression(italic("Uvigerina") * " Mg/Ca T sensitivity (" * alpha[2] * ")")) #prior
lined(rmod$BUGSoutput$sims.list$a[,2], col="red") #posterior

#Now the coretop calibration for Uvigerina
d = read.csv("U_mgca_calib.csv")

#Assign seawater Mg/Ca estimates and uncertainty 
d$MgCa_sw = rep(5.2, nrow(d))
d$MgCa_sw_sd = rep(0.03, nrow(d))

#Setup
parameters = c("a", "mgca_pre")
rdat = list(t_m = d$BWT, t_sd = d$BWT_sd, mgca_sw_m = d$MgCa_sw, mgca_sw_sd = d$MgCa_sw_sd, mgca = d$MgCa)

#Run it
set.seed(proc.time()[3])
rmod <- jags(model.file = "mgca_calib_model.R", parameters.to.save = parameters, 
             data = rdat, inits = NULL, 
             n.chains=3, n.iter = 50000, n.burnin = 500, n.thin = 10)

lined(rmod$BUGSoutput$sims.list$a[,2], col="red", lty=2)
