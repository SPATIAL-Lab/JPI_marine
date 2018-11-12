mg_mod = " model {

#Data model for seawater MgCa observations

for(i in 1:length(MgCa_sw)){
MgCa_sw[i] ~ dnorm(MgCa_sw_m[MgCa_sw.age.ind[i]], 1 / MgCa_sw.sd[i] ^ 2)

}

#System model for MgCa_sw timeseries


for(i in 2:nmgca.ages){
  MgCa_sw_m[i] = MgCa_sw_m[i-1] * (MgCa_sw_m.eps[i] + 1)
  
  #MgCa_sw_m.eps[i] ~ dunif(MgCa_sw_m.eps[i - 1] * MgCa_sw_m.eps.ac - 0.01, MgCa_sw_m.eps[i - 1] * MgCa_sw_m.eps.ac + 0.01)
  MgCa_sw_m.eps[i] ~ dnorm(MgCa_sw_m.eps[i - 1] * MgCa_sw_m.eps.ac, MgCa_sw_m.pre)
}

#MgCa_sw_m.eps[1] ~ dunif(-0.01, 0.01)
MgCa_sw_m.eps[1] ~ dnorm(0, MgCa_sw_m.pre)
MgCa_sw_m[1] ~ dunif(MgCa_sw_m.init.min, MgCa_sw_m.init.max)
MgCa_sw_m.init.min = 1
MgCa_sw_m.init.max = 2

#Priors on MgCa_sw model parameters  

MgCa_sw_m.eps.ac ~ dunif(0.9, 1)

MgCa_sw_m.pre ~ dgamma(MgCa_sw_m.pre.shp, MgCa_sw_m.pre.rate)
MgCa_sw_m.pre.shp = 1000
MgCa_sw_m.pre.rate = 0.01

}"

parameters = c("MgCa_sw_m", "MgCa_sw_m.pre", "MgCa_sw_m.eps.ac", "MgCa_sw_m.eps")

##Set up timeseries for MgCa_sw modeling
mgca_ts.min = 110
mgca_ts.max = 0
mgca_ts.step = 1
mgca_ts.ages = seq(mgca_ts.min, mgca_ts.max, -mgca_ts.step)
mgca_ts.len = length(mgca_ts.ages)

##Read in paleo-seawater MgCa data
d_mgca_sw = read.csv("mgca_sw.csv")

##Age index for seawater MgCa samples
mgca_sw_age.ind = round((mgca_ts.min - d_mgca_sw$Age) / mgca_ts.step) + 1

##Data to pass to BUGS model
dat = list(nmgca.ages = mgca_ts.len, MgCa_sw.age.ind = mgca_sw_age.ind, MgCa_sw = d_mgca_sw$MgCa, MgCa_sw.sd = d_mgca_sw$Sigma)

mg_post = jags(model.file = textConnection(mg_mod), parameters.to.save = parameters, 
                      data = dat, inits = NULL, n.chains=3, n.iter = 50000, 
                      n.burnin = 2000, n.thin = 10) 

plot(-10, 0, xlab="Age", ylab ="Seawater Mg/Ca", xlim=c(0,100), ylim=c(0.8,6))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(mgca_ts.ages, mg_post$BUGSoutput$sims.list$MgCa_sw_m[i,], col = rgb(0,0,0, 0.01))
}

save(mg_post, file = "mg_post.RData")

