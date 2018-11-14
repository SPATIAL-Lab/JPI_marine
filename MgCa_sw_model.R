#####Run MgCa of seawater separately

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)
library(xlsx)

setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/")

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

##Parameters to save
parameters = c("MgCa_sw_m", "MgCa_sw_m.pre", "MgCa_sw_m.eps.ac", "MgCa_sw_m.eps")

##Run it - <2 min for 250k samples
pt = proc.time()
n.iter = 500000
n.burnin = 20000
post.mg = do.call(jags.parallel, list(model.file = "mg_model.R", parameters.to.save = parameters, 
                      data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
                      n.burnin = n.burnin, n.thin = floor((n.iter-n.burnin)/5000))) 
proc.time() - pt

save(post.mg, file = "post_mg.RData")

sl = post.mg$BUGSoutput$sims.list
sims = nrow(sl$MgCa_sw_m)
View(post.mg$BUGSoutput$summary)

plot(-10, 0, xlab="Age", ylab ="Seawater Mg/Ca", xlim=c(0,100), ylim=c(0.8,6))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(mgca_ts.ages, sl$MgCa_sw_m[i,], col = rgb(0,0,0, 0.01))
}


