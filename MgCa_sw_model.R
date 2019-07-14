#####Run MgCa of seawater separately

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)
library(xlsx)

setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/code/")

##Read in paleo-seawater MgCa data and set up ages vector
d_mgca_sw = read.csv("mgca_sw.csv")
mgca_ages = seq(0,80)
mgca_ages = c(mgca_ages, d_mgca_sw$Age)
mgca_ages = unique(mgca_ages)
mgca_ages = sort(mgca_ages, decreasing = TRUE)
mgca_ages.len = length(mgca_ages)
mgca_sw_age.ind = match(d_mgca_sw$Age, mgca_ages)

##Data to pass to BUGS model
dat = list(nmgca.ages = mgca_ages.len, mgca.ages = mgca_ages, 
           MgCa_sw.age.ind = mgca_sw_age.ind, 
           MgCa_sw = d_mgca_sw$MgCa, MgCa_sw.sd = d_mgca_sw$Sigma)

##Parameters to save
parameters = c("MgCa_sw_m", "MgCa_sw_m.pre", "MgCa_sw_m.eps.ac", "MgCa_sw_m.eps")

##Run it - <4 min for 500k samples
pt = proc.time()
n.iter = 500000
n.burnin = 20000
n.thin = floor((n.iter-n.burnin)/5000)
post.mg = do.call(jags.parallel, list(model.file = "mg_model.R", parameters.to.save = parameters, 
                      data = dat, inits = NULL, n.chains=3, n.iter = n.iter, 
                      n.burnin = n.burnin, n.thin = n.thin) )
proc.time() - pt

save(post.mg, file = "post_mg.RData")

sl = post.mg$BUGSoutput$sims.list
su = post.mg$BUGSoutput$summary
sims = nrow(sl$MgCa_sw_m)
View(su)

plot(-10, 0, xlab="Age (Ma)", ylab ="Seawater Mg/Ca", xlim=c(0,80), ylim=c(0.8,6))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(mgca_ages, sl$MgCa_sw_m[i,], col = rgb(0,0,0, 0.01))
}
lines(mgca_ages, su[1:mgca_ages.len, 5], col="red")
lines(mgca_ages, su[1:mgca_ages.len, 3], col="red", lty=3)
lines(mgca_ages, su[1:mgca_ages.len, 7], col="red", lty=3)
points(d_mgca_sw$Age, d_mgca_sw$MgCa, pch=21, bg="white")

