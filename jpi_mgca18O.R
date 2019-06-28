#####
#JPI for Mg/Ca and d18O using CRW timeseries models
#####

#####
#Preliminaries
#####

##Load libraries
library(rjags)
library(R2jags)
library(xlsx)

##My local working directories
setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/code/")
setwd("C:/Users/u0133977/Dropbox/HypoMirror/JPI_marine/code/")

##Functions for plotting and data prep
source("helpers.R")

#####
#Run inversion for different data sets, first site 806 data
#####

##Prepare site 806 data
d = prep.lear()

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "a", "MgCa_calib.pre", "b", "d18O_calib.pre", "d18O_calib.pre.2",
               "MgCa_sw_m", "MgCa_sw_m.pre", "MgCa_sw_m.eps.ac")

##Data to pass to the model
dat = list(nages = d$ts.len, ages = d$ts.ages, nmgca.ages = d$mgca_ts.len, mgca.ages = d$mgca.ages,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d$d_d18O_calib$BWT, d18O_calib.bwt.sd = d$d_d18O_calib$BWT_sd, d18O_calib = d$d_d18O_calib$d18O_f.sw,
           MgCa_sw.age.ind = d$mgca_sw_age.ind, MgCa_sw = d$d_mgca_sw$MgCa, MgCa_sw.sd = d$d_mgca_sw$Sigma,
           MgCa.age.ind = d$mgca_age.ind.all, MgCa = d$d_mgca$MgCa, 
           d18O.age.ind = d$o_age.ind, d18O = d$d_o$d18O)

##Run the inversion - ~15 hours for 1M samples

#Clock time
t1 = proc.time()

#Some parameters for the sampler
set.seed(t1[3])
n.iter = 1.5e6
n.burnin = 1e5
n.thin = floor((n.iter - n.burnin) / 5000)

#Run the MCMC
post.lear = do.call(jags.parallel, list(model.file = "split_temporal_lear.R", parameters.to.save = parameters, 
                                         data = dat, n.chains=3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

##Save the posterior samples + info
save(post.lear, file = "post_lear.RData")

#####
##Now Shackelton site U1385
#####

##Prepare the site U1385 data
d = prep.birn()

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "a", "MgCa_calib.pre", "b", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = d$ts.len, ages = d$ts.ages,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d$d_d18O_calib$BWT, d18O_calib.bwt.sd = d$d_d18O_calib$BWT_sd, d18O_calib = d$d_d18O_calib$d18O_f.sw,
           MgCa_sw.neo = d$mgca_sw_neo,
           MgCa.age.ind = d$mgca_age.ind, MgCa = d$d_mgca$MgCa, 
           d18O.age.ind = d$o_age.ind, d18O = d$d_o$d18O)

##Run the inversion - 50 min for 500k samples

#Start time
t1 = proc.time()

#Some parameters for the sampler
set.seed(t1[3])
n.iter = 50000
n.burnin = 10000
n.thin = floor(n.iter-n.burnin)/5000

#Run it
post.birn = do.call(jags.parallel, list(model.file = "split_temporal_birn.R", parameters.to.save = parameters, 
                                        data = dat, n.chains=3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

##Save posterior samples + info
save(post.birn, file = "post_birn.RData")

#####
##Now site 1123
#####

##Prep the site 1123 data
d = prep.elder()

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "a", "MgCa_calib.pre", "b", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = d$ts.len, ages = d$ts.ages,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d$d_d18O_calib$BWT, d18O_calib.bwt.sd = d$d_d18O_calib$BWT_sd, d18O_calib = d$d_d18O_calib$d18O_f.sw,
           MgCa_sw.neo = d$mgca_sw_neo,
           MgCa.age.ind = d$mgca_age.ind, MgCa = d$d_mgca$MgCa, 
           d18O.age.ind = d$o_age.ind, d18O = d$d_o$d18O)

##Run the inversion - 50 min for 500k samples

#Start time
t1 = proc.time()

#Some parameters for the sampler
set.seed(t1[3])
n.iter = 500000
n.burnin = 10000
n.thin = floor(n.iter-n.burnin)/5000

#Run it
post.elder = do.call(jags.parallel, list(model.file = "split_temporal_elder.R", parameters.to.save = parameters, 
                                         data = dat, n.chains=3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

##Save the posterior samples + info
save(post.elder, file = "post_elder.RData")

#####
##Now sites U1385 and 1123 together
#####

##Prep all the data
d = prep.multi()

##Parameters to be saved
parameters = c("d18O_sw.b", "BWT.b", "d18O_sw.e", "BWT.e", 
               "BWT.b.eps.ac", "BWT.b.pre", "d18O_sw.b.eps.ac", "d18O_sw.b.pre", 
               "BWT.e.eps.ac", "BWT.e.pre", "d18O_sw.e.eps.ac", "d18O_sw.e.pre", 
               "a", "MgCa_calib.pre", "b.c", "b.u", "d18O_calib.c.pre", "d18O_calib.u.pre")

##Data to pass to the model
dat = list(nages = d$ts.len, d$ts.ages, 
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.u.bwt.m = d$d_d18O_calib.u$BWT, d18O_calib.u.bwt.sd = d$d_d18O_calib.u$BWT_sd, d18O_calib.u = d$d_d18O_calib.u$d18O_f.sw,
           d18O_calib.c.bwt.m = d$d_d18O_calib.c$BWT, d18O_calib.c.bwt.sd = d$d_d18O_calib.c$BWT_sd, d18O_calib.c = d$d_d18O_calib.c$d18O_f.sw,
           MgCa_sw.neo = d$mgca_sw_neo,
           MgCa.age.ind.b = d$mgca_age.ind.b, MgCa.b = d$d_mgca.b$MgCa, 
           d18O.age.ind.b = d$o_age.ind.b, d18O.b = d$d_o.b$d18O,
           MgCa.age.ind.e = d$mgca_age.ind.e, MgCa.e = d$d_mgca.e$MgCa,
           d18O.age.ind.e = d$o_age.ind.e, d18O.e = d$d_o.e$d18O)

##Run the inversion ~ 1 hour for 100k sims

#Start time
t1 = proc.time()

#Some parameters for the sampler
set.seed(t1[3])
n.iter = 750000
n.burnin = 10000
n.thin = floor((n.iter - n.burnin) / 5000)

#Run it
post.multi = do.call(jags.parallel, list(model.file = "split_temporal_multi.R", parameters.to.save = parameters, 
                                         data = dat, n.chains=3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1

##Save posterior samples + info
save(post.multi, file="post_multi.RData")

#####
#Bayesian inversion of Mg/Ca proxy calibration data comparing the core-top and
#down-core approaches of Elderfield et al. (2010)
#####

#First the down-core method
#Read d18O calibration data needed for down-core constraint
d = read.csv("U_d18O_calib.csv")
d = d[is.na(d$Ignore),]

#Now the down-core data themselves
#These are the 1123 and CHAT 1K coretop data plus 1123 data from 25 to 15 ka
#from SI of Elderfield et al (2012)
LGM = read.csv("1123_LGM.csv")
HOL = read.csv("1123_HOL.csv")

#Organize the down-core records
MgCa_LGM = LGM[!is.na(LGM$MgCa),]
d18O_LGM = LGM[!is.na(LGM$d18O),]
MgCa_HOL = HOL[!is.na(HOL$MgCa),]
d18O_HOL = HOL[!is.na(HOL$d18O),]

#Parameters to save
parameters = c("a", "MgCa_calib.pre", "b", "d18O_calib.pre"," D_d18O_sw_LGM", 
               "D_BWT_LGM", "BWT_HOL")

#Data to pass to the model
rdat = list(d18O_calib.bwt.m = d$BWT, d18O_calib.bwt.sd = d$BWT_sd, 
            d18O_calib = d$d18O_f.sw, MgCa_LGM = MgCa_LGM$MgCa, d18O_LGM = d18O_LGM$d18O,
            MgCa_HOL = MgCa_HOL$MgCa, d18O_HOL = d18O_HOL$d18O)

#Run it
set.seed(proc.time()[3])
post.mgca.dc <- jags(model.file = "mgca_dc_calib_model.R", parameters.to.save = parameters, 
                     data = rdat, inits = NULL, 
                     n.chains=3, n.iter = 50000, n.burnin = 500, n.thin = 10)

#Save the results
save(post.mgca.dc, file = "post_mgca_dc.RData")

#Now the core-top calibration
#Read the data
d = read.csv("U_mgca_calib.csv")

#Assign seawater Mg/Ca estimates and uncertainty, these are all modern data 
d$MgCa_sw = rep(5.2, nrow(d))
d$MgCa_sw_sd = rep(0.03, nrow(d))

#Parameters to save
parameters = c("a", "mgca_pre")

#Data to pass to the model
rdat = list(t_m = d$BWT, t_sd = d$BWT_sd, mgca_sw_m = d$MgCa_sw, mgca_sw_sd = d$MgCa_sw_sd, mgca = d$MgCa)

#Run it
set.seed(proc.time()[3])
post.mgca.ct <- jags(model.file = "mgca_calib_model.R", parameters.to.save = parameters, 
                     data = rdat, inits = NULL, 
                     n.chains=3, n.iter = 50000, n.burnin = 500, n.thin = 10)

#Save the output
save(post.mgca.ct, file = "post_mgca_ct.RData")
