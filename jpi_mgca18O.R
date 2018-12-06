#####
#JPI for Mg/Ca and d18O using CRW timeseries models
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
#Running inversion for different data sets, first Lear 2015 data
#####

d = prep.lear()

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "a", "MgCa_calib.pre", "b", "d18O_calib.pre", "d18O_calib.pre.2",
               "MgCa_sw_m", "MgCa_sw_m.pre", "MgCa_sw_m.eps.ac")

##Data to pass to the model
dat = list(nages = d$ts.len, nmgca.ages = d$mgca_ts.len,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d$d_d18O_calib$BWT, d18O_calib.bwt.sd = d$d_d18O_calib$BWT_sd, d18O_calib = d$d_d18O_calib$d18O_f.sw,
           MgCa_sw.age.ind = d$mgca_sw_age.ind, MgCa_sw = d$d_mgca_sw$MgCa, MgCa_sw.sd = d$d_mgca_sw$Sigma,
           MgCa.age.ind = d$mgca_age.ind.all, MgCa = d$d_mgca$MgCa, 
           d18O.age.ind = d$o_age.ind, d18O = d$d_o$d18O)

##Run the inversion - 15 hours for 1M samples
t1 = proc.time()
set.seed(t1[3])
n.iter = 1.5e6
n.burnin = 1e5
n.thin = floor((n.iter - n.burnin) / 5000)
post.lear = do.call(jags.parallel, list(model.file = "split_temporal_lear.R", parameters.to.save = parameters, 
                                         data = dat, n.chains=3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))
proc.time() - t1

save(post.lear, file = "post_lear.RData")

#####
##Now Shackelton site U1385
#####

d = prep.birn()

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "a", "MgCa_calib.pre", "b", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = d$ts.len,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d$d_d18O_calib$BWT, d18O_calib.bwt.sd = d$d_d18O_calib$BWT_sd, d18O_calib = d$d_d18O_calib$d18O_f.sw,
           MgCa_sw.neo = d$mgca_sw_neo,
           MgCa.age.ind = d$mgca_age.ind, MgCa = d$d_mgca$MgCa, 
           d18O.age.ind = d$o_age.ind, d18O = d$d_o$d18O)

##Run the inversion - 50 min for 500k samples
t1 = proc.time()
set.seed(t1[3])
n.iter = 500000
n.burnin = 10000
n.thin = floor(n.iter-n.burnin)/5000
post.birn = do.call(jags.parallel, list(model.file = "split_temporal_birn.R", parameters.to.save = parameters, 
                                        data = dat, n.chains=3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))
proc.time() - t1

save(post.birn, file = "post_birn.RData")

#####
##Now the Elderfield record site 1123
#####

d = prep.elder()

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "a", "MgCa_calib.pre", "b", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = d$ts.len,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d$d_d18O_calib$BWT, d18O_calib.bwt.sd = d$d_d18O_calib$BWT_sd, d18O_calib = d$d_d18O_calib$d18O_f.sw,
           MgCa_sw.neo = d$mgca_sw_neo,
           MgCa.age.ind = d$mgca_age.ind, MgCa = d$d_mgca$MgCa, 
           d18O.age.ind = d$o_age.ind, d18O = d$d_o$d18O)

##Run the inversion - 50 min for 500k samples
t1 = proc.time()
set.seed(t1[3])
n.iter = 500000
n.burnin = 10000
n.thin = floor(n.iter-n.burnin)/5000
post.elder = do.call(jags.parallel, list(model.file = "split_temporal_elder.R", parameters.to.save = parameters, 
                                         data = dat, n.chains=3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))
proc.time() - t1

save(post.elder, file = "post_elder.RData")

#####
##Now Shackelton site and 1123 together
#####

d = prep.multi()

##Parameters to be saved
parameters = c("d18O_sw.b", "BWT.b", "d18O_sw.e", "BWT.e", 
               "BWT.b.eps.ac", "BWT.b.pre", "d18O_sw.b.eps.ac", "d18O_sw.b.pre", 
               "BWT.e.eps.ac", "BWT.e.pre", "d18O_sw.e.eps.ac", "d18O_sw.e.pre", 
               "a", "MgCa_calib.pre", "b.c", "b.u", "d18O_calib.c.pre", "d18O_calib.u.pre")

##Data to pass to the model
dat = list(nages = d$ts.len,
           MgCa_calib.bwt.m = d$d_mgca_calib$BWT, MgCa_calib.bwt.sd = d$d_mgca_calib$BWT_sd, MgCa_calib = d$d_mgca_calib$MgCa,
           d18O_calib.u.bwt.m = d$d_d18O_calib.u$BWT, d18O_calib.u.bwt.sd = d$d_d18O_calib.u$BWT_sd, d18O_calib.u = d$d_d18O_calib.u$d18O_f.sw,
           d18O_calib.c.bwt.m = d$d_d18O_calib.c$BWT, d18O_calib.c.bwt.sd = d$d_d18O_calib.c$BWT_sd, d18O_calib.c = d$d_d18O_calib.c$d18O_f.sw,
           MgCa_sw.neo = d$mgca_sw_neo,
           MgCa.age.ind.b = d$mgca_age.ind.b, MgCa.b = d$d_mgca.b$MgCa, 
           d18O.age.ind.b = d$o_age.ind.b, d18O.b = d$d_o.b$d18O,
           MgCa.age.ind.e = d$mgca_age.ind.e, MgCa.e = d$d_mgca.e$MgCa,
           d18O.age.ind.e = d$o_age.ind.e, d18O.e = d$d_o.e$d18O)

##Run the inversion ~ 1 hour for 100k sims
t1 = proc.time()
set.seed(t1[3])
n.iter = 750000
n.burnin = 10000
n.thin = floor((n.iter - n.burnin) / 5000)
post.multi = do.call(jags.parallel, list(model.file = "split_temporal_multi.R", parameters.to.save = parameters, 
                                         data = dat, n.chains=3, n.iter = n.iter, 
                                         n.burnin = n.burnin, n.thin = n.thin))
proc.time() - t1

save(post.multi, file="post_multi.RData")


#####
##Plots for the manuscript
#####

##Some preliminaries and shorthand for plotting ease
load("post_lear.RData")
load("post_multi.RData")
load("post_birn.RData")
load("post_elder.RData")

##These first plots use the Lear data, so prep it using the function above
d = prep.lear()

#Shorthand
sl = post.lear$BUGSoutput$sims.list
su = post.lear$BUGSoutput$summary

#Get some indicies
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))

##Figure 3: the modeled timeseries
png("../Figure03.png", units="in", width=5, height=5, res=300)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(0,18), ylim=c(-3,9), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 5], col="red")
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 3], col="red", lty=3)
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 7], col="red", lty=3)
tp = d$d_mgca[order(d$d_mgca$Age.Ma), "Age.Ma"]
points(tp, rep(-3, nrow(d$d_mgca)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(0,18), ylim=c(2,-1.5))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(d$ts.ages, su[d18O.start:(d18O.start+d$ts.len-1), 5], col="red")
lines(d$ts.ages, su[d18O.start:(d18O.start+d$ts.len-1), 3], col="red", lty=3)
lines(d$ts.ages, su[d18O.start:(d18O.start+d$ts.len-1), 7], col="red", lty=3)
op = d$d_o[order(d$d_o$Age.Ma),"Age.Ma"]
points(op, rep(2, nrow(d$d_o)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

##Figure 2: The Mg/Ca sw time series
png("../Figure02.png", units="in", width=5, height=2.75, res=300)
par(mar=c(4,4,1,1), cex=0.85)
plot(-10, 0, xlab="Age (Ma)", ylab ="Seawater Mg/Ca", xlim=c(0,80), ylim=c(0.8,5.5))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$mgca_ts.ages, sl$MgCa_sw_m[i,], col = rgb(0,0,0, 0.01))
}
lines(d$mgca_ts.ages, su[MgCa.start:(MgCa.start+d$mgca_ts.len-1), 5], col="red")
lines(d$mgca_ts.ages, su[MgCa.start:(MgCa.start+d$mgca_ts.len-1), 3], col="red", lty=3)
lines(d$mgca_ts.ages, su[MgCa.start:(MgCa.start+d$mgca_ts.len-1), 7], col="red", lty=3)

#curve fit from Lear 15 for comparison
ages = seq(0, 55, 1)
vals = 5.2 - 0.238 * ages + 0.00661 * ages^2 - 6.66e-5 * ages^3
lines(ages, vals, lty=2)

#points showing Mg/Ca proxy obs and distribution of proxy and calib data
points(d$d_mgca_sw$Age, d$d_mgca_sw$MgCa, pch=21, bg = "white")
points(d$d_mgca$Age.Ma, rep(0.8, nrow(d$d_mgca)), pch=21, bg = "black")
calib_ages = d$d_mgca_calib$Age
calib_ages = calib_ages[calib_ages>0]
points(calib_ages, rep(0.8, length(calib_ages)), pch=21, bg = "grey")

dev.off()

##Figure 5: Prior/posterior plots for calibration parameters
png("../Figure05.png", res = 300, units = "in", width = 8, height = 4)
layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 2, ncol = 4, byrow=TRUE))
par(mai=c(0.5,0.5,0.1,0.1))
xoff = 2.3

#MgCa calibration parms
plotd(sl$a[,1], col="red")
lined(rnorm(100000, 1.5, 0.1))
title(xlab=expression(paste(alpha[1])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

plotd(sl$a[,2], col="red", ylab="")
lined(rnorm(100000, 0.1, 0.01))
title(xlab=expression(paste(alpha[2])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

plotd(sl$a[,3], col="red", ylab="")
lined(rnorm(100000, -0.02, 0.03))
title(xlab=expression(paste(alpha[3])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(c)")

plotd(sqrt(1/(sl$MgCa_calib.pre)), col="red", xlim=c(0.05,0.25), ylab="")
lined(sqrt(1/(rgamma(100000, 2, 1/30))))
title(xlab=expression(paste(sigma["MgCaf"])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(d)")

#d18O calibration parms
plotd(sl$b[,1], col="red")
lined(rnorm(100000, 3.32, 0.02))
title(xlab=expression(paste(beta[1])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(e)")

plotd(sl$b[,2], col="red", ylab="")
lined(rnorm(100000, -0.237, 0.01))
title(xlab=expression(paste(beta[2])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(f)")

plotd(sl$b[,3], col="red", ylab="")
lined(rnorm(100000, 0.001, 0.0005))
title(xlab=expression(paste(beta[3])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(g)")

plotd(sqrt(1/(sl$d18O_calib.pre)), col="red", xlim=c(0.05,0.7), ylab="")
lined(sqrt(1/(sl$d18O_calib.pre.2)), col="red", lty=2)
lined(sqrt(1/(rgamma(100000, 3, 1/30))))
lined(sqrt(1/(rgamma(100000, 6, 1))), lty=2)
title(xlab=expression(paste(sigma[paste(delta, "18Of")])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(h)")

dev.off()

##Figure 7: Parameter covariance
png("../Figure07.png", res=300, units="in", width=6, height=4)
layout(matrix(c(1,2,3,4,5,6), nrow=2, byrow=TRUE))
par(mar=c(4,4,0.4,0.4))

smoothScatter(sl$a[,1], sl$a[,2], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")
title(xlab=expression(paste(alpha[1])), line = 2.5)
title(ylab=expression(paste(alpha[2])), line = 2.5)

smoothScatter(sl$a[,1], sl$a[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")
title(xlab=expression(paste(alpha[1])), line = 2.5)
title(ylab=expression(paste(alpha[3])), line = 2.5)

smoothScatter(sl$a[,2], sl$a[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(c)")
title(xlab=expression(paste(alpha[2])), line = 2.5)
title(ylab=expression(paste(alpha[3])), line = 2.5)

smoothScatter(sl$b[,1], sl$b[,2], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(d)")
title(xlab=expression(paste(beta[1])), line = 2.5)
title(ylab=expression(paste(beta[2])), line = 2.5)

smoothScatter(sl$b[,1], sl$b[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(e)")
title(xlab=expression(paste(beta[1])), line = 2.5)
title(ylab=expression(paste(beta[3])), line = 2.5)

smoothScatter(sl$b[,2], sl$b[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(f)")
title(xlab=expression(paste(beta[2])), line = 2.5)
title(ylab=expression(paste(beta[3])), line = 2.5)

dev.off()

##Correlation matricies for the plots above
a.cor = cor(sl$a)
b.cor = cor(sl$b)

##Figure 9 showing 2-d posterior density of environmental time series parameters 

#First calculate the difference of each BWT and d18O_sw values relative to
#the 18Ma value for that posterior draw
D_BWT = matrix(double(), nrow = 15000, ncol = 361)
D_d18O_sw = matrix(double(), nrow = 15000, ncol = 361)
for(i in 1:15000){
  for(j in  1:361){
    D_BWT[i,j] = sl$BWT[i,j] - sl$BWT[i,1]
    D_d18O_sw[i,j] = sl$d18O_sw[i,j] - sl$d18O_sw[i,1]
  }
}

#Now calculate the medians of the values obtained above
D_BWT.m = double()
D_d18O_sw.m = double()
for(i in 1:361){
  D_BWT.m[i] = median(D_BWT[,i])
  D_d18O_sw.m[i] = median(D_d18O_sw[,i])
}

#2-d kernal density for the differences
library(MASS)
DKE = kde2d(D_BWT, D_d18O_sw, h=c(1, 0.5), n=50)

#Plot it
png("../Figure09.png", res=300, units="in", width=3, height=3)
par(mar=c(5,5,0.5,0.5), cex=0.75)
smoothScatter(D_BWT, D_d18O_sw, xlab=expression(Delta*"BWT ("*degree*" C)"),
              ylab = expression(Delta*delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
              xlim = c(-9,4.5), ylim = c(2.75,-1.5), col="white")
contour(DKE, add=TRUE, drawlabels=FALSE, col="grey")
pal = heat.colors(361)
points(D_BWT.m, D_d18O_sw.m, pch=19, col=pal, cex=0.15)

dev.off()

##Figure 8: 9 panel posterior density plots for time series model parameters

#Set it up
png("../Figure08.png", res=300, units="in", width = 5.6, height = 5.2)
layout(matrix(seq(1,9), nrow=3, byrow = TRUE), 
       widths=c(lcm(2*2.54), lcm(1.8*2.54), lcm(1.8*2.54)),
       heights=c(lcm(1.6*2.54), lcm(1.6*2.54), lcm(2*2.54)))

#Load Lear results
sl = post.lear$BUGSoutput$sims.list

#BWT and d18O timeseries autocorrelation
par(mai = c(0.1,0.5,0.1,0.1))
plotd(sl$BWT.eps.ac, col="red", ylab="", ylim=c(0,20), xlim=c(0,0.8), axes=FALSE)
axis(1, labels = FALSE)
axis(2)
box()
lined(sl$d18O_sw.eps.ac, col="red", lty=2)
lines(c(0,0.4), c(2.5,2.5))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

#BWT TS SD
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/sl$BWT.pre), col="red", ylab="", xlim=c(0.2,0.5), axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/rgamma(100000, 20, 2)))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

#d18O TS SD
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/(sl$d18O_sw.pre)), col="red", ylab="", xlim = c(0.075, 0.21), lty=2, axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/(rgamma(100000, 10, 1/5))), lty=2)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(c)")

#load multi
sl = post.multi$BUGSoutput$sims.list

#BWT and d18O timeseries autocorrelation
par(mai = c(0.1,0.5,0.1,0.1))
plotd(sl$BWT.b.eps.ac, col="red", ylim=c(0,6), xlim=c(0,0.8), axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sl$d18O_sw.b.eps.ac, col="red", lty=2)
lines(c(0,0.8), c(2.5,2.5))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(d)")

#BWT TS SD
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/sl$BWT.b.pre), col="red", ylab="", xlim=c(0.2,0.5), axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/rgamma(100000, 20, 2)))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(e)")

#d18O TS SD
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/(sl$d18O_sw.b.pre)), col="red", ylab="", xlim = c(0.075, 0.21), lty=2, axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/(rgamma(100000, 10, 1/5))), lty=2)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(f)")

#BWT and d18O timeseries autocorrelation
par(mai = c(0.5,0.5,0.1,0.1))
plotd(sl$BWT.e.eps.ac, col="red", ylab="", ylim=c(0,6), xlim=c(0,0.8))
lined(sl$d18O_sw.e.eps.ac, col="red", lty=2)
lines(c(0,0.8), c(2.5,2.5))
title(xlab=expression(phi), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(g)")

#BWT TS SD
par(mai = c(0.5,0.3,0.1,0.1))
plotd(sqrt(1/sl$BWT.e.pre), col="red", ylab="", xlim=c(0.2,0.5))
lined(sqrt(1/rgamma(100000, 20, 2)))
title(xlab=expression(paste(sigma ["BWT"])), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(h)")

#d18O TS SD
par(mai = c(0.5,0.3,0.1,0.1))
plotd(sqrt(1/(sl$d18O_sw.e.pre)), col="red", ylab="", xlim = c(0.075, 0.21), lty=2)
lined(sqrt(1/(rgamma(100000, 10, 1/5))), lty=2)
title(xlab=expression(sigma [paste(delta, "18Osw")]), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(i)")

dev.off()

##Figure 4: Site 1123 and U1385 posterior time series

#These use data from Multi run
d = prep.multi()
db = read.csv("birner_2016_interp.csv")
de = read.csv("elderfield_2012_interp.csv")

#Shorthand
sl = post.multi$BUGSoutput$sims.list
su = post.multi$BUGSoutput$summary

#Get some indicies
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))

#Setup
png("../Figure04.png", units="in", width=5, height=5, res=300)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)

#First panel w BWT
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1239,1315), ylim=c(-3,8), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$BWT.b[i,], col = rgb(0.5,0,0, 0.01))
}
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$BWT.e[i,], col = rgb(0,0,1, 0.01))
}
lines(d$ts.ages, su[BWT.b.start:(BWT.b.start + d$ts.len - 1), 5], col="red")
lines(d$ts.ages, su[BWT.b.start:(BWT.b.start + d$ts.len - 1), 3], col="red", lty=3)
lines(d$ts.ages, su[BWT.b.start:(BWT.b.start + d$ts.len - 1), 7], col="red", lty=3)
lines(d$ts.ages, su[BWT.e.start:(BWT.e.start + d$ts.len - 1), 5], col=rgb(0,0,0.7))
lines(d$ts.ages, su[BWT.e.start:(BWT.e.start + d$ts.len - 1), 3], col=rgb(0,0,0.7), lty=3)
lines(d$ts.ages, su[BWT.e.start:(BWT.e.start + d$ts.len - 1), 7], col=rgb(0,0,0.7), lty=3)

#add reconstructions from original papers 
arrows(1279.45, 3.5-1.74, 1279.45, 3.5+1.74, code=3, angle=90, length=0.06, col="red")
points(db$Age_ka, db$BWT, pch=21, cex=0.5, bg="red")
arrows(1273.7, -0.543-2, 1273.7, -0.543+2, code=3, angle=90, length=0.06, col=rgb(0,0,0.7))
points(de$Age.ka, de$BWT, pch=21, cex=0.5, bg=rgb(0,0,0.7))

#panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (ka)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1239,1315), ylim=c(1.5,-0.6))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$d18O_sw.b[i,], col = rgb(0.5,0,0, 0.01))
}
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$d18O_sw.e[i,], col = rgb(0,0,1, 0.01))
}
lines(d$ts.ages, su[d18O.b.start:(d18O.b.start+d$ts.len-1), 5], col="red")
lines(d$ts.ages, su[d18O.b.start:(d18O.b.start+d$ts.len-1), 3], col="red", lty=3)
lines(d$ts.ages, su[d18O.b.start:(d18O.b.start+d$ts.len-1), 7], col="red", lty=3)
lines(d$ts.ages, su[d18O.e.start:(d18O.e.start+d$ts.len-1), 5], col=rgb(0,0,0.7))
lines(d$ts.ages, su[d18O.e.start:(d18O.e.start+d$ts.len-1), 3], col=rgb(0,0,0.7), lty=3)
lines(d$ts.ages, su[d18O.e.start:(d18O.e.start+d$ts.len-1), 7], col=rgb(0,0,0.7), lty=3)

#add reconstructions from original papers 
arrows(1279.45, 0.45-0.46, 1279.45, 0.45+0.46, code=3, angle=90, length=0.06, col="red")
lines(db$Age_ka, db$d18O_sw, lty=2, col="red")
arrows(1273.7, -0.079-0.4, 1273.7, -0.079+0.4, code=3, angle=90, length=0.06, col=rgb(0,0,0.7))
lines(de$Age.ka, de$d18O_sw, lty=2, col=rgb(0,0,0.7))

#panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

##Figure 10, 2 panels, first shows change in BWT values as a function of time

#using Lear analysis
d = prep.lear()

#Shorthand
sl = post.lear$BUGSoutput$sims.list
su = post.lear$BUGSoutput$summary

#Get some indicies
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))

#This one calculates change relative to modern
BWT.delta = matrix(rep(0, sims * (d$ts.len)), nrow = sims, ncol = (d$ts.len))
for(i in 1:sims){
  for(j in 1:(d$ts.len)){ BWT.delta[i,j] = sl$BWT[i,j] - sl$BWT[i,d$ts.len]} 
}

#This gets the zero change value from the emperical CDF of the change time series 
BWT.delta.p = double()
for(j in 1:d$ts.len){
  tst = ecdf(BWT.delta[,j])
  BWT.delta.p[j] = tst(0)
}

#This gets the CDF value of the modern median from each time steps ECDF 
trad = su[(BWT.start + d$ts.len - 1), 5]
BWT.p = double()
for(j in 1:d$ts.len){
  tst = ecdf(sl$BWT[,j])
  BWT.p[j] = tst(trad)
}

#Now a plot showing probabilities on delta T relative to modern
png("../Figure10.png", res=300, units="in", width = 5, height = 5.5)
layout(matrix(c(1,2), nrow = 2))
par(mar = c(4,4.5,1,4), cex = 0.85)
plot(-10, 0, xlab = "", ylab = "", xlim=c(0.05,2), ylim=c(-2.5,2))
title(xlab = "Age (Ma)", line = 2.75)
title(ylab = expression("BWT ("*degree*" C)"), line = 2.75)
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 5], col="red")
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 3], col="red", lty=3)
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 7], col="red", lty=3)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

par(new=TRUE)
plot(d$ts.ages[1:(d$ts.len-1)], BWT.delta.p[1:(d$ts.len-1)], xlim=c(0.05,2), ylim=c(5e-3, 1), type="l", 
     axes=FALSE, log="y", xlab="", ylab="", col="blue")
lines(d$ts.ages[1:(d$ts.len-1)], BWT.p[1:(d$ts.len-1)], lty=3, col="blue")
lines(c(-1,6), c(0.05,0.05), lty=2, col="blue")
lines(c(-1,6), c(0.01,0.01), lty=2, col="blue")
axis(4, at=c(5e-5,5e-4,5e-3,5e-2,5e-1))
mtext(side = 4, "Zero change probability", line = 2.75)

##Now switch over to the multi analysis
d = prep.multi()

#Shorthand
sl = post.multi$BUGSoutput$sims.list
su = post.multi$BUGSoutput$summary

#Get some indicies
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))

##Calculate paired difference between d18O_sw samples for two records
#Get quantile value for zero difference
d18O_sw.delta = sl$d18O_sw.b - sl$d18O_sw.e
d18O_sw.ptiles = matrix(double(), ncol = ncol(d18O_sw.delta), nrow = 4)
for(j in 1:ncol(d18O_sw.delta)){
  d18O_sw.ptiles[1:3,j] = quantile(d18O_sw.delta[,j], c(0.025,0.5,0.975))
  tst = ecdf(d18O_sw.delta[,j])
  d18O_sw.ptiles[4,j] = tst(0)
}

##Plot difference in d18Osw and zero difference probs
plot(-10, 0, xlab = "Age (ka)", ylab = expression(Delta*delta^{18}*"O"[sw]*" (U1385 - 1123)"), 
     xlim=c(1239,1315), ylim=c(-0.5,2))
for(i in seq(1, sims, by = max(floor(sims / 1500),1))){
  lines(d$ts.ages, d18O_sw.delta[i,], col = rgb(0,0,0, 0.01))
}
lines(d$ts.ages, d18O_sw.ptiles[2,], col="red")
lines(d$ts.ages, d18O_sw.ptiles[1,], col="red", lty=3)
lines(d$ts.ages, d18O_sw.ptiles[3,], col="red", lty=3)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

par(new = TRUE)
plot(d$ts.ages, pmax(pmin(d18O_sw.ptiles[4,],1-d18O_sw.ptiles[4,]), 5e-4), type="l", log="y", 
     axes = FALSE, xlim=c(1239,1315), ylim=c(5e-4,5e-1), xlab="", ylab="", col="blue")
axis(4, at=c(5e-4,5e-3,5e-2,5e-1))
lines(c(1230,1320), c(0.05,0.05), lty=2, col="blue")
lines(c(1230,1320), c(0.01,0.01), lty=2, col="blue")

##Calculate paired difference between d18O_sw samples for two records
#Get quantile value for zero difference
d18O_sw.delta = post.birn$BUGSoutput$sims.list$d18O_sw - post.elder$BUGSoutput$sims.list$d18O_sw
d18O_sw.ptiles = matrix(double(), ncol = ncol(d18O_sw.delta), nrow = 4)
for(j in 1:ncol(d18O_sw.delta)){
  d18O_sw.ptiles[1:3,j] = quantile(d18O_sw.delta[,j], c(0.025,0.5,0.975))
  tst = ecdf(d18O_sw.delta[,j])
  d18O_sw.ptiles[4,j] = tst(0)
}
lines(d$ts.ages, pmax(pmin(d18O_sw.ptiles[4,],1-d18O_sw.ptiles[4,]), 5e-4), lty=3, col="blue")
mtext(side = 4, "Zero difference probability", line = 2.75)

dev.off()

##Figure 6 compares constraints of downcore 1123 data and coretop data 
##on Uvigerina Mg/Ca T sensitivity; downcore modeled on Elderfield et al, 2010

#Read d18O calibration data needed for downcore constraint
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

#needed for plotting
source("helpers.R")

#First plot prior and results of downcore analysis
png("../Figure06.png", res=300, units="in", width=3, height=3)
par(mar=c(5,5,0.5,0.5), cex=0.75)
plotd(runif(100000, -0.05, 0.35), ylim=c(0,100), 
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

dev.off()




#####
##The rest is exploratory code, few useful stats, etc.
#####


##Get some stats
#95%CI width for environmental records
#Lear
prep.lear()

#Shorthand
sl = post.lear$BUGSoutput$sims.list
su = post.lear$BUGSoutput$summary

#Get some indicies
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))

su.diff = as.double(su[,7] - su[,3])
mean(su.diff[BWT.start:(BWT.start+ts.len-1)])
mean(su.diff[d18O.start:(d18O.start+ts.len-1)])

#Multi
prep.multi()

#Shorthand
sl = post.multi$BUGSoutput$sims.list
su = post.multi$BUGSoutput$summary

#Get some indicies
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))

su.diff = as.double(su[,7] - su[,3])
mean(su.diff[BWT.b.start:(BWT.b.start+ts.len-1)])
mean(su.diff[d18O.b.start:(d18O.b.start+ts.len-1)])
mean(su.diff[BWT.e.start:(BWT.e.start+ts.len-1)])
mean(su.diff[d18O.e.start:(d18O.e.start+ts.len-1)])

#Shorthand
sl = post.birn$BUGSoutput$sims.list
su = post.birn$BUGSoutput$summary

#Show summary
View(su)

#Get some indicies
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))

##A couple of standard plots of the modeled timeseries
png("T_18O_bir.png", units="in", width=5, height=5, res=300)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1239,1320), ylim=c(-0.5,6.5), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 5], col="red")
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 3], col="red", lty=3)
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 7], col="red", lty=3)
tp = d_mgca[order(d_mgca$Age_ka), "Age_ka"]
points(tp, rep(-0.5, nrow(d_mgca)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (ka)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1239,1320), ylim=c(1.6,0))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages, sl$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 5], col="red")
lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 3], col="red", lty=3)
lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 7], col="red", lty=3)
op = d_o[order(d_o$Age_ka),"Age_ka"]
points(op, rep(1.6, nrow(d_o)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

dev.off()


#Shorthand
sl = post.elder$BUGSoutput$sims.list
su = post.elder$BUGSoutput$summary

#Show summary
View(su)

#Get some indicies
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))

##A couple of standard plots of the modeled timeseries
png("T_18O_eld.png", units="in", width=5, height=5, res=300)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1239,1315), ylim=c(-3,4), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 5], col="red")
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 3], col="red", lty=3)
lines(d$ts.ages, su[BWT.start:(BWT.start + d$ts.len - 1), 7], col="red", lty=3)
tp = d$d_mgca[order(d$d_mgca$Age_ka), "Age_ka"]
points(tp, rep(-3, nrow(d$d_mgca)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (ka)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1239,1315), ylim=c(1.0,-0.6))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages, sl$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(d$ts.ages, su[d18O.start:(d18O.start+d$ts.len-1), 5], col="red")
lines(d$ts.ages, su[d18O.start:(d18O.start+d$ts.len-1), 3], col="red", lty=3)
lines(d$ts.ages, su[d18O.start:(d18O.start+d$ts.len-1), 7], col="red", lty=3)
op = d$d_o[order(d$d_o$Age_ka),"Age_ka"]
points(op, rep(1, nrow(d$d_o)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

dev.off()


