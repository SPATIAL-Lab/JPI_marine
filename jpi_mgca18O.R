#####
#JPI for Mg/Ca and d18O using RW timeseries models
#####

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)
library(xlsx)

setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/")
setwd("C:/Users/u0133977/Dropbox/HypoMirror/JPI_marine/")

source("helpers.R")

##Set up timeseries for d18O_sw and BWT modeling
ts.min = 18
ts.max = 0
ts.step = 0.05
ts.ages = seq(ts.min, ts.max, -ts.step)
ts.len = length(ts.ages)

##Prep the foram data, first read
d = read.csv("Lear_combined.csv")

##Now split out the d18O data and strip one outlier
d_o = d[!is.na(d$d18O),]
d_o = d_o[d_o$Sample.ID != "806B 47-5 38-43",]
#Timeseries index for each d18O sample
o_age.ind = round((ts.min - d_o$Age.Ma) / ts.step) + 1

##Now split out the MgCa data and get TS index
d_mgca = d[!is.na(d$MgCa),]
mgca_age.ind = round((ts.min - d_mgca$Age.Ma) / ts.step) + 1

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

##Add indicies for seawater MgCa TS to MgCa foram data
mgca_age.ind.sw = round((mgca_ts.min - d_mgca$Age.Ma) / mgca_ts.step) + 1
mgca_age.ind.all = matrix(c(mgca_age.ind, mgca_age.ind.sw), ncol = 2)

##Read in MgCa calibration dataset
d_mgca_calib = read.csv("mgca_calib.csv")

##Age index for MgCa calibration samples
mgca_calib_age.ind = round((mgca_ts.min - d_mgca_calib$Age) / mgca_ts.step) + 1

##Read in d18O calibration dataset
d_d18O_calib = read.csv("C_comp.csv")
d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "lc", "MgCa_calib.pre", "a", "d18O_calib.pre", "d18O_calib.pre.2",
               "MgCa_sw_m", "MgCa_sw_m.pre", "MgCa_sw_m.eps.ac")

##Data to pass to the model
dat = list(nages = ts.len, nmgca.ages = mgca_ts.len,
           MgCa_calib.bwt.m = d_mgca_calib$BWT, MgCa_calib.bwt.sd = d_mgca_calib$BWT_sd, MgCa_calib = d_mgca_calib$MgCa,
           d18O_calib.bwt.m = d_d18O_calib$Temperature_C, d18O_calib.bwt.sd = rep(0.2,nrow(d_d18O_calib)), d18O_calib = d_d18O_calib$C.SW_d18O,
           MgCa_sw.age.ind = mgca_sw_age.ind, MgCa_sw = d_mgca_sw$MgCa, MgCa_sw.sd = d_mgca_sw$Sigma,
           MgCa.age.ind = mgca_age.ind.all, MgCa = d_mgca$MgCa, 
           d18O.age.ind = o_age.ind, d18O = d_o$d18O)

##Run the inversion
t1 = proc.time()
set.seed(t1[3])
post2 = jags.parallel(model.file = "split_temporal.R", parameters.to.save = parameters, 
             data = dat, n.chains=3, n.iter = 1000000, 
             n.burnin = 100000, n.thin = 100) 
proc.time() - t1

#save(post2, file = "full_post.RData")
#load("full_post.RData")

#Shorthand
sl = post2$BUGSoutput$sims.list
su = post2$BUGSoutput$summary

#Show summary
View(su)

#Get some indicies
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))

##A couple of standard plots of the modeled timeseries
png("T_18O_full.png", units="in", width=5, height=5, res=300)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(0,18), ylim=c(-3,11), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 5], col="red")
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 3], col="red", lty=3)
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 7], col="red", lty=3)
#lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 4], col="red", lty=2)
#lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 6], col="red", lty=2)
tp = d_mgca[order(d_mgca$Age.Ma), "Age.Ma"]
points(tp, rep(-3, nrow(d_mgca)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(0,18), ylim=c(2,-1.5))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages, sl$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 5], col="red")
lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 3], col="red", lty=3)
lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 7], col="red", lty=3)
#lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 4], col="red", lty=2)
#lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 6], col="red", lty=2)
op = d_o[order(d_o$Age.Ma),"Age.Ma"]
points(op, rep(2, nrow(d_o)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

dev.off()

png("MgCa_sw_full.png", units="in", width=5, height=2.75, res=300)
par(mar=c(4,4,1,1), cex=0.85)
plot(-10, 0, xlab="Age (Ma)", ylab ="Seawater Mg/Ca", xlim=c(0,100), ylim=c(0.8,6))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(mgca_ts.ages, sl$MgCa_sw_m[i,], col = rgb(0,0,0, 0.01))
}
lines(mgca_ts.ages, su[MgCa.start:(MgCa.start+mgca_ts.len-1), 5], col="red")
lines(mgca_ts.ages, su[MgCa.start:(MgCa.start+mgca_ts.len-1), 3], col="red", lty=3)
lines(mgca_ts.ages, su[MgCa.start:(MgCa.start+mgca_ts.len-1), 7], col="red", lty=3)
#lines(mgca_ts.ages, su[MgCa.start:(MgCa.start+mgca_ts.len-1), 4], col="red", lty=2)
#lines(mgca_ts.ages, su[MgCa.start:(MgCa.start+mgca_ts.len-1), 6], col="red", lty=2)
points(d_mgca_sw$Age, d_mgca_sw$MgCa, pch=21, bg = "white")
points(d_mgca$Age.Ma, rep(1, nrow(d_mgca)), pch=21, bg = "black")
points(d_mgca_calib$Age, rep(1, nrow(d_mgca_calib)), pch=21, bg = "grey")
dev.off()

##Plotting priors and posteriors

#####
##Big old calibration figure
png("calibration.png", res = 300, units = "in", width = 8, height = 4)
layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 2, ncol = 4, byrow=TRUE))
par(mai=c(0.5,0.5,0.1,0.1))
xoff = 2.3

#MgCa calibration parms
plotd(sl$lc[,1], col="red")
lined(rnorm(100000, 1.5, 0.1))
title(xlab=expression(paste(alpha[1])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

plotd(sl$lc[,2], col="red", ylab="")
lined(rnorm(100000, 0.1, 0.01))
title(xlab=expression(paste(alpha[2])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

plotd(sl$lc[,3], col="red", ylab="")
lined(rnorm(100000, -0.02, 0.03))
title(xlab=expression(paste(alpha[3])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "C")

plotd(sqrt(1/(sl$MgCa_calib.pre)), col="red", xlim=c(0.05,0.25), ylab="")
lined(sqrt(1/(rgamma(100000, 2, 1/30))))
title(xlab=expression(paste(sigma["MgCaf"])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "D")

#d18O calibration parms
plotd(sl$a[,1], col="red")
lined(rnorm(100000, 3.32, 0.02))
title(xlab=expression(paste(beta[1])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "E")

plotd(sl$a[,2], col="red", ylab="")
lined(rnorm(100000, -0.237, 0.01))
title(xlab=expression(paste(beta[2])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "F")

plotd(sl$a[,3], col="red", ylab="")
lined(rnorm(100000, 0.001, 0.0005))
title(xlab=expression(paste(beta[3])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "G")

plotd(sqrt(1/(sl$d18O_calib.pre)), col="red", xlim=c(0.05,0.7), ylab="")
lined(sqrt(1/(sl$d18O_calib.pre.2)), col="red", lty=2)
lined(sqrt(1/(rgamma(100000, 3, 1/30))))
lined(sqrt(1/(rgamma(100000, 6, 1))), lty=2)
title(xlab=expression(paste(sigma[paste(delta, "18Of")])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "H")

dev.off()

#This is the Pleistocene foram d18O sd estimate
plotd(sqrt(1/(sl$d18O_calib.pre.2)), col="red", ylab="", xlim=c(0.25,0.7))
lined(sqrt(1/(rgamma(100000, 6, 1))))
title(xlab=expression(paste(sigma[paste(delta, "18Of")])), line = xoff)

#####
#Let's look at covariance for these parms
png("calibration_covar.png", res=300, units="in", width=6, height=4)
layout(matrix(c(1,2,3,4,5,6), nrow=2, byrow=TRUE))
par(mar=c(5,5,0.2,0.2))

smoothScatter(sl$lc[,1], sl$lc[,2], xlab=expression(paste(alpha[1])), ylab=expression(paste(alpha[2])))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

smoothScatter(sl$lc[,1], sl$lc[,3], xlab=expression(paste(alpha[1])), ylab=expression(paste(alpha[3])))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

smoothScatter(sl$lc[,2], sl$lc[,3], xlab=expression(paste(alpha[2])), ylab=expression(paste(alpha[3])))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "C")

smoothScatter(sl$a[,1], sl$a[,2], xlab=expression(paste(beta[1])), ylab=expression(paste(beta[2])))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "D")

smoothScatter(sl$a[,1], sl$a[,3], xlab=expression(paste(beta[1])), ylab=expression(paste(beta[3])))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "E")

smoothScatter(sl$a[,2], sl$a[,3], xlab=expression(paste(beta[2])), ylab=expression(paste(beta[3])))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "F")

dev.off()

lc.cor = cor(sl$lc)
a.cor = cor(sl$a)

smoothScatter(sl$BWT.pre, sl$BWT.eps.ac)
smoothScatter(sl$d18O_sw.pre, sl$d18O_sw.eps.ac)
smoothScatter(sl$MgCa_sw_m.pre, sl$MgCa_sw_m.eps.ac)

##This is kind of cool, BWT and d18Osw covariation across all sims
png("env_covar.png", res=300, units="in", width=3, height=3)
par(mar=c(5,5,0.5,0.5), cex=0.75)
smoothScatter(sl$BWT, sl$d18O_sw, xlim = c(-3, 9), ylim = c(2, -1.5),
              xlab=expression("BWT ("*degree*" C)"),
              ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"))
lines(c(5.9,5.2,5.9), c(-0.6,0.1,1.1), col="grey")
arrows(5.9,1.1,-0.2,0.35, length=0.05, col="grey")

#this calculates the average correlation of d18O_sw and BWT w/in individual timesteps
BWT.norm = sl$BWT
d18O_sw.norm = sl$d18O_sw
for(j in 1:ts.len){
  BWT.mean = mean(sl$BWT[,j])
  d18O_sw.mean = mean(sl$d18O_sw[,j])
  for(i in 1:sims){
    BWT.norm[i,j] = sl$BWT[i,j] - BWT.mean
    d18O_sw.norm[i,j] = sl$d18O_sw[i,j] - d18O_sw.mean
  }
}
BWT.norm = as.vector(BWT.norm)
d18O_sw.norm = as.vector(d18O_sw.norm)
m = cor(BWT.norm, d18O_sw.norm)

#add it as an arrow
arrows(-2.5, 1.25, -1, 1.25+m*1.15, length=0.05, code = 3)

dev.off()

##Lets look at the change in BWT values as a function of time
#This one calculates change relative to modern
BWT.delta = matrix(rep(0, sims * (ts.len)), nrow = sims, ncol = (ts.len))
for(i in 1:sims){
  for(j in 1:(ts.len)){ BWT.delta[i,j] = sl$BWT[i,j] - sl$BWT[i,ts.len]} 
}

#This gets the zero change value from the emperical CDF of the change time series 
BWT.delta.p = double()
for(j in 1:ts.len){
  tst = ecdf(BWT.delta[,j])
  BWT.delta.p[j] = tst(0)
}

#This gets the CDF value of the modern median from each time steps ECDF 
trad = su[(BWT.start + ts.len - 1), 5]
BWT.p = double()
for(j in 1:ts.len){
  tst = ecdf(sl$BWT[,j])
  BWT.p[j] = tst(trad)
}

#Now a plot showing probabilities on delta T relative to modern
png("deltaT.png", res=300, units="in", width = 5, height = 2.75)
par(mar = c(4,4,1,4), cex = 0.85)
plot(-10, 0, xlab = "", ylab = "", xlim=c(0.05,2), ylim=c(-2.5,2))
title(xlab = "Age (Ma)", line = 2.75)
title(ylab = expression("BWT ("*degree*" C)"), line = 2.75)
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 5], col="red")
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 3], col="red", lty=3)
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 7], col="red", lty=3)
#lines(c(-1,6), rep(trad, 2), col="red", lty=2)
#arrows(1.1, trad, 1.1, -3, length = 0.1, lty=2, col="red")

par(new=TRUE)
plot(ts.ages[1:(ts.len-1)], BWT.delta.p[1:(ts.len-1)], xlim=c(0.05,2), ylim=c(5e-3, 1), type="l", 
     axes=FALSE, log="y", xlab="", ylab="")
lines(ts.ages[1:(ts.len-1)], BWT.p[1:(ts.len-1)], col="red")
lines(c(-1,6), c(0.05,0.05), lty=2)
axis(4, at=c(5e-5,5e-4,5e-3,5e-2,5e-1))
mtext(side = 4, "Zero change probability", line = 2.75)
#arrows(1.02, 5e-2, 1.02, 1e-4, length=0.1, lty=2)

dev.off()

##Same thing now relative to 15.1Ma (ts index 59)
BWT.delta = matrix(rep(0, sims * (ts.len)), nrow = sims, ncol = (ts.len))
for(i in 1:sims){
  for(j in 1:(ts.len)){ BWT.delta[i,j] = sl$BWT[i,j] - sl$BWT[i,59]} 
}

#This gets the zero change value from the emperical CDF of the change time series 
BWT.delta.p = double()
for(j in 1:ts.len){
  tst = ecdf(BWT.delta[,j])
  BWT.delta.p[j] = tst(0)
}

#This gets the CDF value of the 15.1 median from each time steps ECDF 
trad = su[(BWT.start + 59 - 1), 5]
BWT.p = double()
for(j in 1:ts.len){
  tst = ecdf(sl$BWT[,j])
  BWT.p[j] = tst(trad)
}

#Now a plot showing probabilities on delta T relative to 15.1
png("deltaT_Mio.png", res=300, units="in", width = 5, height = 2.75)
par(mar = c(4,4,1,4), cex = 0.85)
plot(-10, 0, xlab = "", ylab = "", xlim=c(14,15.1), ylim=c(4,8))
title(xlab = "Age (Ma)", line = 2.75)
title(ylab = expression("BWT ("*degree*" C)"), line = 2.75)
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 5], col="red")
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 3], col="red", lty=3)
lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 7], col="red", lty=3)
#lines(c(-1,6), rep(trad, 2), col="red", lty=2)
#arrows(1.1, trad, 1.1, -3, length = 0.1, lty=2, col="red")

par(new=TRUE)
plot(ts.ages[60:ts.len], 1-BWT.delta.p[60:ts.len], xlim=c(14,15.1), ylim=c(2e-3, 1), type="l", 
     axes=FALSE, log="y", xlab="", ylab="")
lines(ts.ages[60:ts.len], 1-BWT.p[60:ts.len], col="red")
lines(c(13,16), c(0.05,0.05), lty=2)
axis(4, at=c(5e-5,5e-4,5e-3,5e-2,5e-1))
mtext(side = 4, "Zero change probability", line = 2.75)
#arrows(1.02, 5e-2, 1.02, 1e-4, length=0.1, lty=2)

dev.off()


#####
##Time series properties

#Set it up
png("timeseries.png", res=300, units="in", width = 6, height = 2)
layout(matrix(c(1,2,3), ncol=3))
par(mai = c(0.5,0.5,0.1,0.1))

#BWT and d18O timeseries autocorrelation
plotd(sl$BWT.eps.ac, col="red", ylim=c(0,20))
lined(sl$d18O_sw.eps.ac, col="red", lty=2)
lines(c(0,0.4), c(2.5,2.5))
title(xlab=expression(phi), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

#BWT TS SD
plotd(sqrt(1/sl$BWT.pre), col="red", ylab="")
lined(sqrt(1/rgamma(100000, 20, 2)))
title(xlab=expression(paste(sigma ["BWT"])), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

#d18O TS SD
plotd(sqrt(1/(sl$d18O_sw.pre)), col="red", ylab="", xlim = c(0.075, 0.21), lty=2)
lined(sqrt(1/(rgamma(100000, 10, 1/5))), lty=2)
title(xlab=expression(sigma [paste(delta, "18Osw")]), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "C")

dev.off()

#####
#Covariance

smoothScatter(sl$BWT.pre, sl$BWT.eps.ac)
smoothScatter(sl$d18O_sw.pre, sl$d18O_sw.eps.ac)


#MgCasw timeseries autocorrelation
plotd(sl$MgCa_sw_m.eps.ac, xlab = "TS autocorrelation")

#MgCasw TS SD
plotd(sqrt(1/(sl$MgCa_sw_m.pre)), col="red")
lined(sqrt(1/(rgamma(100000, 1000, 0.01))))

##These are first stabs at plots showing derivatives for BWT, d18O_sw
dv = double()
plot(-1,-1,xlim=c(min(ts.ages),max(ts.ages)), ylim=c(-1.5,1.5))
for(i in seq(1, sims, by = max(floor(sims / 3000),1))){
  for(j in 1:(ts.len - 1)){
    dv[j] = sl$BWT[i,j+1] - sl$BWT[i,j]
  }
  lines(ts.ages[1:(ts.len - 1)], dv, col = rgb(0,0,0, 0.01))
}
for(j in 1:(ts.len - 1)){
  dv[j] = su[BWT.start+j,5] - su[BWT.start+j-1,5]
}
lines(ts.ages[1:(ts.len - 1)], dv, col="red")

dv = double()
plot(-1,-1,xlim=c(min(ts.ages),max(ts.ages)), ylim=c(-0.5,0.5))
for(i in seq(1, sims, by = max(floor(sims / 3000),1))){
  for(j in 1:(ts.len - 1)){
    dv[j] = sl$d18O_sw[i,j+1] - sl$d18O_sw[i,j]
  }
  lines(ts.ages[1:(ts.len - 1)], dv, col = rgb(0,0,0, 0.01))
}
for(j in 1:(ts.len - 1)){
  dv[j] = su[d18O.start+j,5] - su[d18O.start+j-1,5]
}
lines(ts.ages[1:(ts.len - 1)], dv, col="red")

#####
###Now let's try to Shackelton site
d.i = read.table("Birner_2016/datasets/339-U1385_isotope_toRead.tab", sep = "\t", header = TRUE)
d.e = read.table("Birner_2016/datasets/339-U1385_Mg-Ca_toRead.tab", sep = "\t", header = TRUE)
plot(d.e$Age..ka.BP., d.e$U..peregerina.Mg.Ca..mmol.mol.)
plot(d.i$Age..ka.BP., d.i$C..wuellerstorfi.d18O..per.mil.PDB.)

##Set up timeseries for d18O_sw and BWT modeling
ts.min.s = 1320
ts.max.s = 1235
ts.step.s = 1
ts.ages.s = seq(ts.min.s, ts.max.s, -ts.step.s)
ts.len.s = length(ts.ages.s)

#prep the d18O data and add age indicies
d_o.s = d.i[!is.na(d.i$C..wuellerstorfi.d18O..per.mil.PDB.), ]
o_age.ind.s = round((ts.min.s - d_o.s$Age..ka.BP.) / ts.step.s) + 1

#prep the MgCa data and add age indicies
d_mgca.s = d.e[!is.na(d.e$U..peregerina.Mg.Ca..mmol.mol.),]
mgca_age.ind.s = round((ts.min.s - d_mgca.s$Age..ka.BP.) / ts.step.s) + 1

#U spp calibration data from Elderfield 2012 compilation
d_mgca_calib.u = read.csv("u_mgca_calib.csv")

##Read in d18O calibration dataset
d_d18O_calib = read.csv("C_comp.csv")
d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

#get distributions of sw Mg/Ca from long model - not needed for U spp dataset
load("mg_post.RData") #uses output from MgCa_sw_model.R
plotd(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,52:58])
mgca_sw_m.paleo = mean(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,52:58])
mgca_sw_sd.paleo = sd(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,52:58])
mgca_sw_paleo = c(mgca_sw_m.paleo, mgca_sw_sd.paleo)

plotd(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_m.neo = mean(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_sd.neo = sd(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "lc", "MgCa_calib.pre", "a", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = ts.len.s,
           MgCa_calib.bwt.m = d_mgca_calib.u$BWT, MgCa_calib.bwt.sd = rep(0.2, nrow(d_mgca_calib.u)), MgCa_calib = d_mgca_calib.u$MgCa,
           d18O_calib.bwt.m = d_d18O_calib$Temperature_C, d18O_calib.bwt.sd = rep(0.2,nrow(d_d18O_calib)), d18O_calib = d_d18O_calib$C.SW_d18O,
           MgCa_sw.neo = mgca_sw_neo,
           MgCa.age.ind = mgca_age.ind.s, MgCa = d_mgca.s$U..peregerina.Mg.Ca..mmol.mol., 
           d18O.age.ind = o_age.ind.s, d18O = d_o.s$C..wuellerstorfi.d18O..per.mil.PDB.)

##Run the inversion
t1 = proc.time()
set.seed(t1[3])
post = jags.parallel(model.file = "split_temporal_birk.R", parameters.to.save = parameters, 
                      data = dat, n.chains=3, n.iter = 250000, 
                      n.burnin = 5000, n.thin = 25) 
proc.time() - t1

#Shorthand
sl = post$BUGSoutput$sims.list
su = post$BUGSoutput$summary

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
     xlim=c(1230,1320), ylim=c(-1,6), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, sl$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages.s, su[BWT.start:(BWT.start + ts.len.s - 1), 5], col="red")
lines(ts.ages.s, su[BWT.start:(BWT.start + ts.len.s - 1), 3], col="red", lty=3)
lines(ts.ages.s, su[BWT.start:(BWT.start + ts.len.s - 1), 7], col="red", lty=3)
#lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 4], col="red", lty=2)
#lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 6], col="red", lty=2)
tp = d_mgca.s[order(d_mgca.s$Age..ka.BP.), "Age..ka.BP."]
points(tp, rep(-1, nrow(d_mgca.s)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1230,1320), ylim=c(1.5,-0.25))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, sl$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages.s, su[d18O.start:(d18O.start+ts.len.s-1), 5], col="red")
lines(ts.ages.s, su[d18O.start:(d18O.start+ts.len.s-1), 3], col="red", lty=3)
lines(ts.ages.s, su[d18O.start:(d18O.start+ts.len.s-1), 7], col="red", lty=3)
#lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 4], col="red", lty=2)
#lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 6], col="red", lty=2)
op = d_o.s[order(d_o.s$Age..ka.BP.),"Age..ka.BP."]
points(op, rep(1.5, nrow(d_o.s)), pch=21, bg = "white")
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

dev.off()

###Now the Elderfield big record ~ turns out this is too heavy for full simulation

d = read.csv("elderfield_2012.csv")

##Set up timeseries for d18O_sw and BWT modeling
ts.min.s = 1550
ts.max.s = 5
ts.step.s = 5
ts.ages.s = seq(ts.min.s, ts.max.s, -ts.step.s)
ts.len.s = length(ts.ages.s)

#prep the d18O data and add age indicies
d_o.s = d[!is.na(d$d18O), ]
o_age.ind.s = round((ts.min.s - d_o.s$Age_ka) / ts.step.s) + 1

#prep the MgCa data and add age indicies
d_mgca.s = d[!is.na(d$MgCa),] 
mgca_age.ind.s = round((ts.min.s - d_mgca.s$Age_ka) / ts.step.s) + 1

#U spp calibration data from Elderfield 2012 compilation
d_mgca_calib.u = read.csv("u_mgca_calib.csv")

##Read in d18O calibration dataset
d_d18O_calib = read.csv("C_comp.csv")
d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

#get distributions of sw Mg/Ca from long model
load("mg_post.RData") #uses output from MgCa_sw_model.R

mgca_sw_m.neo = mean(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_sd.neo = sd(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.pre", "d18O_sw.eps.ac", "d18O_sw.pre", 
               "lc", "MgCa_calib.pre", "a", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = ts.len.s,
           MgCa_calib.bwt.m = d_mgca_calib.u$BWT, MgCa_calib.bwt.sd = rep(0.2, nrow(d_mgca_calib.u)), MgCa_calib = d_mgca_calib.u$MgCa,
           d18O_calib.bwt.m = d_d18O_calib$Temperature_C, d18O_calib.bwt.sd = rep(0.2,nrow(d_d18O_calib)), d18O_calib = d_d18O_calib$C.SW_d18O,
           MgCa_sw.neo = mgca_sw_neo,
           MgCa.age.ind = mgca_age.ind.s, MgCa = d_mgca.s$MgCa, 
           d18O.age.ind = o_age.ind.s, d18O = d_o.s$d18O)

##Run the inversion
t1 = proc.time()
set.seed(t1[3])
post = jags(model.file = "split_temporal_elderfield.R", parameters.to.save = parameters, 
                     data = dat, n.chains=3, n.iter = 2500, 
                     n.burnin = 500, n.thin = 25) 
proc.time() - t1

#Shorthand
sl = post$BUGSoutput$sims.list
su = post$BUGSoutput$summary

#Show summary
View(su)


#####
###Now Shackelton site and 1123 together
d.b.i = read.table("Birner_2016/datasets/339-U1385_isotope_toRead.tab", sep = "\t", header = TRUE)
d.b.e = read.table("Birner_2016/datasets/339-U1385_Mg-Ca_toRead.tab", sep = "\t", header = TRUE)

d.e = read.csv("elderfield_2012.csv")
d.e = d.e[d.e$Age_ka > 1235,]
d.e = d.e[d.e$Age_ka < 1320,]

##Set up timeseries for d18O_sw and BWT modeling
ts.min.s = 1320
ts.max.s = 1235
ts.step.s = 1
ts.ages.s = seq(ts.min.s, ts.max.s, -ts.step.s)
ts.len.s = length(ts.ages.s)

#prep the d18O data and add age indicies
d_o.b = d.b.i[!is.na(d.b.i$C..wuellerstorfi.d18O..per.mil.PDB.), ]
o_age.ind.b = round((ts.min.s - d_o.b$Age..ka.BP.) / ts.step.s) + 1
d_o.e = d.e[!is.na(d.e$d18O), ]
o_age.ind.e = round((ts.min.s - d_o.e$Age_ka) / ts.step.s) + 1

#prep the MgCa data and add age indicies
d_mgca.b = d.b.e[!is.na(d.b.e$U..peregerina.Mg.Ca..mmol.mol.),]
mgca_age.ind.b = round((ts.min.s - d_mgca.b$Age..ka.BP.) / ts.step.s) + 1
d_mgca.e = d.e[!is.na(d.e$MgCa),] 
mgca_age.ind.e = round((ts.min.s - d_mgca.e$Age_ka) / ts.step.s) + 1

#U spp calibration data from Elderfield 2012 compilation
d_mgca_calib.u = read.csv("u_mgca_calib.csv")

##Read in d18O calibration dataset
d_d18O_calib = read.csv("C_comp.csv")
d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

#get distribution of sw Mg/Ca for early Pleistocene from long model
load("mg_post.RData") #uses output from MgCa_sw_model.R
mgca_sw_m.neo = mean(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_sd.neo = sd(mg_post$BUGSoutput$sims.list$MgCa_sw_m[,110])
mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)

##Parameters to be saved
parameters = c("d18O_sw.b", "BWT.b", "d18O_sw.e", "BWT.e", 
               "BWT.b.eps.ac", "BWT.b.pre", "d18O_sw.b.eps.ac", "d18O_sw.b.pre", 
               "BWT.e.eps.ac", "BWT.e.pre", "d18O_sw.e.eps.ac", "d18O_sw.e.pre", 
               "lc", "MgCa_calib.pre", "a", "d18O_calib.pre")

##Data to pass to the model
dat = list(nages = ts.len.s,
           MgCa_calib.bwt.m = d_mgca_calib.u$BWT, MgCa_calib.bwt.sd = rep(0.2, nrow(d_mgca_calib.u)), MgCa_calib = d_mgca_calib.u$MgCa,
           d18O_calib.bwt.m = d_d18O_calib$Temperature_C, d18O_calib.bwt.sd = rep(0.2,nrow(d_d18O_calib)), d18O_calib = d_d18O_calib$C.SW_d18O,
           MgCa_sw.neo = mgca_sw_neo,
           MgCa.age.ind.b = mgca_age.ind.b, MgCa.b = d_mgca.b$U..peregerina.Mg.Ca..mmol.mol., 
           d18O.age.ind.b = o_age.ind.b, d18O.b = d_o.b$C..wuellerstorfi.d18O..per.mil.PDB.,
           MgCa.age.ind.e = mgca_age.ind.e, MgCa.e = d_mgca.e$MgCa,
           d18O.age.ind.e = o_age.ind.e, d18O.e = d_o.e$d18O)

##Run the inversion ~ 1 hour
t1 = proc.time()
set.seed(t1[3])
post = jags.parallel(model.file = "split_temporal_multi.R", parameters.to.save = parameters, 
                     data = dat, n.chains=3, n.iter = 100000, 
                     n.burnin = 5000, n.thin = 10) 
proc.time() - t1


#Shorthand
sl = post$BUGSoutput$sims.list
su = post$BUGSoutput$summary

#Show summary
View(su)

#Get some indicies
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))

##A couple of standard plots of the modeled timeseries
png("T_18O_multi.png", units="in", width=5, height=5, res=300)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1240,1315), ylim=c(-3,8), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, sl$BWT.b[i,], col = rgb(0.5,0,0, 0.01))
}
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, sl$BWT.e[i,], col = rgb(0,0,0.5, 0.01))
}
lines(ts.ages.s, su[BWT.b.start:(BWT.b.start + ts.len.s - 1), 5], col="red")
lines(ts.ages.s, su[BWT.b.start:(BWT.b.start + ts.len.s - 1), 3], col="red", lty=3)
lines(ts.ages.s, su[BWT.b.start:(BWT.b.start + ts.len.s - 1), 7], col="red", lty=3)
lines(ts.ages.s, su[BWT.e.start:(BWT.e.start + ts.len.s - 1), 5], col=rgb(0.2,0.2,1))
lines(ts.ages.s, su[BWT.e.start:(BWT.e.start + ts.len.s - 1), 3], col=rgb(0.2,0.2,1), lty=3)
lines(ts.ages.s, su[BWT.e.start:(BWT.e.start + ts.len.s - 1), 7], col=rgb(0.2,0.2,1), lty=3)
#lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 4], col="red", lty=2)
#lines(ts.ages, su[BWT.start:(BWT.start + ts.len - 1), 6], col="red", lty=2)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "A")

#Second panel for seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1240,1315), ylim=c(1.75,0))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, sl$d18O_sw.b[i,], col = rgb(0.5,0,0, 0.01))
}
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, sl$d18O_sw.e[i,], col = rgb(0,0,0.5, 0.01))
}
lines(ts.ages.s, su[d18O.b.start:(d18O.b.start+ts.len.s-1), 5], col="red")
lines(ts.ages.s, su[d18O.b.start:(d18O.b.start+ts.len.s-1), 3], col="red", lty=3)
lines(ts.ages.s, su[d18O.b.start:(d18O.b.start+ts.len.s-1), 7], col="red", lty=3)
lines(ts.ages.s, su[d18O.e.start:(d18O.e.start+ts.len.s-1), 5], col=rgb(0.2,0.2,1))
lines(ts.ages.s, su[d18O.e.start:(d18O.e.start+ts.len.s-1), 3], col=rgb(0.2,0.2,1), lty=3)
lines(ts.ages.s, su[d18O.e.start:(d18O.e.start+ts.len.s-1), 7], col=rgb(0.2,0.2,1), lty=3)
#lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 4], col="red", lty=2)
#lines(ts.ages, su[d18O.start:(d18O.start+ts.len-1), 6], col="red", lty=2)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "B")

dev.off()

BWT.delta = sl$BWT.b - sl$BWT.e
BWT.ptiles = matrix(double(), ncol = ncol(BWT.delta), nrow = 4)
for(j in 1:ncol(BWT.delta)){
  BWT.ptiles[1:3,j] = quantile(BWT.delta[,j], c(0.025,0.5,0.975))
  tst = ecdf(BWT.delta[,j])
  BWT.ptiles[4,j] = tst(0)
}

d18O_sw.delta = sl$d18O_sw.b - sl$d18O_sw.e
d18O_sw.ptiles = matrix(double(), ncol = ncol(d18O_sw.delta), nrow = 4)
for(j in 1:ncol(d18O_sw.delta)){
  d18O_sw.ptiles[1:3,j] = quantile(d18O_sw.delta[,j], c(0.025,0.5,0.975))
  tst = ecdf(d18O_sw.delta[,j])
  d18O_sw.ptiles[4,j] = tst(0)
}

plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1240,1315), ylim=c(0,8))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, BWT.delta[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages.s, BWT.ptiles[2,], col="red")
lines(ts.ages.s, BWT.ptiles[1,], col="red", lty=3)
lines(ts.ages.s, BWT.ptiles[3,], col="red", lty=3)
par(new = TRUE)
plot(ts.ages.s, pmax(BWT.ptiles[4,],1e-4), type="l", log="y", axes = FALSE, 
     ylim=c(5e-3,5e-1), xlab="", ylab="")
axis(4, at=c(5e-3,5e-2,5e-1))

plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1240,1315), ylim=c(-1,1))
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(ts.ages.s, d18O_sw.delta[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages.s, d18O_sw.ptiles[2,], col="red")
lines(ts.ages.s, d18O_sw.ptiles[1,], col="red", lty=3)
lines(ts.ages.s, d18O_sw.ptiles[3,], col="red", lty=3)
par(new = TRUE)
plot(ts.ages.s, pmin(d18O_sw.ptiles[4,],1-d18O_sw.ptiles[4,]), type="l", log="y", 
     axes = FALSE, xlim=c(1240,1315), ylim=c(5e-3,5e-1), xlab="", ylab="")
axis(4, at=c(5e-3,5e-2,5e-1))
lines(c(1230,1320), c(0.05,0.05), lty=2)
