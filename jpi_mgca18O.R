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
d_mgca_sw = read.csv("mgca_sw.txt")

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
lines(ts.ages, su[BWT.start:ts.len, 5], col="red")
lines(ts.ages, su[BWT.start:ts.len, 3], col="red", lty=3)
lines(ts.ages, su[BWT.start:ts.len, 7], col="red", lty=3)
#lines(ts.ages, su[BWT.start:ts.len, 4], col="red", lty=2)
#lines(ts.ages, su[BWT.start:ts.len, 6], col="red", lty=2)
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
lines(c(5.9,5.2,5.9,3.3), c(-0.6,0.1,1.1,0.5), col="grey")
arrows(3.3,0.5,-0.2,0.35, length=0.05, col="grey")
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


#save(post2, file = "full_post.RData")

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

###Now let's try to Shackelton site
d.i = read.table("Birner_2016/datasets/339-U1385_isotope_toRead.tab", sep = "\t", header = TRUE)
d.e = read.table("Birner_2016/datasets/339-U1385_Mg-Ca_toRead.tab", sep = "\t", header = TRUE)
plot(d.e$Age..ka.BP., d.e$BWT..Â.C.)
