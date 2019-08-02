#####
##Plots and statistics for JPI manuscript submitted to CoP
#####

#####
#Preliminaries
#####

##My local working directories
setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/code/")
setwd("C:/Users/u0133977/Dropbox/HypoMirror/JPI_marine/code/")

#These are the posterior samples, needed for all of the below
load("post_lear.RData")
load("post_multi.RData")
load("post_birn.RData")
load("post_elder.RData")

#Functions for plotting and data prep, also used throughout
source("helpers.R")

#These first plots use the Lear data, so prep it using the function from helpers.R
d = prep.lear()

##This is the data on which the original interpretations were based
#Data from Lear et al. 2003 and 2015, 'equlibrium offset' of +0.64 removed from
#2003 d18O data. Data only for core levels w/ both Mg/Ca and d18O. d18O values have
#been averaged for levels with multiple observations
dl = read.csv("Lear_combined_interp.csv")
dl_mgca = dl[order(dl$Age.Ma),]
dl_d18O = dl[!is.na(dl$d18O),]

#Shorthand
sl = post.lear$BUGSoutput$sims.list
su = post.lear$BUGSoutput$summary

#Get some indicies used to parse data in su
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))

#####
##Figure 1b: Example segment of posterior sample from 806
#####

#Data frame holding example parameters and variables
ts.ind = round(runif(1, 1, sims))
age.ind = c(99, 120)
ts.df = data.frame(Age = d$ts.ages[age.ind[1]:age.ind[2]], 
                   BWT = sl$BWT[ts.ind, age.ind[1]:age.ind[2]], 
                   BWT_eps = (sl$BWT[ts.ind, age.ind[1]:age.ind[2]] - 
                     sl$BWT[ts.ind, (age.ind[1]-1):(age.ind[2]-1)]))

#Data frame holding MgCa_sw time series, then match these and extrac to ts.df
MgCa_sw.df = data.frame(Age = d$mgca.ages, MgCa_sw = sl$MgCa_sw_m[ts.ind,])
MgCa_sw.merge = MgCa_sw.df[match(ts.df$Age, MgCa_sw.df$Age), 2]
ts.df = cbind(ts.df, MgCa_sw = MgCa_sw.merge)

#Subset of levels w MgCa proxy obs, where MgCa_f will be simulated
ts.df.sub = ts.df[!is.na(ts.df$MgCa_sw),]

#Apply proxy model equation to forward model foram values
ts.df.sub$MgCa_m = sl$a[ts.ind, 1] + sl$a[ts.ind, 2] * ts.df.sub$BWT * ts.df.sub$MgCa_sw ^ sl$a[ts.ind, 3] 
MgCa_sd = sqrt(1/sl$MgCa_calib.pre[ts.ind])

#Set up plots
setEPS()
postscript("../Figure01b.eps", width = 5, height = 7)
layout(matrix(c(1,2,3),1,3), widths = c(lcm(1.5*2.54), lcm(1*2.54), lcm(1.5*2.54)),
       heights = c(lcm(16), lcm(16), lcm(16)))
par(cex = 1.2)

#First panel, epsilon BWT
par(mai=c(0.5,0.5,0,0))
plot(ts.df$BWT_eps, ts.df$Age, ylim=c(14, 13.0), pch=21, bg="white", axes=FALSE)
points(ts.df.sub$BWT_eps, ts.df.sub$Age, pch=21, bg="dark grey")
axis(1)
axis(2)

#Second panel, BWT
par(mai=c(0.5,0,0,0))
plot(ts.df$BWT, ts.df$Age, ylim=c(14, 13.0), pch=21, bg="white", axes=FALSE)
points(ts.df.sub$BWT, ts.df.sub$Age, pch=21, bg="dark grey")
axis(1)

#Third panel, proxy record
plot(0, 0, xlab="Foraminiferal Mg/Ca", ylab="", xlim=c(1.4,2.6), ylim=c(14, 13.0), axes=FALSE)
axis(1)

for(i in nrow(ts.df.sub):1){
  dp = density(rnorm(100000, ts.df.sub$MgCa_m[i], MgCa_sd))
  polygon(c(dp$x, dp$x[1]), ts.df.sub$Age[i] - c(dp$y, dp$y[1])/50, col="light grey", lty=0)
  lines(dp$x, ts.df.sub$Age[i] - dp$y / 50, col="dark grey")
}
points(ts.df.sub$MgCa_m, ts.df.sub$Age, pch=21, bg="dark grey")
points(dl_mgca$MgCa[dl_mgca$Age.Ma>13.1], dl_mgca$Age.Ma[dl_mgca$Age.Ma>13.1], pch=21, bg="red")

dev.off()

#####
##Figure 3: Site 806 BWT and d18Osw time series
#####

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Make space and layout
png("../Figure03.png", units="in", width=5, height=5, res=600)
layout(matrix(c(1,2,3), 3, 1), heights = c(lcm(2.1*2.54), lcm(0.2*2.54), lcm(2.7*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)

##Set up first panel for BWT
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(0,18), ylim=c(-3,9), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()

#Plot 2500 representative JPI posterior samples
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$BWT[i, ts.ind], col = rgb(0,0,0, 0.01))
}

#Plot data distribution
tp = d$d_mgca[order(d$d_mgca$Age.Ma), "Age.Ma"]
points(tp, rep(-3, nrow(d$d_mgca)), pch=21, cex=0.7, bg = "white")

#Recalculate the original interpretations from Lear15 using the LS-NBB calibration equations
#Seawater Mg/Ca based on Lear curve fit, using a +/- 0.5 mol/mol envelope
dl_mgca$mgca_sw = 5.2 - 0.238 * dl_mgca$Age.Ma + 0.00661 * dl_mgca$Age.Ma^2 - 6.66e-5 * dl_mgca$Age.Ma^3
dl_mgca$mgca_sw.max =  dl_mgca$mgca_sw + 0.5
dl_mgca$mgca_sw.min =  dl_mgca$mgca_sw - 0.5

#Four combos of endmember calibration coefficients from table 2 and seawater values 
dl_mgca$BWT = (dl_mgca$MgCa / dl_mgca$mgca_sw ^ (-0.009) - 1.47) / 0.098
dl_mgca$BWT.1 = (dl_mgca$MgCa / dl_mgca$mgca_sw.max ^ (-0.007) - 1.47) / 0.098
dl_mgca$BWT.2 = (dl_mgca$MgCa / dl_mgca$mgca_sw.min ^ (-0.007) - 1.47) / 0.098
dl_mgca$BWT.3 = (dl_mgca$MgCa / dl_mgca$mgca_sw.max ^ (-0.016) - 1.49) / 0.099
dl_mgca$BWT.4 = (dl_mgca$MgCa / dl_mgca$mgca_sw.min ^ (-0.016) - 1.49) / 0.099

#Get the highest and lowest values, add the +/- 1 degree envelope used by L15
dl_mgca$BWT.max = pmax(dl_mgca$BWT.1, dl_mgca$BWT.2, dl_mgca$BWT.3, dl_mgca$BWT.4) + 1
dl_mgca$BWT.min = pmin(dl_mgca$BWT.1, dl_mgca$BWT.2, dl_mgca$BWT.3, dl_mgca$BWT.4) - 1

#Plot the L15 values
lines(dl_mgca$Age.Ma, dl_mgca$BWT, col=rgb(0.2,0.4,1))
lines(dl_mgca$Age.Ma, dl_mgca$BWT.max, col=rgb(0.2,0.4,1), lty=3)
lines(dl_mgca$Age.Ma, dl_mgca$BWT.min, col=rgb(0.2,0.4,1), lty=3)
points(dl_mgca$Age.Ma, dl_mgca$BWT, pch = 20, cex=0.25)

#Median and CIs for JPI
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Second panel with geologic ages
#Setup
par(mai=c(0.05,1,0,0.2))
plot(-10, 0, xlab = "", ylab = "", xlim=c(0,18), ylim=c(0,1), axes = FALSE)

#Add ages
rect(0, par("usr")[3], 2.58, par("usr")[4], col=rgb(0.99,0.99,0.5))
text(2.58/2, (par("usr")[4] - par("usr")[3])/2, "Q", cex=0.8)
rect(2.58, par("usr")[3], 5.33, par("usr")[4], col = rgb(0.9, 0.9, 0.3))
text(2.58 + (5.33-2.58)/2, (par("usr")[4] - par("usr")[3])/2, "Pliocene", cex=0.8)
rect(5.33, par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(0.8, 0.8, 0.1))
text(5.33 + (par("usr")[2]-5.33)/2, (par("usr")[4] - par("usr")[3])/2, "Miocene", cex=0.8)

##Third panel for seawater d18O
#Setup
par(mai=c(1,1,0,0.2))
plot(-10, 0, xlab = "Age (Ma)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(0,18), ylim=c(1.3,-1.5))

#Plot 2500 representative JPI posterior samples
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$d18O_sw[i, ts.ind], col = rgb(0,0,0, 0.01))
}

#Plot data distribution
op = d$d_o[order(d$d_o$Age.Ma),"Age.Ma"]
points(op, rep(1.3, nrow(d$d_o)), pch=21, cex=0.7, bg = "white")

#Calculate d18Osw using original methods
dl_d18O = dl_mgca[!is.na(dl_mgca$d18O),]
dl_d18O$d18Osw = dl_d18O$d18O - (-0.245 * dl_d18O$BWT + 0.0011 * dl_d18O$BWT + 3.58) + 0.27
dl_d18O$d18Osw.max = dl_d18O$d18O - (-0.245 * dl_d18O$BWT.min + 0.0011 * dl_d18O$BWT.min + 3.58) + 0.27
dl_d18O$d18Osw.min = dl_d18O$d18O - (-0.245 * dl_d18O$BWT.max + 0.0011 * dl_d18O$BWT.max + 3.58) + 0.27

#Plot L15 values
lines(dl_d18O$Age.Ma, dl_d18O$d18Osw, col=rgb(0.2,0.4,1))
lines(dl_d18O$Age.Ma, dl_d18O$d18Osw.max, col=rgb(0.2,0.4,1), lty=3)
lines(dl_d18O$Age.Ma, dl_d18O$d18Osw.min, col=rgb(0.2,0.4,1), lty=3)
points(dl_d18O$Age.Ma, dl_d18O$d18Osw, pch=20, cex=0.25)

#Add JPI median and CIs
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

#####
#Figure 2: The Mg/Ca sw time series
#####

#Get index values for base time series
ts.ind = match(d$mgca.ts.ages.base, d$mgca.ages)

#Make space and layout
png("../Figure02.png", units="in", width=5, height=2.75, res=600)
par(mar=c(4,4,1,1), cex=0.85)

##Panel setup
plot(-10, 0, xlab="Age (Ma)", ylab ="Seawater Mg/Ca", xlim=c(0,80), ylim=c(0.8,5.5))

#Plot 2500 representative JPI posterior samples
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$mgca.ages[ts.ind], sl$MgCa_sw_m[i,ts.ind], col = rgb(0,0,0, 0.01))
}

#curve fit from L15 for comparison
ages = seq(0, 55, 1)
vals = 5.2 - 0.238 * ages + 0.00661 * ages^2 - 6.66e-5 * ages^3
lines(ages, vals, col=rgb(0.2,0.4,1))

#Add JPI median and 95% CIs
lines(d$mgca.ages[ts.ind], su[MgCa.start + ts.ind, 5], col="red")
lines(d$mgca.ages[ts.ind], su[MgCa.start + ts.ind, 3], col="red", lty=3)
lines(d$mgca.ages[ts.ind], su[MgCa.start + ts.ind, 7], col="red", lty=3)

#Points showing Mg/Ca proxy obs and distribution of proxy and calib data
points(d$d_mgca_sw$Age, d$d_mgca_sw$MgCa, pch=21, bg = "white")
points(d$d_mgca$Age.Ma, rep(0.8, nrow(d$d_mgca)), pch=21, bg = "black")
calib_ages = d$d_mgca_calib$Age
calib_ages = calib_ages[calib_ages>0]
points(calib_ages, rep(0.8, length(calib_ages)), pch=21, bg = "grey")

dev.off()

######
#Figure 5: Prior/posterior plots for calibration parameters
#####

#Make space and layout
png("../Figure05.png", res = 600, units = "in", width = 8, height = 4)
layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 2, ncol = 4, byrow=TRUE))
par(mai=c(0.5,0.5,0.1,0.1))
xoff = 2.3  #Used in positioning axis labels

##Panel 1 - MgCa parameter alpha 1
#Posterior
plotd(sl$a[,1], col="red") 

#Prior
lined(rnorm(100000, 1.5, 0.1)) 

#Axis and panel labels
title(xlab=expression(paste(alpha[1])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Panel 2 - MgCa parameter alpha 2
plotd(sl$a[,2], col="red", ylab="")
lined(rnorm(100000, 0.1, 0.01))
title(xlab=expression(paste(alpha[2])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

##Panel 3 - MgCa parameter alpha 3
plotd(sl$a[,3], col="red", ylab="")
lined(rnorm(100000, -0.02, 0.03))
title(xlab=expression(paste(alpha[3])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(c)")

##Panel 4 - MgCa precision parameter
plotd(sqrt(1/(sl$MgCa_calib.pre)), col="red", xlim=c(0.05,0.25), ylab="")
lined(sqrt(1/(rgamma(100000, 2, 1/30))))
title(xlab=expression(paste(sigma["MgCaf"])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(d)")

##Panel 5 - d18O parameter beta 1
plotd(sl$b[,1], col="red")
lined(rnorm(100000, 3.32, 0.02))
title(xlab=expression(paste(beta[1])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(e)")

##Panel 6 - d18O parameter beta 2
plotd(sl$b[,2], col="red", ylab="")
lined(rnorm(100000, -0.237, 0.01))
title(xlab=expression(paste(beta[2])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(f)")

##Panel 7 - d18O parameter beta 3
plotd(sl$b[,3], col="red", ylab="")
lined(rnorm(100000, 0.001, 0.0005))
title(xlab=expression(paste(beta[3])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(g)")

##Panel 8 - d18O precision parameter
#Posterior for pre-800 ka
plotd(sqrt(1/(sl$d18O_calib.pre)), col="red", xlim=c(0.05,0.7), ylab="")

#Posterior for post-800 ka
lined(sqrt(1/(sl$d18O_calib.pre.2)), col="red", lty=2)

#Prior for pre-800 ka
lined(sqrt(1/(rgamma(100000, 3, 1/30))))

#Prior for post-800 ka
lined(sqrt(1/(rgamma(100000, 6, 1))), lty=2)
title(xlab=expression(paste(sigma[paste(delta, "18Of")])), line = xoff)
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(h)")

dev.off()

#####
#Figure 6: Parameter covariance
#####

#Make space and layout
png("../Figure06.png", res=600, units="in", width=6, height=4)
layout(matrix(c(1,2,3,4,5,6), nrow=2, byrow=TRUE))
par(mar=c(4,4,0.4,0.4))

##Panel 1 - MgCa parameters alpha 1 and 2
#Density plot
smoothScatter(sl$a[,1], sl$a[,2], xlab="", ylab="")

#Panel label
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

#Axis titles
title(xlab=expression(paste(alpha[1])), line = 2.5)
title(ylab=expression(paste(alpha[2])), line = 2.5)

##Panel 2 - MgCa parameters alpha 1 and 3
smoothScatter(sl$a[,1], sl$a[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")
title(xlab=expression(paste(alpha[1])), line = 2.5)
title(ylab=expression(paste(alpha[3])), line = 2.5)

##Panel 3 - MgCa parameters alpha 2 and 3
smoothScatter(sl$a[,2], sl$a[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(c)")
title(xlab=expression(paste(alpha[2])), line = 2.5)
title(ylab=expression(paste(alpha[3])), line = 2.5)

##Panel 4 - d18O parameters beta 1 and 2
smoothScatter(sl$b[,1], sl$b[,2], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(d)")
title(xlab=expression(paste(beta[1])), line = 2.5)
title(ylab=expression(paste(beta[2])), line = 2.5)

##Panel 5 - d18O parameters beta 1 and 3
smoothScatter(sl$b[,1], sl$b[,3], xlab="", ylab="")
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(e)")
title(xlab=expression(paste(beta[1])), line = 2.5)
title(ylab=expression(paste(beta[3])), line = 2.5)

##Panel 6 - d18O parameters beta 2 and 3
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

#####
#Figure 9: 2-d posterior density of environmental time series parameters 
#####

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Subset the posterior to retain only base ts steps
BWT = sl$BWT[, ts.ind]
d18Osw = sl$d18O_sw[, ts.ind]

#First calculate the difference of each BWT and d18O_sw value relative to
#the 18Ma value for that posterior draw
#Set aside space
nr = nrow(BWT)
nc = ncol(BWT)
D_BWT = matrix(double(), nrow = nr, ncol = nc)
D_d18O_sw = matrix(double(), nrow = nr, ncol = nc)

#Calculate differences
for(i in 1:nr){
  for(j in  1:nc){
    D_BWT[i,j] = BWT[i,j] - BWT[i,1]
    D_d18O_sw[i,j] = d18Osw[i,j] - d18Osw[i,1]
  }
}

#Now calculate the medians of the values obtained above
#Set aside space
D_BWT.m = double()
D_d18O_sw.m = double()

#Get the medians
for(i in 1:nc){
  D_BWT.m[i] = median(D_BWT[,i])
  D_d18O_sw.m[i] = median(D_d18O_sw[,i])
}

#Matrix gets aweful big so subsample prior to KDE
subs = seq(1, length(D_BWT), by=floor(length(D_BWT)/500000))

#2-d kernal density for the differences
library(MASS)
DKE = kde2d(D_BWT[subs], D_d18O_sw[subs], h=c(0.75, 0.3), n=100)

#Correlation for different periods
lm15 = lm(D_d18O_sw.m[d$ts.ages.base > 15] ~ D_BWT.m[d$ts.ages.base > 15])
lm6 =  lm(D_d18O_sw.m[d$ts.ages.base > 6 & d$ts.ages.base < 14] ~ 
            D_BWT.m[d$ts.ages.base > 6 & d$ts.ages.base < 14])

#Make space and layout
png("../Figure09.png", res=600, units="in", width=3, height=3)
par(mar=c(5,5,0.5,0.5), cex=0.75)

##Density plot
smoothScatter(D_BWT, D_d18O_sw, xlab=expression(Delta*"BWT ("*degree*" C)"),
              ylab = expression(Delta*delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
              xlim = c(-9,4.5), ylim = c(-1,2), col="white")

#Add contours
contour(DKE, add=TRUE, drawlabels=FALSE, col="grey")

#Colors for plotting the medians
pal = heat.colors(nc)

#Lines showing correlation of medians
#abline(lm15, lty=3, col="grey")
#abline(lm6, lty=3, col="grey")

#Plot the medians
points(D_BWT.m, D_d18O_sw.m, pch=19, col=pal, cex=0.20)

#Legend for the median colorscale
rect(-9.2, 2.05, -3.6, 1.55, col="grey")
points(seq(-8.9,-3.9,length.out = nc), rep(1.65, nc), col=pal, cex=0.6)
a = c(18,12,6,0)
xs = seq(-8.8, -3.9, length.out = 4)
text(xs, 1.75, a, cex=0.7)
text(-6.4, 1.93, "Age (Ma)", cex=0.7)

dev.off()

#####
#Figure 7: 9-panel posterior density plots for time series model parameters
#####

#Make space and layout
png("../Figure07.png", res=600, units="in", width = 5.6, height = 5.2)
layout(matrix(seq(1,9), nrow=3, byrow = TRUE), 
       widths=c(lcm(2*2.54), lcm(1.8*2.54), lcm(1.8*2.54)),
       heights=c(lcm(1.6*2.54), lcm(1.6*2.54), lcm(2*2.54)))

#Load site 806 results - this plot uses both site 806 and multi results
sl = post.lear$BUGSoutput$sims.list

##Panel 1 - BWT and d18O time series error autocorrelation
#Margins
par(mai = c(0.1,0.5,0.1,0.1))

#Plot posterior for BWT
plotd(sl$BWT.eps.ac, col="red", ylab="", ylim=c(0,25), xlim=c(0,1), axes=FALSE)

#Add axes
axis(1, labels = FALSE)
axis(2)
box()

#Plot posterior for d18Osw
lined(sl$d18O_sw.eps.ac, col="red", lty=2)

#Add prior
lines(c(0,1), c(1,1))

#Panel label
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Panel 2 - BWT time series error standard deviation
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/sl$BWT.pre), col="red", ylab="", xlim=c(0.04,0.14), ylim=c(0,60), axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/rgamma(100000, 20, 0.1)))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

##Panel 3 - d18O time series error standard deviation
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/(sl$d18O_sw.pre)), col="red", ylab="", xlim = c(0.012, 0.03), ylim=c(0,250), lty=2, axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/(rgamma(100000, 30, 0.01))), lty=2)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(c)")

#Load multi-site 1123 and U1385 results
sl = post.multi$BUGSoutput$sims.list

## Panel 4 - BWT and d18O time series error autocorrelation, site U1385
par(mai = c(0.1,0.5,0.1,0.1))
plotd(sl$BWT.b.eps.ac, col="red", ylim=c(0,30), xlim=c(0,1), axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sl$d18O_sw.b.eps.ac, col="red", lty=2)
lines(c(0,1), c(1,1))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(d)")

##Panel 5 - BWT time series error standard deviation, site U1385
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/sl$BWT.b.pre), col="red", ylab="", xlim=c(0.04,0.14), axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/rgamma(100000, 20, 0.1)))
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(e)")

##Panel 6 - d18O time series error standard deviation, site U1385
par(mai = c(0.1,0.3,0.1,0.1))
plotd(sqrt(1/(sl$d18O_sw.b.pre)), col="red", ylab="", xlim = c(0.012, 0.03), lty=2, axes=FALSE)
axis(1, labels=FALSE)
axis(2)
box()
lined(sqrt(1/(rgamma(100000, 30, 0.01))), lty=2)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(f)")

##Panel 7 - BWT and d18O time series error autocorrelation, site 1123
par(mai = c(0.5,0.5,0.1,0.1))
plotd(sl$BWT.e.eps.ac, col="red", ylab="", ylim=c(0,30), xlim=c(0,1))
lined(sl$d18O_sw.e.eps.ac, col="red", lty=2)
lines(c(0,1), c(1,1))
title(xlab=expression(phi), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(g)")

##Panel 8 - BWT time series error standard deviation, site 1123
par(mai = c(0.5,0.3,0.1,0.1))
plotd(sqrt(1/sl$BWT.e.pre), col="red", ylab="", xlim=c(0.04,0.14))
lined(sqrt(1/rgamma(100000, 20, 0.1)))
title(xlab=expression(paste(sigma ["BWT"])), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(h)")

##Panel 9 - d18O time series error standard deviation, site 1123
par(mai = c(0.5,0.3,0.1,0.1))
plotd(sqrt(1/(sl$d18O_sw.e.pre)), col="red", ylab="", xlim = c(0.012, 0.03), lty=2)
lined(sqrt(1/(rgamma(100000, 30, 0.01))), lty=2)
title(xlab=expression(sigma [paste(delta, "18Osw")]), line=xoff)
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(i)")

dev.off()

#####
#Figure 4: Site 1123 and U1385 posterior time series
#####

#These use data from multi run, load it
d = prep.multi()

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Get and recalibrate original interpreted records using our MgCa slope (0.068)
#These data are directly from the authors' SI and include both raw proxy measurements
#and their values interpreted using the Elderfield et al (2010) down-core calibration
#Only levels w/ both Mg/Ca and d18O data have been retained
db = read.csv("birner_2016_interp.csv")
db$BWT.recal = (db$MgCa - 1) / 0.068
db$d18O_sw.recal = db$d18O - 3.95 - (-0.25 * db$BWT.recal)
de = read.csv("elderfield_2012_interp.csv")
de$BWT.recal = (de$MgCa - 1) / 0.068
de$d18O_sw.recal = de$d18O - 3.95 - (-0.25 * de$BWT.recal)

#Shorthand
sl = post.multi$BUGSoutput$sims.list
su = post.multi$BUGSoutput$summary

#Get some indicies for parsing data in su
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))

#Make space and layout
png("../Figure04.png", units="in", width=5, height=5, res=600)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)

##Panel 1 showing BWT
#Setup
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1239,1315), ylim=c(-3,8), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()

#Plot 2500 representative samples from JPI posterior, first site U1385
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$BWT.b[i, ts.ind], col = rgb(0.5,0,0, 0.01))
}

#Now site 1123
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$BWT.e[i, ts.ind], col = rgb(0,0,1, 0.01))
}

#Add reconstructions from original papers 
#U1385 95% CI based on standard deviation reported in the paper
arrows(1279.45, 5.16-1.74, 1279.45, 5.16+1.74, code=3, angle=90, length=0.06, col="red")

#U1385 re-calibrated published record
lines(db$Age_ka, db$BWT.recal, col="red")
points(db$Age_ka, db$BWT.recal, pch=20, cex=0.25)

#1123 95% CI and re-calibrated published record
arrows(1273.7, -0.79-2, 1273.7, -0.79+2, code=3, angle=90, length=0.06, col=rgb(0,0,0.7))
lines(de$Age.ka, de$BWT.recal, col=rgb(0,0,0.7))
points(de$Age.ka, de$BWT.recal, pch=20, cex=0.25)

#Add JPI medians and 95% CIs
lines(d$ts.ages.base, su[BWT.b.start + ts.ind - 1, 5])
lines(d$ts.ages.base, su[BWT.b.start + ts.ind - 1, 3], lty=3)
lines(d$ts.ages.base, su[BWT.b.start + ts.ind - 1, 7], lty=3)
lines(d$ts.ages.base, su[BWT.e.start + ts.ind - 1, 5])
lines(d$ts.ages.base, su[BWT.e.start + ts.ind - 1, 3], lty=3)
lines(d$ts.ages.base, su[BWT.e.start + ts.ind - 1, 7], lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Panel 2 - seawater d18O
#Setup
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (ka)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1239,1315), ylim=c(1.5,-0.6))

#Plot 2500 representative samples from JPI posteriors
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$d18O_sw.b[i, ts.ind], col = rgb(0.5,0,0, 0.01))
}
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$d18O_sw.e[i, ts.ind], col = rgb(0,0,1, 0.01))
}

#Add reconstructions from original papers
#U1385 95% CI
arrows(1279.45, 0.87-0.46, 1279.45, 0.87+0.46, code=3, angle=90, length=0.06, col="red")

#U1385 record
lines(db$Age_ka, db$d18O_sw.recal, col="red")
points(db$Age_ka, db$d18O_sw.recal, pch=20, cex=0.25)

#1123 95% CI
arrows(1273.7, -0.14-0.4, 1273.7, -0.14+0.4, code=3, angle=90, length=0.06, col=rgb(0,0,0.7))

#1123 record
lines(de$Age.ka, de$d18O_sw.recal, col=rgb(0,0,0.7))
points(de$Age.ka, de$d18O_sw.recal, pch=20, cex=0.25)

#Add JPI medians and 95% CIs
lines(d$ts.ages.base, su[d18O.b.start + ts.ind - 1, 5])
lines(d$ts.ages.base, su[d18O.b.start + ts.ind - 1, 3], lty=3)
lines(d$ts.ages.base, su[d18O.b.start + ts.ind - 1, 7], lty=3)
lines(d$ts.ages.base, su[d18O.e.start + ts.ind - 1, 5])
lines(d$ts.ages.base, su[d18O.e.start + ts.ind - 1, 3], lty=3)
lines(d$ts.ages.base, su[d18O.e.start + ts.ind - 1, 7], lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

#####
#Figure 8: Statistical inference with JPI
#####

#Using site 806 analysis, load it
d = prep.lear()

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Shorthand
sl = post.lear$BUGSoutput$sims.list
su = post.lear$BUGSoutput$summary

#Get some indicies used for parsing data in su
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))

#Calculates BWT change relative to modern in each posterior sample
#Make space
BWT.delta = matrix(rep(0, sims * (length(ts.ind))), nrow = sims)

#Calculate the differences
for(i in 1:sims){
  for(j in 1:length(ts.ind)){ BWT.delta[i,j] = sl$BWT[i,ts.ind[j]] - sl$BWT[i,d$ts.len]} 
}

#Now get the zero change value from the emperical CDF of the change time series 
#Make space
BWT.delta.p = double()

#Calculate the CDF and find quantile of zero value in it
for(j in 1:length(ts.ind)){
  tst = ecdf(BWT.delta[,j]) #This creates a function representing the CDF
  BWT.delta.p[j] = tst(0) #This gets the zero value
}

#Get the CDF value of the modern median from the JPI posterior sample 'slice'
#at each time step. This represents the classical comparison, made across rather
#than within the posterior samples
#First get the modern (Age zero) median
trad = su[(BWT.start + d$ts.len - 1), 5]

#Make space
BWT.p = double()

#Calculate the CDF and find quantile value of modern median in it
for(j in 1:length(ts.ind)){
  tst = ecdf(sl$BWT[,ts.ind[j]])
  BWT.p[j] = tst(trad)
}

#Make space and layout for the figure
png("../Figure08.png", res=600, units="in", width = 5, height = 5.5)
layout(matrix(c(1,2), nrow = 2))
par(mar = c(4,4.5,1,4), cex = 0.85)

##Panel 1 showing the difference relative to modern results
#Setup
plot(-10, 0, xlab = "", ylab = "", xlim=c(0.05,2), ylim=c(-2.5,2))
title(xlab = "Age (Ma)", line = 2.75)
title(ylab = expression("BWT ("*degree*" C)"), line = 2.75)

#Plot 2500 representative samples of the JPI posterior for reference
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$BWT[i, ts.ind], col = rgb(0,0,0, 0.01))
}

#Add JPI median and 95% CIs
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

#Now add the CDF values on the right-hand y-axis
#Setup
par(new=TRUE)

#First plot the within-samples result
plot(d$ts.ages.base[1:length(ts.ind)-1], BWT.delta.p[1:length(ts.ind)-1], xlim=c(0.05,2), ylim=c(5e-3, 1), type="l", 
     axes=FALSE, log="y", xlab="", ylab="", col="blue")

#Add the between-samples result
lines(d$ts.ages.base[1:length(ts.ind)-1], BWT.p[1:length(ts.ind)-1], lty=3, col="blue")

#Lines showing exceedence probabilities of 95 and 99%
lines(c(-1,6), c(0.05,0.05), lty=2, col="blue")
lines(c(-1,6), c(0.01,0.01), lty=2, col="blue")

#Add second y-axis and label
axis(4, at=c(5e-5,5e-4,5e-3,5e-2,5e-1))
mtext(side = 4, "Zero change probability", line = 2.75)

##Now switch over to the multi-site 1123 and U1385 analysis, load the results
d = prep.multi()

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Shorthand
sl = post.multi$BUGSoutput$sims.list
su = post.multi$BUGSoutput$summary

#Get some indicies used in parsing data from su
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))

#Calculate within-sample difference between d18O_sw for the two records
d18O_sw.delta = sl$d18O_sw.b[,ts.ind] - sl$d18O_sw.e[,ts.ind]

#Get quantile value for zero difference
#Make space
d18O_sw.ptiles = matrix(double(), ncol = ncol(d18O_sw.delta), nrow = 4)

#Calculate the difference and get the quantile value for zero difference 
for(j in 1:ncol(d18O_sw.delta)){
  d18O_sw.ptiles[1:3,j] = quantile(d18O_sw.delta[,j], c(0.025,0.5,0.975))
  tst = ecdf(d18O_sw.delta[,j])
  d18O_sw.ptiles[4,j] = tst(0)
}

##Panel 2 - Difference in d18Osw and zero difference probs
#Make space and layout
plot(-10, 0, xlab = "Age (ka)", ylab = expression(Delta*delta^{18}*"O"[sw]*" (U1385 - 1123)"), 
     xlim=c(1239,1315), ylim=c(-0.5,2))

#Plot 1500 representative samples from the JPI posterior
for(i in seq(1, sims, by = max(floor(sims / 1500),1))){
  lines(d$ts.ages.base, d18O_sw.delta[i,], col = rgb(0,0,0, 0.01))
}

#Add JPI median and 95% CIs
lines(d$ts.ages.base, d18O_sw.ptiles[2,], col="red")
lines(d$ts.ages.base, d18O_sw.ptiles[1,], col="red", lty=3)
lines(d$ts.ages.base, d18O_sw.ptiles[3,], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

#Now add the CDF values on the right-hand y-axis
#Setup
par(new = TRUE)

#First plot the within-sample zero difference probabilites
#Use of nested min/max - min gives 2-sided test (e.g., 0.99 plots as 0.01),
#max plots off-scale values as flat line at bottom of panel
plot(d$ts.ages.base, pmax(pmin(d18O_sw.ptiles[4,],1-d18O_sw.ptiles[4,]), 5e-4), type="l", 
     log="y", axes = FALSE, xlim=c(1239,1315), ylim=c(5e-4,5e-1), xlab="", ylab="", 
     col="blue")

#Calculate among-sample difference between d18O_sw for two records
#This is based on the independent U1385 and 1123 JPI results, so comparisons
#are made between random samples

d.b = prep.birn()
d.e = prep.elder()

#Get index values for base time series
ts.ind.b = match(d.b$ts.ages.base, d.b$ts.ages)
ts.ind.e = match(d.e$ts.ages.base, d.e$ts.ages)

#Calculate the differences
d18O_sw.delta = post.birn$BUGSoutput$sims.list$d18O_sw[,ts.ind.b] - 
  post.elder$BUGSoutput$sims.list$d18O_sw[,ts.ind.e]

#Get quantile value for zero difference
#Make space
d18O_sw.ptiles = matrix(double(), ncol = ncol(d18O_sw.delta), nrow = 4)

#Calculate the quantile values
for(j in 1:ncol(d18O_sw.delta)){
  d18O_sw.ptiles[1:3,j] = quantile(d18O_sw.delta[,j], c(0.025,0.5,0.975))
  tst = ecdf(d18O_sw.delta[,j])
  d18O_sw.ptiles[4,j] = tst(0)
}

#Plot the results
lines(d$ts.ages.base, pmax(pmin(d18O_sw.ptiles[4,],1-d18O_sw.ptiles[4,]), 5e-4), lty=3, 
      col="blue")

#Lines showing exceedence probabilities of 95 and 99%
lines(c(1230,1320), c(0.05,0.05), lty=2, col="blue")
lines(c(1230,1320), c(0.01,0.01), lty=2, col="blue")

#Add second y-axis and title
axis(4, at=c(5e-4,5e-3,5e-2,5e-1))
mtext(side = 4, "Zero difference probability", line = 2.75)

dev.off()

#####
#Supplemental Figure 1: Constraints of down-core and core-top data 
#on Uvigerina Mg/Ca T sensitivity, after Elderfield et al. (2010)
#####

#Data needed for plotting
load("post_mgca_dc.RData")
load("post_mgca_ct.RData")

#Make space and layout
png("../SI_Figure01.png", res=600, units="in", width=3, height=3)
par(mar=c(5,5,0.5,0.5), cex=0.75)

##Plot prior
plotd(runif(100000, -0.05, 0.35), ylim=c(0,100), 
      xlab = expression(italic("Uvigerina") * " Mg/Ca T sensitivity (" * alpha[2] * ")"))

#Add posterior distribution from down-core analysis
lined(post.mgca.dc$BUGSoutput$sims.list$a[,2], col="red") #posterior

#Add posterior distribution  from core-top analysis
lined(post.mgca.ct$BUGSoutput$sims.list$a[,2], col="red", lty=2)

dev.off()

#####
#Supplementary Figure 2: Results of single-site JPI for site U1385
#####

#Load data for site U1385
d = prep.birn()

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Original interpreted data to plot for comparison
db = read.csv("birner_2016_interp.csv")

#Shorthand
sl = post.birn$BUGSoutput$sims.list
su = post.birn$BUGSoutput$summary

#Get some indicies used to parse data in su
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))

#Make space and layout
png("../SI_Figure02.png", units="in", width=5, height=5, res=600)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)

##Set up panel 1 - BWT
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1239,1315), ylim=c(-0.5,7), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()

#Plot 2500 representative samples from JPI posterior
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$BWT[i,ts.ind], col = rgb(0,0,0, 0.01))
}

#Add original interpretation of U1385 record based on down-core calibration 
lines(db$Age_ka, db$BWT, col=rgb(0,0,0.5))

#Add JPI median and 95% CIs
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Set up panel 2 - seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (ka)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1239,1315), ylim=c(1.5, -0.2))

#Plot 2500 representative samples from JPI posterior
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$d18O_sw[i, ts.ind], col = rgb(0,0,0, 0.01))
}

#Add original authors' interpretation
lines(db$Age_ka, db$d18O_sw, col=rgb(0,0,0.5))

#Add JPI median and 95% CIs
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

#####
#Supplementary Figure 3: Results of single-site JPI for site 1123
#####

#Load data for site 1123
d = prep.elder()

#Get index values for base time series
ts.ind = match(d$ts.ages.base, d$ts.ages)

#Original interpreted data to plot for comparison
de = read.csv("elderfield_2012_interp.csv")

#Shorthand
sl = post.elder$BUGSoutput$sims.list
su = post.elder$BUGSoutput$summary

#Get some indicies used to parse data in su
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))

#Make space and layout
png("../SI_Figure03.png", units="in", width=5, height=5, res=600)
layout(matrix(c(1,2), 2, 1), heights = c(lcm(2.1*2.54), lcm(2.9*2.54)))
par(mai=c(0.2,1,0.2,0.2), cex=0.85)

##Set up panel 1 - BWT
plot(-10, 0, xlab = "", ylab = expression("BWT ("*degree*" C)"),
     xlim=c(1239,1315), ylim=c(-3,4), axes = FALSE)
axis(1, labels=FALSE)
axis(2)
box()

#Plot 2500 represenatative samples from JPI posterior
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$BWT[i, ts.ind], col = rgb(0,0,0, 0.01))
}

#Add original authors' interpretation based on donw-core calibration
lines(de$Age.ka, de$BWT, col=rgb(0,0,0.5))

#Add JPI median and 95% CIs
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[BWT.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Set up second panel - seawater d18O
par(mai=c(1,1,0.2,0.2))
plot(-10, 0, xlab = "Age (ka)", ylab = expression(delta^{18}*"O"[sw]*" (\u2030, VSMOW)"), 
     xlim=c(1239,1315), ylim=c(1.0,-0.6))

#Plot 2500 represenatative samples from JPI posterior
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  lines(d$ts.ages.base, sl$d18O_sw[i, ts.ind], col = rgb(0,0,0, 0.01))
}

#Add original authors' interpretation
lines(de$Age.ka, de$d18O_sw, col=rgb(0,0,0.5))

#Add JPI median and 95% CIs
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 5], col="red")
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 3], col="red", lty=3)
lines(d$ts.ages.base, su[d18O.start + ts.ind - 1, 7], col="red", lty=3)

#Panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/25
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

#####
#Supplementary Figure 4: Site 806 calibration data and posterior draws
#####

d = prep.lear()

#Shorthand
sl = post.lear$BUGSoutput$sims.list
sims = nrow(sl$BWT)

#Make space and layout
png("../SI_Figure04.png", units="in", width=8, height=4, res=600)
layout(matrix(c(1,2), ncol=2), widths = c(4,4), heights = c(4,4))
par(mai=c(1,1,0.2,0.2), cex=0.85)

##Panel 1: Mg/Ca
#Set up plot
plot(d$d_mgca_calib$BWT, d$d_mgca_calib$MgCa, xlab = expression("BWT ("*degree*" C)"), 
     ylab = expression("Mg/Ca"[f]))

#Sequence of BWT values for plotting 
bwts = seq(-5, 20, by=2)

#Plot 2500 represenatative samples from JPI posterior
#Used fixed value of 3.5 for Mg/Ca_sw
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  mgcas = (sl$a[i,1] + sl$a[i,2] * bwts) * 3.5 ^ sl$a[i,3]
  lines(bwts, mgcas, col = rgb(0,0,0, 0.01))
}

#Calculate median model parameters
a1 = median(sl$a[,1])
a2 = median(sl$a[,2])
a3 = median(sl$a[,3])

#Add median values showing Mg/Ca sensitivity
lines(bwts, (a1 + a2 * bwts) * 3.5 ^ a3, col="red")
lines(bwts, (a1 + a2 * bwts) * 5.5 ^ a3, col="red", lty=2)
lines(bwts, (a1 + a2 * bwts) * 1.5 ^ a3, col="red", lty=3)

points(d$d_mgca_calib$BWT, d$d_mgca_calib$MgCa, pch=21, bg="white")

#panel label
xl = par("usr")[1]+(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(a)")

##Panel 2: d18O
#Set up plot
plot(d$d_d18O_calib$BWT, d$d_d18O_calib$d18O_f.sw, 
     xlab = expression("BWT ("*degree*" C)"), 
     ylab = expression(Delta*delta^{18}*"O"[f]*" (\u2030)"))

#Plot 2500 represenatative samples from JPI posterior
#Used fixed value of 3.5 for Mg/Ca_sw
for(i in seq(1, sims, by = max(floor(sims / 2500),1))){
  mgcas = sl$b[i,1] + sl$b[i,2] * bwts + sl$b[i,3]* bwts^2
  lines(bwts, mgcas, col = rgb(0,0,0, 0.01))
}

#Calculate median model parameters
b1 = median(sl$b[,1])
b2 = median(sl$b[,2])
b3 = median(sl$b[,3])

#Add median values showing Mg/Ca sensitivity
lines(bwts, b1 + b2 * bwts + b3 * bwts^2, col="red")

#replot data on top
points(d$d_d18O_calib$BWT, d$d_d18O_calib$d18O_f.sw, pch=21, bg="white")

#panel label
xl = par("usr")[2]-(par("usr")[2]-par("usr")[1])/15
yl = par("usr")[4]-(par("usr")[4]-par("usr")[3])/15
text(xl, yl, "(b)")

dev.off()

#####
#Calculate 95% CI width for environmental records
#####

##Site 806

#Shorthand
sl = post.lear$BUGSoutput$sims.list
su = post.lear$BUGSoutput$summary

#Get some indicies used to parse data in su
sims = nrow(sl$BWT)
BWT.start = match("BWT[1]", row.names(su))
d18O.start = match("d18O_sw[1]", row.names(su))
MgCa.start = match("MgCa_sw_m[1]", row.names(su))
ts.len = ncol(sl$BWT)

#95% CI width
su.diff = as.double(su[,7] - su[,3])

#Average widths for BWT and d18O
mean(su.diff[BWT.start:(BWT.start+ts.len-1)])
mean(su.diff[d18O.start:(d18O.start+ts.len-1)])

##Multi-site analysis for sites 1123 and U1385

#Shorthand
sl = post.multi$BUGSoutput$sims.list
su = post.multi$BUGSoutput$summary

#Get some indicies used to parse data in su
sims = nrow(sl$BWT.b)
BWT.b.start = match("BWT.b[1]", row.names(su))
d18O.b.start = match("d18O_sw.b[1]", row.names(su))
BWT.e.start = match("BWT.e[1]", row.names(su))
d18O.e.start = match("d18O_sw.e[1]", row.names(su))
ts.len = ncol(sl$BWT.e)

#95% CI width
su.diff = as.double(su[,7] - su[,3])

#Average widths for each variable and record
mean(su.diff[BWT.b.start:(BWT.b.start+ts.len-1)])
mean(su.diff[d18O.b.start:(d18O.b.start+ts.len-1)])
mean(su.diff[BWT.e.start:(BWT.e.start+ts.len-1)])
mean(su.diff[d18O.e.start:(d18O.e.start+ts.len-1)])
