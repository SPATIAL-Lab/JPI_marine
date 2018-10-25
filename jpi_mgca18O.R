#####

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)
library(xlsx)

setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/")
source("lear_model.R")

##Analysis of Lear 2015 data
d = read.xlsx("Lear_2015_data.xlsx", sheetIndex = 1)
d = d[!is.na(d$d18O),]

parameters = c("d18O.sw", "BWT", "MgCa.sw", "le")
set.seed(1395)

posts = list()

for(i in 1:nrow(d)){
  dat = list(MgCa = d$MgCa[i], d18O = d$d18O[i], Age = d$Age.Ma[i])
  post = jags(model.file = textConnection(lear), parameters.to.save = parameters, 
               data = dat, inits = NULL, n.chains=3, n.iter = 50000, 
               n.burnin = 1000, n.thin = 25)  
  posts[[i]] = post
}

ptiles = list()
for(i in 1:nrow(d)){
  post.df = as.data.frame(posts[[i]]$BUGSoutput$sims.list)
  
  #gather up some stats on the posterior
  ptiles[[i]] =  sapply(post.df, quantile, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))
}

ptiles[[1]]  

nr = nrow(ptiles[[1]])
nc = ncol(ptiles[[1]])

ptiles.df = data.frame(matrix(ncol = nr*nc, nrow = 0))
for(i in 1:nc){
  for(j in 1:nr){
    colnames(ptiles.df)[j+(i-1)*nr] = paste0(names(ptiles[[1]][1,][i]), "_", names(ptiles[[1]][,1][j]))
    for(k in 1:length(ptiles)){
      ptiles.df[k, j+(i-1)*nr] = ptiles[[k]][j,i]
    }
  }
}
ptiles.df$Age = d$Age

jpeg("learout.jpeg", res=300, units="in", width=8, height=8)
layout(matrix(c(1,2), 2, 1))
par(mar = c(5,5,1,1))
plot(ptiles.df$Age, ptiles.df$`BWT_50%`, type="l", ylim=c(4,11), xlab = "Age (Ma)", ylab="Bottom water temp")
polygon(c(ptiles.df$Age, rev(ptiles.df$Age)), c(ptiles.df$`BWT_2.5%`, rev(ptiles.df$`BWT_97.5%`)), col="grey70", lty=0)
lines(ptiles.df$Age, ptiles.df$`BWT_50%`)
lines(ptiles.df$Age, ptiles.df$`BWT_2.5%`, lty=2)
lines(ptiles.df$Age, ptiles.df$`BWT_97.5%`, lty=2)

par(mar = c(5,5,1,1))
plot(ptiles.df$Age, ptiles.df$`d18O.sw_50%`, type="l", ylim=c(-1,1.5), xlab = "Age (Ma)", 
     ylab=expression(paste("Seawater ",delta^{18}, "O (\u2030)")))
polygon(c(ptiles.df$Age, rev(ptiles.df$Age)), c(ptiles.df$`d18O.sw_2.5%`, rev(ptiles.df$`d18O.sw_97.5%`)), col="grey70", lty=0)
lines(ptiles.df$Age, ptiles.df$`d18O.sw_50%`)
lines(ptiles.df$Age, ptiles.df$`d18O.sw_2.5%`, lty=2)
lines(ptiles.df$Age, ptiles.df$`d18O.sw_97.5%`, lty=2)

dev.off()

#####V1 of autoregressive model

d = read.xlsx("Lear_2015_data.xlsx", sheetIndex = 1)
d = d[!is.na(d$d18O),]

parameters = c("d18O_sw", "BWT")
set.seed(1395)

ts.min = 18
ts.max = 11
ts.step = 0.1
ts.ages = seq(18, 11, -0.1)
ts.len = length(ts.ages)

age.ind = round((ts.min - d$Age.Ma) / ts.step) + 1

dat = list(nages = ts.len, nobs = nrow(d), MgCa = d$MgCa, d18O = d$d18O, Age = d$Age.Ma, age.ind = age.ind)

source("lear_temporal.R")

post = jags(model.file = textConnection(lear_AR), parameters.to.save = parameters, 
              data = dat, inits = NULL, n.chains=3, n.iter = 50000, 
              n.burnin = 1000, n.thin = 25)  

plot(0, 0, xlab="Age", ylab ="Temperature", xlim=c(11,18), ylim=c(4,11))
for(i in 1:nrow(post$BUGSoutput$sims.list$BWT)){
  lines(ts.ages, post$BUGSoutput$sims.list$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, post$BUGSoutput$summary[1:71, 5], col="red")
lines(ts.ages, post$BUGSoutput$summary[1:71, 3], col="red", lty=3)
lines(ts.ages, post$BUGSoutput$summary[1:71, 7], col="red", lty=3)
points(d$Age.Ma, rep(4, nrow(d)), pch=21, bg = "white")

plot(0, 0, xlab="Age", ylab ="Seawater d18O", xlim=c(11,18), ylim=c(-1,1.5))
for(i in 1:nrow(post$BUGSoutput$sims.list$d18O_sw)){
  lines(ts.ages, post$BUGSoutput$sims.list$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, post$BUGSoutput$summary[72:142, 5], col="red")
lines(ts.ages, post$BUGSoutput$summary[72:142, 3], col="red", lty=3)
lines(ts.ages, post$BUGSoutput$summary[72:142, 7], col="red", lty=3)
points(d$Age.Ma, rep(-1, nrow(d)), pch=21, bg = "white")

###version 2 splits MgCa and d18O data, adds swMgCa model 

set.seed(19395)

##Set up timeseries for d18O_sw and BWT modeling
ts.min = 18
ts.max = 11
ts.step = 0.1
ts.ages = seq(ts.min, ts.max, -ts.step)
ts.len = length(ts.ages)

##Prep the foram data, first read
d = read.xlsx("Lear_2015_data.xlsx", sheetIndex = 1)

##Now split out the d18O data and strip one outlier
d_o = d[!is.na(d$d18O),]
d_o = d_o[d_o$Depth.m != 442.48,]
#Timeseries index for each d18O sample
o_age.ind = round((ts.min - d_o$Age.Ma) / ts.step) + 1

##Now split out the MgCa data and get TS index
d_mgca = d[!is.na(d$MgCa),]
mgca_age.ind = round((ts.min - d_mgca$Age.Ma) / ts.step) + 1

##Read in paleo-seawater MgCa data
d_mgca_sw = read.csv("mgca_sw.txt")

##Read in MgCa calibration dataset
d_mgca_calib = read.csv("mgca_calib.csv")

##Parameters to be saved
parameters = c("d18O_sw", "BWT", "BWT.eps.ac", "BWT.var", "MgCa_sw.b", "lc",
               "d18O_sw.eps.ac", "d18O_sw.var")

##Data to pass to BUGS model
dat = list(nages = ts.len, age.old = ts.min, age.int = ts.step,
           MgCa_calib.age = d_mgca_calib$Age, MgCa_calib.bwt = d_mgca_calib$BWT, MgCa_calib = d_mgca_calib$MgCa,
           MgCa_sw.age = d_mgca_sw$Age, MgCa_sw = d_mgca_sw$MgCa, MgCa_sw.sd = d_mgca_sw$Sigma,
           MgCa.age.ind = mgca_age.ind, MgCa = d_mgca$MgCa, 
           d18O.age.ind = o_age.ind, d18O = d_o$d18O)

##Here's the BUGS code
source("split_temporal.R")

##Run the inversion
t1 = proc.time()
post2 = jags(model.file = textConnection(split_AR), parameters.to.save = parameters, 
            data = dat, inits = NULL, n.chains=3, n.iter = 50000, 
            n.burnin = 1000, n.thin = 25)  
proc.time() - t1

##A couple of standard plots of the modeled timeseries
plot(0, 0, xlab="Age", ylab ="Temperature", xlim=c(11,18), ylim=c(4,11))
for(i in 1:nrow(post2$BUGSoutput$sims.list$BWT)){
  lines(ts.ages, post2$BUGSoutput$sims.list$BWT[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, post2$BUGSoutput$summary[1:ts.len, 5], col="red")
lines(ts.ages, post2$BUGSoutput$summary[1:ts.len, 3], col="red", lty=3)
lines(ts.ages, post2$BUGSoutput$summary[1:ts.len, 7], col="red", lty=3)
points(d_mgca$Age.Ma, rep(4, nrow(d_mgca)), pch=21, bg = "white")

plot(0, 0, xlab="Age", ylab ="Seawater d18O", xlim=c(11,18), ylim=c(-1,1.5))
for(i in 1:nrow(post2$BUGSoutput$sims.list$d18O_sw)){
  lines(ts.ages, post2$BUGSoutput$sims.list$d18O_sw[i,], col = rgb(0,0,0, 0.01))
}
lines(ts.ages, post2$BUGSoutput$summary[(ts.len+7):(ts.len*2+6), 5], col="red")
lines(ts.ages, post2$BUGSoutput$summary[(ts.len+7):(ts.len*2+6), 3], col="red", lty=3)
lines(ts.ages, post2$BUGSoutput$summary[(ts.len+7):(ts.len*2+6), 7], col="red", lty=3)
points(d_o$Age.Ma, rep(-1, nrow(d_o)), pch=21, bg = "white")


###Now let's try to Shackelton site
d.i = read.table("Birner_2016/datasets/339-U1385_isotope_toRead.tab", sep = "\t", header = TRUE)
d.e = read.table("Birner_2016/datasets/339-U1385_Mg-Ca_toRead.tab", sep = "\t", header = TRUE)
plot(d.e$Age..ka.BP., d.e$BWT..�.C.)