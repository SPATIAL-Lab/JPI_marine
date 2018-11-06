#####
#Old version of the Lear analysis using discrete measurements
#####

library(R2OpenBUGS)
library(coda)
library(rjags)
library(R2jags)
library(xlsx)

setwd("C:/Users/gjbowen/Dropbox/HypoMirror/JPI_marine/")
setwd("C:/Users/u0133977/Dropbox/HypoMirror/JPI_marine/")
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