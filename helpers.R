#####
#Helper functions for JPI plotting etc
#####

##Make a density plot
plotd = function(x, col="black", lty=1, main="", xlab="", ylab="Density", 
                 xlim=c(min(x),max(x)), ylim=NULL, axes = TRUE) {
  if(is.null(ylim)){
    plot(density(x), col=col, lty=lty, main=main, xlab=xlab, ylab=ylab, xlim=xlim, axes = axes)  
  } else {
    plot(density(x), col=col, lty=lty, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes = axes)
  }
  
}

##Add line to density plot
lined = function(x, col="black", lty=1) {
  lines(density(x), col=col, lty=lty)
}

##Prep data for Lear analysis
prep.lear = function(){
  
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
  mgca_ts.min = 80
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
  d_mgca_calib = read.csv("O_mgca_calib.csv")
  
  ##Age index for MgCa calibration samples
  mgca_calib_age.ind = round((mgca_ts.min - d_mgca_calib$Age) / mgca_ts.step) + 1
  
  ##Read in d18O calibration dataset
  d_d18O_calib = read.csv("C_d18O_calib.csv")
  d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

  return(list(ts.ages = ts.ages, ts.len = ts.len, mgca_ts.ages = mgca_ts.ages,
              mgca_ts.len = mgca_ts.len, d_mgca = d_mgca, d_o = d_o, d_mgca_sw = d_mgca_sw,
              d_mgca_calib = d_mgca_calib, d_d18O_calib = d_d18O_calib,
              o_age.ind = o_age.ind, mgca_age.ind.all = mgca_age.ind.all,
              mgca_calib_age.ind = mgca_calib_age.ind, mgca_sw_age.ind = mgca_sw_age.ind))  
}

##Load data for site U1385 analysis
prep.birn = function(){
  
  ##Read proxy data
  d = read.csv("birner_2016.csv")
  
  ##Set up timeseries for d18O_sw and BWT modeling
  ts.min = 1320
  ts.max = 1235
  ts.step = 1
  ts.ages = seq(ts.min, ts.max, -ts.step)
  ts.len = length(ts.ages)
  
  #prep the d18O data and add age indicies
  d_o = d[!is.na(d$d18O), ]
  o_age.ind = round((ts.min - d_o$Age_ka) / ts.step) + 1
  
  #prep the MgCa data and add age indicies
  d_mgca = d[!is.na(d$MgCa),]
  mgca_age.ind = round((ts.min - d_mgca$Age_ka) / ts.step) + 1
  
  #U spp calibration data from Elderfield 2012 compilation
  d_mgca_calib = read.csv("U_mgca_calib.csv")
  
  ##Read in d18O calibration dataset
  d_d18O_calib = read.csv("C_d18O_calib.csv")
  d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]
  
  #get distributions of sw Mg/Ca from long model
  load("post_mg.RData") #uses output from MgCa_sw_model.R
  mgca_sw_m.neo = mean(post.mg$BUGSoutput$sims.list$MgCa_sw_m[,80])
  mgca_sw_sd.neo = sd(post.mg$BUGSoutput$sims.list$MgCa_sw_m[,80])
  mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)
  
  return(list(ts.ages = ts.ages, ts.len = ts.len, d_mgca = d_mgca, d_o = d_o, 
              d_mgca_calib = d_mgca_calib, d_d18O_calib = d_d18O_calib,
              o_age.ind = o_age.ind, mgca_age.ind = mgca_age.ind,
              mgca_sw_neo = mgca_sw_neo))  
}

##Prep data for site 1123 analysis
prep.elder = function(){
  #read data
  d = read.csv("elderfield_2012.csv")
  
  #subset desired age range
  d = d[d$Age_ka > 1235,]
  d = d[d$Age_ka < 1320,]
  
  ##Set up timeseries for d18O_sw and BWT modeling
  ts.min = 1320
  ts.max = 1235
  ts.step = 1
  ts.ages = seq(ts.min, ts.max, -ts.step)
  ts.len = length(ts.ages)
  
  #prep the d18O data and add age indicies
  d_o = d[!is.na(d$d18O), ]
  o_age.ind = round((ts.min - d_o$Age_ka) / ts.step) + 1
  
  #prep the MgCa data and add age indicies
  d_mgca = d[!is.na(d$MgCa),] 
  mgca_age.ind = round((ts.min - d_mgca$Age_ka) / ts.step) + 1
  
  #U spp calibration data from Elderfield 2012 compilation
  d_mgca_calib = read.csv("U_mgca_calib.csv")
  
  ##Read in d18O calibration dataset from Marchitto 2014 compilation
  d_d18O_calib = read.csv("U_d18O_calib.csv")
  d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]
  
  #Constraints based on LGM downcore, not used
#  LGM = c(-0.24, 1.7) #D_MgCa, D_d18O
  
  #get distributions of sw Mg/Ca from long model
  load("post_mg.RData") #uses output from MgCa_sw_model.R
  mgca_sw_m.neo = mean(post.mg$BUGSoutput$sims.list$MgCa_sw_m[,80])
  mgca_sw_sd.neo = sd(post.mg$BUGSoutput$sims.list$MgCa_sw_m[,80])
  mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)
  
  return(list(ts.ages = ts.ages, ts.len = ts.len, d_mgca = d_mgca, d_o = d_o, 
              d_mgca_calib = d_mgca_calib, d_d18O_calib = d_d18O_calib,
              o_age.ind = o_age.ind, mgca_age.ind = mgca_age.ind,
              mgca_sw_neo = mgca_sw_neo))  

}

##Prep data for multi-site analysis
prep.multi = function(){
  d.b = read.csv("birner_2016.csv")
  
  d.e = read.csv("elderfield_2012.csv")
  d.e = d.e[d.e$Age_ka > 1235,]
  d.e = d.e[d.e$Age_ka < 1320,]
  
  ##Set up timeseries for d18O_sw and BWT modeling
  ts.min = 1320
  ts.max = 1235
  ts.step = 1
  ts.ages = seq(ts.min, ts.max, -ts.step)
  ts.len = length(ts.ages)
  
  #prep the d18O data and add age indicies
  d_o.b = d.b[!is.na(d.b$d18O), ]
  o_age.ind.b = round((ts.min - d_o.b$Age_ka) / ts.step) + 1
  d_o.e = d.e[!is.na(d.e$d18O), ]
  o_age.ind.e = round((ts.min - d_o.e$Age_ka) / ts.step) + 1
  
  #prep the MgCa data and add age indicies
  d_mgca.b = d.b[!is.na(d.b$MgCa),]
  mgca_age.ind.b = round((ts.min - d_mgca.b$Age_ka) / ts.step) + 1
  d_mgca.e = d.e[!is.na(d.e$MgCa),] 
  mgca_age.ind.e = round((ts.min - d_mgca.e$Age_ka) / ts.step) + 1
  
  #U spp calibration data from Elderfield 2012 compilation
  d_mgca_calib = read.csv("U_mgca_calib.csv")
  
  ##Read in d18O calibration datasets for both genera
  d_d18O_calib.u = read.csv("U_d18O_calib.csv")
  d_d18O_calib.u = d_d18O_calib.u[is.na(d_d18O_calib.u$Ignore),]
  
  d_d18O_calib.c = read.csv("C_d18O_calib.csv")
  d_d18O_calib.c = d_d18O_calib.c[is.na(d_d18O_calib.c$Ignore),]
  
  #get distribution of sw Mg/Ca for early Pleistocene from long model
  load("post_mg.RData") #uses output from MgCa_sw_model.R
  mgca_sw_m.neo = mean(post.mg$BUGSoutput$sims.list$MgCa_sw_m[,80])
  mgca_sw_sd.neo = sd(post.mg$BUGSoutput$sims.list$MgCa_sw_m[,80])
  mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)
  
  return(list(ts.ages = ts.ages, ts.len = ts.len, d_mgca.b = d_mgca.b, d_o.b = d_o.b, 
              d_mgca.e = d_mgca.e, d_o.e = d_o.e, d_mgca_calib = d_mgca_calib, 
              d_d18O_calib.u = d_d18O_calib.u, d_d18O_calib.c = d_d18O_calib.c,
              o_age.ind.b = o_age.ind.b, mgca_age.ind.b = mgca_age.ind.b,
              o_age.ind.e = o_age.ind.e, mgca_age.ind.e = mgca_age.ind.e,
              mgca_sw_neo = mgca_sw_neo))  
  
}

