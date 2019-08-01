#####
#Helper functions for JPI plotting and data manipulation
#####

#####
#Make a density plot
#####

plotd = function(x, col="black", lty=1, main="", xlab="", ylab="Density", 
                 xlim=c(min(x),max(x)), ylim=NULL, axes = TRUE) {
  if(is.null(ylim)){
    plot(density(x), col=col, lty=lty, main=main, xlab=xlab, ylab=ylab, xlim=xlim, axes = axes)  
  } else {
    plot(density(x), col=col, lty=lty, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes = axes)
  }
  
}

#####
#Add line to density plot
#####

lined = function(x, col="black", lty=1) {
  lines(density(x), col=col, lty=lty)
}

#####
#Prep data for Site 806 analysis
#####

prep.lear = function(){
  
  ##Set up timeseries for d18O_sw and BWT modeling
  ts.min = 18
  ts.max = 0
  ts.step = 0.05
  ts.ages.base = seq(ts.min, ts.max, -ts.step)
  
  ##Prep the foram data, first read
  d = read.csv("Lear_combined.csv")
  
  ##Now split out the d18O data and strip one outlier
  d_o = d[!is.na(d$d18O),]
  d_o = d_o[d_o$Sample.ID != "806B 47-5 38-43",]

  #Add O data ages to age vector
  ts.ages = c(ts.ages.base, d_o$Age.Ma)
  
  ##Now split out the MgCa data
  d_mgca = d[!is.na(d$MgCa),]
  
  #Add Mg/Ca data ages to age vector, condense, and sort
  ts.ages = c(ts.ages, d_mgca$Age.Ma)
  ts.ages = unique(ts.ages)
  ts.ages = sort(ts.ages, decreasing = TRUE)
  
  #Timeseries index for each proxy series and length
  ts.len = length(ts.ages)
  o_age.ind = match(d_o$Age.Ma, ts.ages)
  mgca_age.ind = match(d_mgca$Age.Ma, ts.ages)
  
  ##Set up timeseries for MgCa_sw modeling
  mgca_ts.min = 80
  mgca_ts.max = 0
  mgca_ts.step = 1
  mgca_ts.ages = seq(mgca_ts.min, mgca_ts.max, -mgca_ts.step)
  
  ##Read in paleo-seawater MgCa data
  d_mgca_sw = read.csv("mgca_sw.csv")
  
  ##Read in MgCa calibration dataset
  d_mgca_calib = read.csv("O_mgca_calib.csv")
  
  ##Append all sample ages to ts age vector
  mgca_ages = c(mgca_ts.ages, d_mgca_sw$Age, d_mgca$Age.Ma, d_mgca_calib$Age)
  mgca_ages = unique(mgca_ages)
  mgca_ages = sort(mgca_ages, decreasing = TRUE)
  
  ##Age index for seawater MgCa proxy data and ts length
  mgca_ages.len = length(mgca_ages)
  mgca_sw_age.ind = match(d_mgca_sw$Age, mgca_ages)
  
  ##Add age indicies for seawater MgCa TS to MgCa foram data
  mgca_age.ind.sw = match(d_mgca$Age.Ma, mgca_ages)
  mgca_age.ind.all = matrix(c(mgca_age.ind, mgca_age.ind.sw), ncol = 2)
  
  ##Age index for MgCa calibration samples
  mgca_calib_age.ind = match(d_mgca_calib$Age, mgca_ages)
  
  ##Read in d18O calibration dataset
  d_d18O_calib = read.csv("C_d18O_calib.csv")
  d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

  return(list(ts.ages.base = ts.ages.base, ts.ages = ts.ages, ts.len = ts.len, mgca.ages = mgca_ages,
              mgca_ts.len = mgca_ages.len, d_mgca = d_mgca, d_o = d_o, d_mgca_sw = d_mgca_sw,
              d_mgca_calib = d_mgca_calib, d_d18O_calib = d_d18O_calib,
              o_age.ind = o_age.ind, mgca_age.ind.all = mgca_age.ind.all,
              mgca_calib_age.ind = mgca_calib_age.ind, mgca_sw_age.ind = mgca_sw_age.ind))  
}

#####
#Prep data for site U1385 analysis
#####

prep.birn = function(){
  
  ##Read proxy data
  d = read.csv("birner_2016.csv")
  
  ##Set up timeseries for d18O_sw and BWT modeling
  ts.min = 1320
  ts.max = 1235
  ts.step = 1
  ts.ages.base = seq(ts.min, ts.max, -ts.step)

  #prep the d18O and MgCa data
  d_o = d[!is.na(d$d18O), ]
  d_mgca = d[!is.na(d$MgCa),]
  
  #append proxy record ages and condense
  ts.ages = c(ts.ages.base, d_o$Age_ka, d_mgca$Age_ka)
  ts.ages = unique(ts.ages)
  ts.ages = sort(ts.ages, decreasing = TRUE)
  
  #Timeseries index for each proxy series and length
  ts.len = length(ts.ages)
  o_age.ind = match(d_o$Age_ka, ts.ages)
  mgca_age.ind = match(d_mgca$Age_ka, ts.ages)
  
  #U spp calibration data from Elderfield 2012 compilation
  d_mgca_calib = read.csv("U_mgca_calib.csv")
  
  ##Read in d18O calibration dataset
  d_d18O_calib = read.csv("C_d18O_calib.csv")
  d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]
  
  ##Get distributions of sw Mg/Ca from long model
  load("post_lear.RData") #uses output from site 806 JPI analysis
  mgca_sw_m.neo = mean(post.lear$BUGSoutput$sims.list$MgCa_sw_m[,80]) #column 80
  mgca_sw_sd.neo = sd(post.lear$BUGSoutput$sims.list$MgCa_sw_m[,80])  #is 1 Ma
  mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)
  
  return(list(ts.ages.base = ts.ages.base, ts.ages = ts.ages, ts.len = ts.len, d_mgca = d_mgca, d_o = d_o, 
              d_mgca_calib = d_mgca_calib, d_d18O_calib = d_d18O_calib,
              o_age.ind = o_age.ind, mgca_age.ind = mgca_age.ind,
              mgca_sw_neo = mgca_sw_neo))  
}

#####
#Prep data for site 1123 analysis
#####

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
  ts.ages.base = seq(ts.min, ts.max, -ts.step)

  #prep the proxy data
  d_o = d[!is.na(d$d18O), ]
  d_mgca = d[!is.na(d$MgCa),] 

  #append proxy record ages and condense
  ts.ages = c(ts.ages.base, d_o$Age_ka, d_mgca$Age_ka)
  ts.ages = unique(ts.ages)
  ts.ages = sort(ts.ages, decreasing = TRUE)
  
  #Timeseries index for each proxy series and length
  ts.len = length(ts.ages)
  o_age.ind = match(d_o$Age_ka, ts.ages)
  mgca_age.ind = match(d_mgca$Age_ka, ts.ages)

  #U spp calibration data from Elderfield 2012 compilation
  d_mgca_calib = read.csv("U_mgca_calib.csv")
  
  ##Read in d18O calibration dataset from Marchitto 2014 compilation
  d_d18O_calib = read.csv("U_d18O_calib.csv")
  d_d18O_calib = d_d18O_calib[is.na(d_d18O_calib$Ignore),]

  ##Get distributions of sw Mg/Ca from long model
  load("post_lear.RData") #uses output from site 806 JPI analysis
  mgca_sw_m.neo = mean(post.lear$BUGSoutput$sims.list$MgCa_sw_m[,80]) #column 80
  mgca_sw_sd.neo = sd(post.lear$BUGSoutput$sims.list$MgCa_sw_m[,80])  #is 1 Ma
  mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)
  
  return(list(ts.ages.base = ts.ages.base, ts.ages = ts.ages, ts.len = ts.len, d_mgca = d_mgca, d_o = d_o, 
              d_mgca_calib = d_mgca_calib, d_d18O_calib = d_d18O_calib,
              o_age.ind = o_age.ind, mgca_age.ind = mgca_age.ind,
              mgca_sw_neo = mgca_sw_neo))  

}

#####
#Prep data for multi-site analysis
#####

prep.multi = function(){
  
  ##Read site U1385 data
  d.b = read.csv("birner_2016.csv")
  
  ##Read site 1123 data and subset for target time interval
  d.e = read.csv("elderfield_2012.csv")
  d.e = d.e[d.e$Age_ka > 1235,]
  d.e = d.e[d.e$Age_ka < 1320,]
  
  ##Set up timeseries for d18O_sw and BWT modeling
  ts.min = 1320
  ts.max = 1235
  ts.step = 1
  ts.ages.base = seq(ts.min, ts.max, -ts.step)
  
  #prep the d18O data
  d_o.b = d.b[!is.na(d.b$d18O), ]
  d_o.e = d.e[!is.na(d.e$d18O), ]
  
  #prep the MgCa data
  d_mgca.b = d.b[!is.na(d.b$MgCa),]
  d_mgca.e = d.e[!is.na(d.e$MgCa),] 
  
  #append proxy record ages and condense
  ts.ages = c(ts.ages.base, d_o.b$Age_ka, d_o.e$Age_ka, d_mgca.b$Age_ka, d_mgca.e$Age_ka)
  ts.ages = unique(ts.ages)
  ts.ages = sort(ts.ages, decreasing = TRUE)
  
  #Timeseries index for each proxy series and length
  ts.len = length(ts.ages)
  o_age.ind.b = match(d_o.b$Age_ka, ts.ages)
  o_age.ind.e = match(d_o.e$Age_ka, ts.ages)
  mgca_age.ind.b = match(d_mgca.b$Age_ka, ts.ages)
  mgca_age.ind.e = match(d_mgca.e$Age_ka, ts.ages)
  
  #U spp calibration data from Elderfield 2012 compilation
  d_mgca_calib = read.csv("U_mgca_calib.csv")
  
  ##Read in d18O calibration datasets for both genera
  d_d18O_calib.u = read.csv("U_d18O_calib.csv")
  d_d18O_calib.u = d_d18O_calib.u[is.na(d_d18O_calib.u$Ignore),]
  
  d_d18O_calib.c = read.csv("C_d18O_calib.csv")
  d_d18O_calib.c = d_d18O_calib.c[is.na(d_d18O_calib.c$Ignore),]
  
  ##Get distributions of sw Mg/Ca from long model
  load("post_lear.RData") #uses output from site 806 JPI analysis
  mgca_sw_m.neo = mean(post.lear$BUGSoutput$sims.list$MgCa_sw_m[,80]) #column 80
  mgca_sw_sd.neo = sd(post.lear$BUGSoutput$sims.list$MgCa_sw_m[,80])  #is 1 Ma
  mgca_sw_neo = c(mgca_sw_m.neo, mgca_sw_sd.neo)
  
  return(list(ts.ages.base = ts.ages.base, ts.ages = ts.ages, ts.len = ts.len, d_mgca.b = d_mgca.b, d_o.b = d_o.b, 
              d_mgca.e = d_mgca.e, d_o.e = d_o.e, d_mgca_calib = d_mgca_calib, 
              d_d18O_calib.u = d_d18O_calib.u, d_d18O_calib.c = d_d18O_calib.c,
              o_age.ind.b = o_age.ind.b, mgca_age.ind.b = mgca_age.ind.b,
              o_age.ind.e = o_age.ind.e, mgca_age.ind.e = mgca_age.ind.e,
              mgca_sw_neo = mgca_sw_neo))  
  
}

