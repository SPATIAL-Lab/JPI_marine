#####
#Helper functions for JPI plotting etc
#####

##Make a density plot
plotd = function(x, col="black", lty=1, main="", xlab="", ylab="Density", xlim=c(min(x),max(x)), ylim=NULL) {
  if(is.null(ylim)){
    plot(density(x), col=col, lty=lty, main=main, xlab=xlab, ylab=ylab, xlim=xlim)  
  } else {
    plot(density(x), col=col, lty=lty, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  }
  
}

##Add line to density plot
lined = function(x, col="black", lty=1) {
  lines(density(x), col=col, lty=lty)
}
