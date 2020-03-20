myCreateNewFrame <- function(file_name) {
  # !!!!!!!!!!!!! BE CAREFUL NOT TO OVERWRITE EXISTING DATA !!!!!!!!!!!!!!!!!
  
  # Columns:
  # 1 "RUN_ID": Run No.
  # 2 "IAT"  : Integrated Autocorrelation Time
  # 3 "E_end": Energy of Last State
  # 4 "dim"  : Number of Dimensions
  # 5 "rho"  : Particle Density
  # 6 "vol"  : Cuboid Volume
  # 7 "beta" : Inverse Temperature
  # 8 "nIt"  : Number of MCMC iteration steps
  # 9 "nPart": Number of Simulated Particles
  #10 "t"    : Temperature
  #11 X_end  : 
  
  resFrame <- data.frame(
    "RUN_ID"= 0,
    "IAT"   =NA,
    "E_end" =NA,
    "q"     =NA,
    "dim"   =NA,
    "rho"   =NA,
    "vol"   =NA,
    "beta"  =NA,
    "nIt"   =NA,
    "nPart" =NA,
    "t"     =NA
  )
  save(resFrame, file=file_name)
  
  rm(resFrame)
  rm(file_name)
}

myShowContentOverview <- function(file_name) {
  load(file_name)

  cat("This file contains the following parameter sets:\n")
  cat("\t d \t vol \t rho \t beta \t N_beta \t N_samples\n")
  v.dims <- resFrame$dim
  v.dims <- v.dims[!duplicated(v.dims)]
  v.dims <- sort(v.dims[!is.na(v.dims)])
  
  for (c.d in c(1:length(v.dims))) {
    v.vol <- resFrame$vol[which(resFrame$dim==v.dims[c.d])]
    v.vol <- v.vol[!duplicated(v.vol)]
    v.vol <- sort(v.vol[!is.na(v.vol)])
    
    for (c.v in c(1:length(v.vol))) {
      v.rho <- resFrame$rho[which(resFrame$dim==v.dims[c.d] & resFrame$vol==v.vol[c.v])]
      v.rho <- v.rho[!duplicated(v.rho)]
      v.rho <- sort(v.rho[!is.na(v.rho)])
      
      for (c.r in c(1:length(v.rho))) {
        v.beta <- resFrame$beta[which(resFrame$dim==v.dims[c.d] & resFrame$vol==v.vol[c.v] & resFrame$rho==v.rho[c.r])]
        v.beta <- v.beta[!duplicated(v.beta)]
        v.beta <- sort(v.beta[!is.na(v.beta)])
        
        N_beta <- length(v.beta)
        N_samples <- length(which(resFrame$dim==v.dims[c.d] & resFrame$vol==v.vol[c.v] & resFrame$rho==v.rho[c.r]))
        cat("\t",v.dims[c.d],"\t",v.vol[c.v],"\t",v.rho[c.r],"\t \t", N_beta, "\t\t", N_samples ,"\n")
      }
    }  
  }
  
  rm(v.dims, v.vol, v.rho, v.beta, c.d, c.v, c.r, N_beta, N_samples)
}

myGetBetaValues <- function(dim, vol, rho) {
  b <- resFrame$beta[which(resFrame$dim==dim & resFrame$vol==vol & resFrame$rho==rho)]
  b <- b[!duplicated(b)]
  b <- sort(b[!is.na(b)])
  
  # rm(dim, vol, rho)
  return(b)
}

myGetMeanQValues <- function(dim, vol, rho, beta) {
  y.q <- array(NA, dim= c(length(beta), 2)) # [,1]: IAT, [,2]: Standard Error
  # y.qsd <- rep(NA, length(beta))
  
  for (i in c(1:length(beta))) {
    idx <- which(resFrame$dim==d & resFrame$rho==particleDensity & resFrame$vol==vol & resFrame$beta==beta[i])
    
    if (length(idx)!=0) {
      y.q[i,1] <- mean(resFrame$q[idx])
      # y.q[i,2] <- sqrt(var(resFrame$q[idx]) / length(idx))
      y.q[i,2] <- getBootstrapMeanError(resFrame$q[idx])
      # y.qsd[i] <- sd(resFrame$q[idx])
    }
  }
  
  return(list(x=beta, y=y.q[,1], yerr=y.q[,2]))
}

myPlotData <- function(plt.data) {
  x <- plt.data$x
  y <- plt.data$y
  yerr <- plt.data$yerr
  
  xlim <- c(5e-5, 1e5)
  # xlim <- c(1e-3, max(x.beta,na.rm=TRUE)+1)
  ylim <- c(1, max(y+yerr,na.rm=TRUE))
  # ylim <- c(0,10)
  
  yerr[which(yerr < 1e0)] <- NA
  err_low <- y - yerr
  err_high<- y + yerr
  
  log <- "x"
  
  plot(x, y, xlim=xlim, ylim=ylim, type = "l", log=log)
  points(x, y, pch=4, col="black", cex=0.5)
  arrows(x0=x, y0=err_low, x1=x, y1=err_high, angle=90, code=3, length=0.02)
}