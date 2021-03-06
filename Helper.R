# Frame Helper functions:
frameHelper.createNewFrame <- function(file_name) {
  # !!!!!!!!!!!!! BE CAREFUL NOT TO OVERWRITE EXISTING DATA !!!!!!!!!!!!!!!!!
  
  # Columns:
  # 1 "RUN_ID": Run No.
  # 2 "IAT"   : Integrated Autocorrelation Time
  # 3 "E_end" : Energy of Last State
  # 4 "lambda": Distance betw. particles
  # 5 "omega" : Nearest neigbour distance
  # 6 "q"     : Cluster sizes
  # 7 "kappa" : Number of particles in neighbourhood
  # 8 "dim"   : Number of Dimensions
  # 9 "rho"   : Particle Density
  #10 "vol"   : Cuboid Volume
  #11 "beta"  : Inverse Temperature
  #12 "t"     : Temperature
  #13 "nIt"   : Number of MCMC iteration steps
  #14 "nPart" : Number of Simulated Particles
  #15 X_end   : 
  
  resFrame <- data.frame(
    "RUN_ID"= 0,
    "IAT"   =NA,
    "E_end" =NA,
    "lambda"=NA,
    "omega" =NA,
    "q"     =NA,
    "kappa" =NA,
    "dim"   =NA,
    "rho"   =NA,
    "vol"   =NA,
    "beta"  =NA,
    "t"     =NA,
    "nIt"   =NA,
    "nPart" =NA,
    "X_end" =NA
  )
  save(resFrame, file=file_name)
  
  rm(resFrame)
}

frameHelper.showContentOverview <- function(file_name) {
  load(file_name)

  cat("This file contains the following parameter sets:\n")
  if (length(resFrame$RUN_ID)==1) {
    cat("File empty.\n")
    return()
  }
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

frameHelper.getBetaValues <- function(dim, vol, rho) {
  b <- resFrame$beta[which(resFrame$dim==dim & resFrame$vol==vol & resFrame$rho==rho)]
  b <- b[!duplicated(b)]
  b <- sort(b[!is.na(b)])
  return(b)
}

frameHelper.getObservables <- function(dim, vol, rho, beta) {
  
  y.iat <- array(NA, dim= c(length(beta), 2))
  y.q   <- array(NA, dim= c(length(beta), 2))
  y.E   <- array(NA, dim= c(length(beta), 2))
  y.lambda <- array(NA, dim= c(length(beta), 2))
  y.omega  <- array(NA, dim= c(length(beta), 2))
  y.kappa  <- array(NA, dim= c(length(beta), 2))
  
  for (i in c(1:length(beta))) {
    idx <- which(resFrame$dim==dim & resFrame$rho==rho & resFrame$vol==vol & resFrame$beta==beta[i])
    if (length(idx)!=0) {
  
      y.E[i,1] <- mean(resFrame$E_end[idx])
      y.E[i,2] <- getBootstrapMeanError(resFrame$E_end[idx], N_replicas = 100)
    
      qs <- numeric(0)
      for (j in c(1:length(idx))) {
        qs <- c(qs, mean(resFrame$q[[idx[j]]]))
      }
      y.q[i,1] <- mean(qs)
      y.q[i,2] <- getBootstrapMeanError(qs, N_replicas = 100)
      rm(qs)
      
      omegas <- numeric(0)
      for (j in c(1:length(idx))) {
        omegas <- c(omegas, mean(resFrame$omega[[idx[j]]]))
      }
      y.omega[i,1] <- mean(omegas)
      y.omega[i,2] <- getBootstrapMeanError(omegas, N_replicas = 100)
      rm(omegas)
      
      
      lambdas <- numeric(0)
      for (j in c(1:length(idx))) {
        lambdas <- c(lambdas, mean(resFrame$lambda[[idx[j]]]))
      }
      y.lambda[i,1] <- mean(lambdas)
      y.lambda[i,2] <- getBootstrapMeanError(lambdas, N_replicas = 100)
      rm(lambdas)
      
      
      kappas <- numeric(0)
      for (j in c(1:length(idx))) {
        kappas <- c(kappas, mean(resFrame$kappa[[idx[j]]]))
      }
      y.kappa[i,1] <- mean(kappas)
      y.kappa[i,2] <- getBootstrapMeanError(kappas, N_replicas = 100)
      rm(kappas)
    }
  }
  
  return(list(x=beta, 
              y.iat=y.iat[,1] , y.iat_Err=y.iat[,2],
              y.q  =y.q[,1]   , y.q_Err  =y.q[,2],
              y.E  =y.E[,1]   , y.E_Err  =y.E[,2],
              y.lambda=y.lambda[,1] , y.lambda_Err=y.lambda[,2],
              y.omega=y.omega[,1] , y.omega_Err=y.omega[,2],
              y.kappa=y.kappa[,1] , y.kappa_Err=y.kappa[,2]
              )
         )
}


# Drawing and plotting functions:
myPlotData <- function(x, y, dy) {
  xlim <- c(1e-7, 1e7)
  xlim <- c(1e-3, 1e0)
  ylim <- c(0, max(y+dy,na.rm=TRUE))
  plot(NA, NA, xlim=xlim, ylim=ylim, type = "l", log="x")
  draw.DataWithError(x=x, y=y, dy=dy, col="black", bars=TRUE, points=TRUE, lines=TRUE)
}

draw.DataWithError <- function(x, y, dy, col="black", bars=TRUE, points=TRUE, lines=TRUE) {
  units = par(c('usr', 'pin'))
  x_to_inches = with(units, pin[1L]/diff(usr[1:2]))
  y_to_inches = with(units, pin[2L]/diff(usr[3:4]))
  
  dists = sqrt((x_to_inches * dy)**2) 
  id_validArrow = which(dists > .0025)
  
  if (bars==TRUE)  points(x, y, pch=4, col=col, cex=0.5)
  if (lines==TRUE) lines(x, y, col=col, lwd=1)
  if (bars==TRUE)  arrows(x0=x[id_validArrow], y0=y[id_validArrow]-dy[id_validArrow], x1=x[id_validArrow], y1=y[id_validArrow]+dy[id_validArrow], angle=90, code=3, length=0.02, col=col)
}


# Error, Fit and Bootstrapping Functions
getBootstrapMeanError <- function(data, N_replicas=100){
  N <- length(data)
  # N_replicas <- replicaFactor*N
  id <- array(sample.int(n=N, size=N_replicas*N, replace = TRUE), dim = c(N_replicas, N))
  
  replicas <- array(data[id], dim = c(N_replicas, N))
  mean_rep <- apply(X = replicas, MARGIN = 1, FUN = mean)
  st_error <- sd(mean_rep)
  return(st_error)
}

chisqr_PiecewiseLinAndExpFunction <- function(par, x, y, dy) {
  # par[1]=A1, par[2]=X1, par[3]=A2, par[4]=X2
  f <- rep(0,length(x))
  if(par[2]<=0) return(Inf) # prevent the log of a negative number
  if(par[4]<=0) return(Inf) # just return chi2 of infinity
  if(par[4]<=par[2]) return(Inf)
  
  case1 <- which(x <= par[2])
  case2 <- which(x >  par[2] & x < par[4])
  case3 <- which(x >= par[4])
  
  f[case1] <- par[1]
  f[case2] <- (par[1] + (log(x[case2]) - log(par[2])) * (par[3]-par[1])/(log(par[4])-log(par[2])))
  f[case3] <- par[3]
  
  return(sum((y-f)**2/dy**2))
}

fit.ToData <- function(par, x, y, dy, fn){
  dy[which(dy==0)] <- 0.01
  dof <- length(y)-length(par)
  res <- optim(par=par, fn=fn, y=y, dy=dy, x=x)
  return(list(par=res$par, chi2=res$value, pval = 1-pchisq(q=res$value, df = dof), dof=dof))
}

fit.getErrors <- function(par, x, y, dy, fn){
  N <- length(y)
  dy[which(dy==0)] <- 0.01
  N_replicas <- 50
  replicas <- array(data = NA, dim = c(N_replicas, N))
  for (i in c(1:N)) {
    replicas[ ,i] <- rnorm(n = N_replicas, mean=y[i], sd=dy[i])
  }
  #calculate best fit values (minimise chi2)
  best_par <- array(data = NA, dim = c(N_replicas, length(par)+1))
  for (j in c(1:N_replicas)) {
    myfit <- optim(par = par, fn=fn, x=x, y=replicas[j, ], dy=dy)
    best_par[j, c(1:length(par))] <- myfit$par
    best_par[j, length(par)+1] <- myfit$value
  }  
  
  errors <- list(A1err=sd(best_par[ , 1]), X1err=sd(best_par[ ,2]), A2err=sd(best_par[ ,3]), X2err=sd(best_par[ ,4]), chi2err=sd(best_par[ ,5]))
  return(errors)
}

fit.automaticRoutine <- function(startparam, x, y, dy, xrange, fn) {
  id <- which(x > xrange[1] &  x < xrange[2])
  par    <- fit.ToData   (startparam, x[id], y[id], dy[id], fn)
  parerr <- fit.getErrors(par$par, x[id], y[id], dy[id], fn)
  return(list(par=par, parerr=parerr))
}

fit.draw <- function(par, x, col="black", vlines=TRUE, lines=TRUE) {
  if (vlines==TRUE) {
    abline(v=par$par[2], col=col, lwd=1)
    abline(v=par$par[4], col=col, lwd=1)
  }
  if (lines==TRUE) {
    case1 <- which(x <= par$par[2])
    case2 <- which(x >  par$par[2] & x < par$par[4])
    case3 <- which(x >= par$par[4])
    
    f <- numeric(length(x))
    f[case1] <- par$par[1]
    f[case2] <- (par$par[1] + (log(x[case2]) - log(par$par[2])) * (par$par[3]-par$par[1])/(log(par$par[4])-log(par$par[2])))
    f[case3] <- par$par[3]
    
    lines(x, f, col=col, lwd=1)
  }
}














