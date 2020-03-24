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

myGetObservables <- function(dim, vol, rho, beta) {
  
  y.iat <- array(NA, dim= c(length(beta), 2))
  y.q   <- array(NA, dim= c(length(beta), 2))
  y.E   <- array(NA, dim= c(length(beta), 2))
  y.lambda <- array(NA, dim= c(length(beta), 2))
  y.omega  <- array(NA, dim= c(length(beta), 2))
  y.kappa  <- array(NA, dim= c(length(beta), 2))
  
  for (i in c(1:length(beta))) {
    idx <- which(resFrame$dim==dim & resFrame$rho==rho & resFrame$vol==vol & resFrame$beta==beta[i])
    if (length(idx)!=0) {
      # y.iat[i,1] <- mean(resFrame$iat[idx])
      # y.iat[i,2] <- sqrt(var(resFrame$iat[idx]) / length(idx))
      # y.iat[i,2] <- getBootstrapMeanError(resFrame$iat[idx], N_replicas = 10)
      
      y.q[i,1] <- mean(resFrame$q[idx])
      # y.q[i,2] <- sqrt(var(resFrame$q[idx]) / length(idx))
      y.q[i,2] <- getBootstrapMeanError(resFrame$q[idx], N_replicas = 100)

      y.E[i,1] <- mean(resFrame$E_end[idx])
      # y.E[i,2] <- sqrt(var(resFrame$E_end[idx]) / length(idx))
      y.E[i,2] <- getBootstrapMeanError(resFrame$E_end[idx], N_replicas = 100)
    
      omegas <- numeric(0)
      for (j in c(1:length(idx))) {
        omegas <- c(omegas, resFrame$omega[[idx[j]]])
      }
      y.omega[i,1] <- mean(omegas)
      y.omega[i,2] <- sqrt(var(omegas) / length(omegas))
      # y.omega[i,2] <- getBootstrapMeanError(omegas, N_replicas = 100)
      rm(omegas)
      
      
      lambdas <- numeric(0)
      for (j in c(1:length(idx))) {
        lambdas <- c(lambdas, resFrame$lambda[[idx[j]]])
      }
      y.lambda[i,1] <- mean(lambdas)
      y.lambda[i,2] <- sqrt(var(lambdas) / length(lambdas))
      # y.lambda[i,2] <- getBootstrapMeanError(lambdas, N_replicas = 100)
      rm(lambdas)
      
      
      kappas <- numeric(0)
      for (j in c(1:length(idx))) {
        kappas <- c(kappas, resFrame$kappa[[idx[j]]])
      }
      y.kappa[i,1] <- mean(kappas)
      y.kappa[i,2] <- sqrt(var(kappas) / length(kappas))
      # y.kappa[i,2] <- getBootstrapMeanError(kappas, N_replicas = 100)
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

myPlotData <- function(x, y, yerr) {
  xlim <- c(1e-7, 1e7)
  # xlim <- c(1e-3, max(x.beta,na.rm=TRUE)+1)
  ylim <- c(0, max(y+yerr,na.rm=TRUE))
  # ylim <- c(0,10)
  
  # yerr[which(yerr < 0.1e0)] <- NA
  err_low <- y - yerr
  err_high<- y + yerr
  
  log <- "x"
  
  plot(x, y, xlim=xlim, ylim=ylim, type = "l", log=log)
  points(x, y, pch=4, col="black", cex=0.5)
  arrows(x0=x, y0=err_low, x1=x, y1=err_high, angle=90, code=3, length=0.02)
}

myMCMCSwipe <- function(dim, nPart=0, rho=0, vol=0, sigma=1, nIt=5e5, nResamples=1, x.beta, df) {
  cat("Starting swipe of ",nPart," particles (RUN_ID ",run_no,") over ",length(x.beta)," points:\n")
  
  if      (nPart == 0 & vol!=0 & rho!=0) nPart <- ceiling(vol * rho)
  else if (rho == 0 & vol!=0 & nPart!=0) rho <- nPart / vol
  else if (vol == 0 & rho!=0 & nPart!=0) vol <- nPart / rho
  
  l     <- vol**(1/dim)
  
  for (i in c(1:length(x.beta))) {
    cat("i:",i,", ",x.beta[i]," : ")
    
    watch.iter <- Sys.time()
    for (j in c(1:nResamples)) {
      cat(j, ", ")
      
      # Generate start configuration
      X0 <- array(data = runif(nPart*dim, min=-l/2, max=l/2), dim = c(nPart, dim))
      q  <- sample(c(-1,+1), size=nPart, replace=TRUE)
      
      # --------------- Positioned on a regular lattice
      X0 <- array(data = NA, dim = c(nPart, dim))
      nRow  <- ceiling(nPart**(1./dim))
      nCol  <- round  (nPart**(1./dim))

      rx <- vol**(1/2) / nRow / 2
      rx <- 1.5
      ry <- vol**(1/2) / nCol / 2
      ry <- 1.5
      offset <- vol**(1/2) * (1/2 - 1/2/nCol) - vol**(1/2) / 4

      for (iterator in c(1:nPart)) {
        X0[iterator,1] <- (floor((iterator-1) / nCol)) * rx - offset
        X0[iterator,2] <- ((iterator-1) %% nCol      ) * ry - offset
      }
      q <- rep(1, nPart)
      for (iterator in c(2:nPart)) {
        q[iterator] <- q[iterator-1]*-1
      }
      rm(rx,ry,offset)
      
      
      res <- MCMC_CPP(nIt, nPart, vol, 1/x.beta[i], sigma, X0, q) # Make Markov Chain
      
      
      # Calculate ovbservables
      # iat <- IAT_CPP(res$E[c(1e5:nIt)])
      iat <- 0

            # - - - - - - - - - SAVE RESULT TO FRAME
      df <- rbind(
        df, 
        list(
          "RUN_ID" = run_no,
          "IAT"    = iat, 
          "E_end"  = res$E[nIt], 
          "lambda" = I(list(res$lambda)),
          "omega"  = I(list(res$omega)),
          "q"      = res$q,
          "kappa"  = I(list(res$kappa)),
          "dim"    = dim,
          "rho"    = rho,
          "vol"    = vol,
          "beta"   = x.beta[i],
          "t"      = 1/x.beta[i],         
          "nIt"    = nIt,
          "nPart"  = nPart,
          "X_end"  = NA
          )
      )
      # - - - - - - - - - 
    }
    watch.iter <- as.double(Sys.time() - watch.iter, units='secs')
    cat("dur/stp:",watch.iter/nResamples,"\n")
  }
  
  return(df)
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

chisqr <- function(par, x, y, dy){
  f <- par[1]*x + par[2]
  return(sum((y-f)**2/dy**2))
}

myFit <- function(par, x, y, dy){
  myfit <- optim(par=par, fn=chisqr, y=y, dy=dy, x=x)
  return(list(a=myfit$par[1], b=myfit$par[2], chi2=myfit$value))
}

myFitErrors <- function(x, y, dy){
  N <- length(y)
  N_replicas <- 1500
  replicas <- array(data = NA, dim = c(N_replicas, N))
  for (i in c(1:N)) {
    replicas[ ,i] <- rnorm(n = N_replicas, mean=y[i], sd=dy[i])
  }
  #calculate best fit values (minimise chi2)
  par <- c(1,0)
  best_par <- array(data = NA, dim = c(N_replicas, length(par)+1))
  for (j in c(1:N_replicas)) {
    myfit <- optim(par = par, fn =chisqr, x=x, y=replicas[j, ], dy=dy)
    best_par[j, c(1:length(par))] <- myfit$par
    best_par[j, length(par)+1] <- myfit$value
  }  
  
  return(best_par)
}