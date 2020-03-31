library(Rcpp)
library(rgl)           # 3D Plotting
library(tidyverse)     # Remove duplicate elements
library(ggplot2)
library(ggforce)

setwd(dirname(parent.frame(2)$ofile))

show_particles <- function(x, q, V){
  dim <- length(x[1,])
  l <- V**(1/dim)
  
  ipos <- which(q==+1)
  ineg <- which(q==-1)
  posParticles <- x[ipos,]
  negParticles <- x[ineg,]
  
  colours <- numeric(length(q))
  colours[which(q==+1)] <- "red"
  colours[which(q==-1)] <- "blue"
  
  if (dim == 1) {
    yLength <- 0.5
    plot(NA, xlim=c(-l/2, l/2), ylim=c(-1.5,1.5), xlab='x', ylab='', yaxt='n')
    
    if (length(ipos)!=0) {
      arrows(x0=x[ipos,]-0.5,
             y0=0,
             x1=x[ipos,]+0.5,
             y1=0,
             length=0.02, angle=90, code=3, col='red')
      arrows(x0=x[ipos,],
             y0=yLength,
             x1=x[ipos,],
             y1=-yLength,
             length=0.02, angle=90, code=0, col='red')
    }
    if (length(ineg)!=0) {
      arrows(x0=x[ineg,]-0.5,
             y0=0,
             x1=x[ineg,]+0.5,
             y1=0,
             length=0.02, angle=90, code=3, col='blue')
      arrows(x0=x[ineg,],
             y0=yLength,
             x1=x[ineg,],
             y1=-yLength,
             length=0.02, angle=90, code=0, col='blue')
    }    
  } else if (dim == 2) {
    df_pos <- data.frame("x"=x[ipos,1],"y"=x[ipos,2])
    df_neg <- data.frame("x"=x[ineg,1],"y"=x[ineg,2])
    
    gg <- ggplot() 
    gg <- gg + geom_circle(
      data = df_pos,
      mapping = aes(x0 = x,  y0 = y, r = 0.5),
      fill = "red")
    gg <- gg + geom_circle(
      data = df_neg,
      mapping = aes(x0 = x,  y0 = y, r = 0.5),
      fill = "blue")
    
    
    gg <- gg + xlim(c(-l/2, +l/2)*1.05) + ylim(c(-l/2, +l/2)*1.05)
    gg <- gg + coord_fixed()
    gg <- gg + theme_bw()
    
    gg <- gg + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), size=0.5, data=data.frame(x1=-l/2,x2=+l/2,y1=-l/2,y2=-l/2))
    gg <- gg + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), size=0.5, data=data.frame(x1=-l/2,x2=+l/2,y1=+l/2,y2=+l/2))
    gg <- gg + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), size=0.5, data=data.frame(x1=-l/2,x2=-l/2,y1=-l/2,y2=+l/2))
    gg <- gg + geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), size=0.5, data=data.frame(x1=+l/2,x2=+l/2,y1=-l/2,y2=+l/2))
    
    gg <- gg + labs(x = "x-axis", y = "y-axis")
    
    plot(gg)
    
  } else if (dim == 3) {
    ballradius <- 0.5
    
    open3d()
    plot3d(x=NA, y=NA, z=NA, xlim=c(-l/2,l/2),ylim=c(-l/2,l/2),zlim=c(-l/2,l/2), xlab="x-axis", ylab="y-axis", zlab="z-axis")
    spheres3d(x=posParticles[,1], y=posParticles[,2], z=posParticles[,3], col = "red",
              radius=ballradius)
    spheres3d(x=negParticles[,1], y=negParticles[,2], z=negParticles[,3], col = "blue",
              radius=ballradius)
    
    rgl.viewpoint( theta = 240, phi = 10)
    
  }
}

gen_initConfig <- function(type="random", dim, vol, nPart, lc=1.5) {
  if (type!="random" && type!="ordered") {
    cat ("Invalid start configuration type: ",type,"! Setting to random.\n")
    type <- "random"
  }
  
  X0 <- numeric(0)
  q  <- rep(1,nPart)
  
  l <- vol**(1./dim)
  
  if (type == "random") 
  {
    X0 <- array(data = runif(nPart*dim, min=-l/2, max=l/2), dim = c(nPart, dim))
    q  <- c(rep(+1,ceiling(nPart/2)),rep(-1,floor(nPart/2)))
    # q  <- sample(c(-1,+1), size=nPart, replace=TRUE)
  } 
  else if (type == "ordered") 
  {
    X0 <- array(data = NA, dim = c(nPart, dim))
    nRow  <- ceiling(nPart**(1./dim))
    nCol  <- round  (nPart**(1./dim))
    
    offset <- l / 2 * (1 - 1/nCol) - l / 4
    
    if (dim==1) {
      for (i in c(1:nPart)) {
        X0[i,1] <- (i-1) * lc - offset
      }
    } else if (dim == 2) {
      for (i in c(1:nPart)) {
        X0[i,1] <- (floor((i-1) / nCol)) * lc - offset
        X0[i,2] <- ((i-1) %% nCol      ) * lc - offset
      }
    } else if (dim == 3) {
      for (i in c(1:nPart)) {
        a <- floor( (i-1)                         / nRow/nCol)
        b <- floor(((i-1)         - a*nCol*nRow)  / nCol     )
        c <-       ((i-1)- b*nCol - a*nCol*nRow) %% nCol
        
        X0[i,1] <- a * lc - offset
        X0[i,2] <- b * lc - offset
        X0[i,3] <- c * lc - offset
        
        q[i] <- (-1)**(a+b+c)
      }
    } else if (dim == 4) {
      if (nPart!=81) {
        cat("ERROR: Startconfiguration in 4D only available for 81 particles!")
        return()
      }
      for (a in c(0:2)) {
        for (b in c(0:2)) {
          for (c in c(0:2)) {
            for (d in c(0:2)) {
              i <- a * 27 + b * 9 + c * 3 + d
              
              X0[i,1] <- a * lc - offset
              X0[i,2] <- b * lc - offset
              X0[i,3] <- c * lc - offset
              X0[i,4] <- d * lc - offset
              
              q[i] <- (-1)**(a+b+c+d)
            }
          }
        }
      }
    }
    
    
    if (dim==1 | dim==2) {
      for (iterator in c(2:nPart)) {
        q[iterator] <- q[iterator-1]*-1
      }
    }
  }
  
  return(list(X0=X0, q=q))
}

myMCMCSwipe <- function(dim, nPart=0, rho=0, vol=0, sigma=1, sct="random", nIt=5e5, nResamples=1, x.beta, df) {
  cat("Starting swipe of ",nPart," particles (RUN_ID ",run_no,") over ",length(x.beta)," points:\n")
  
  if      (nPart == 0 & vol!=0 & rho!=0) nPart <- ceiling(vol * rho)
  else if (rho == 0 & vol!=0 & nPart!=0) rho <- nPart / vol
  else if (vol == 0 & rho!=0 & nPart!=0) vol <- nPart / rho
  
  l     <- vol**(1/dim)
  watch.records <- rep(NA, length(x.beta)) 
  
  for (i in c(1:length(x.beta))) {
    cat("i:",i,"\t, ",log10(x.beta[i]),"   \t : ")
    
    watch.iter <- Sys.time()
    for (j in c(1:nResamples)) {
      cat(j, ", ")
      
      # Generate start configuration
      init_config <- gen_initConfig(type=sct, dim=dim, vol=vol, nPart=nPart)
      X0 <- init_config$X0
      q  <- init_config$q
      rm(init_config)
      
      
      res <- MCMC_CPP(nIt, nPart, vol, 1/x.beta[i], sigma, N_burnIn=0, X0, q) # Make Markov Chain
      
      
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
          "q"      = I(list(res$q)),
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
    watch.records[i] <- as.double(Sys.time() - watch.iter, units='mins')
    mins <- round(mean(watch.records[c(max(c(1,i-10)),i)])*(length(x.beta)-i) , 2)
    secs <- (mins - floor(mins)) * 60
    cat("ETA= ",floor(mins)," min ",secs," s \n")
  }
  
  return(df)
}

