library(Rcpp)
library(rgl)           # 3D Plotting
library(tidyverse)     # Remove duplicate elements

setwd(dirname(parent.frame(2)$ofile))

show_particles <- function(x, q, V){
  dim <- length(x[1,])
  l <- V**(1/dim)
  
  ipos <- which(q==+1)
  ineg <- which(q==-1)
  posParticles <- x[ipos,]
  negParticles <- x[ineg,]
  
  # pdf("test.pdf")
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
    par(pty="s")
    plot(NA, xlim=c(-l/2, l/2), ylim=c(-l/2,l/2), xlab='x-axis', ylab='y-axis', asp=1)
    
    arrows(x0=posParticles[,1]-0.5,
           y0=posParticles[,2],
           x1=posParticles[,1]+0.5,
           y1=posParticles[,2],
           length=0.01, angle=90, code=0, col='red')
    arrows(x0=posParticles[,1],
           y0=posParticles[,2]-0.5,
           x1=posParticles[,1],
           y1=posParticles[,2]+0.5,
           length=0.01, angle=90, code=0, col='red')
    arrows(x0=negParticles[,1]-0.5,
           y0=negParticles[,2],
           x1=negParticles[,1]+0.5,
           y1=negParticles[,2],
           length=0.01, angle=90, code=0, col='blue')
    arrows(x0=negParticles[,1],
           y0=negParticles[,2]-0.5,
           x1=negParticles[,1],
           y1=negParticles[,2]+0.5,
           length=0.01, angle=90, code=0, col='blue')
    
    lines(c(-l/2,-l/2,l/2,l/2,-l/2),c(-l/2,l/2,l/2,-l/2,-l/2))
  } else if (dim == 3) {
    ballradius <- 0.5
    open3d()
    plot3d(x=NA, y=NA, z=NA, xlim=c(-l/2,l/2),ylim=c(-l/2,l/2),zlim=c(-l/2,l/2), xlab="", ylab="", zlab="")
    spheres3d(x=posParticles[,1], y=posParticles[,2], z=posParticles[,3], col = "red",
              radius=ballradius)
    spheres3d(x=negParticles[,1], y=negParticles[,2], z=negParticles[,3], col = "blue",
              radius=ballradius)
  }
  # dev.off()
}

print("Load MCMC_functions.cpp.\n")

sourceCpp("MCMC_functions.cpp")
sourceCpp("MCMC_Collective.cpp")

print("Loaded MCMC_functions.R.\n")
