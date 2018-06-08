# PKSim Preparation Script for Observed Data
# -----------------------------------------------------------------------------
# Data already published:
# Rozewski DM, Herman SEM, Towns WH, Mahoney E, Stefanovski MR, Shin JD, et al.
# Pharmacokinetics and Tissue Disposition of Lenalidomide in Mice. AAPS J.
# 2012;14(4):872-82.
# -----------------------------------------------------------------------------
# Ready workspace
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list=ls(all=TRUE))
    graphics.off()
    # git.dir <- "E:/Hughes/Git"
    # git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "len_pbpk"
  }
  # Load libraries
  library(reshape2)
  library(GA)
  library(plyr)

  # Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_po.R", sep = "/"))

# Fit fourier curve to mean pooled oral data
# data object is datapo.av
  subdata <- datapo.av[datapo.av$DOSEMGKG == 0.5 & !is.nan(datapo.av$DV),]
# Add concentration data at zero for fitting
  firstrow <- subdata[1, ]
  firstrow$TADNOM <- 0
  firstrow$TIME <- 0
  firstrow$DV <- 0
  subdata <- rbind(firstrow, subdata)
  t <- subdata$TIME
  y <- subdata$DV
  sigma <- 0.1
  pr <- 2 # number of pairs
  n <- length(t) # sample size
  tmp <- rep(NA, 2*pr*n) #construction of a vector to be reshaped to a matrix
  col_p <- matrix(tmp, n, 2*pr) #reshape the vector tmp to a matrix with
  # n rows and 2*pr columns
  for(k in 1:pr) { # create pr columns of cosines followed by
    # pr columns of sines
    for(j in 1:n) {
      col_p[j,k] <- cos(2*pi*t[j]*k/1440)
      col_p[j, pr+k] <- sin(2*pi*t[j]*k/1440)
    }
  }
  ones <- rep(1,n)
  X <- cbind(ones,col_p)
  fit <- lm(y ~ 0 + X)
  coef <- unname(fit$coefficients)
  yhat <- (X%*%coef)[,1]
  
  ts <- seq(0, 180, 2)
  n <- length(ts) # sample size
  tmp <- rep(NA, 2*pr*n) #construction of a vector to be reshaped to a matrix
  col_p <- matrix(tmp, n, 2*pr) #reshape the vector tmp to a matrix with
  # n rows and 2*pr columns
  for(k in 1:pr) { # create pr columns of cosines followed by
    # pr columns of sines
    for(j in 1:n) {
      col_p[j,k] <- cos(2*pi*ts[j]*k/1440)
      col_p[j, pr+k] <- sin(2*pi*ts[j]*k/1440)
    }
  }
  ones <- rep(1,n)
  Xhat <- cbind(ones,col_p)
  yplot <- (Xhat%*%coef)[,1]
  with(subdata, plot(DV ~ TIME))
  lines(yplot ~ ts)
  loglik <- dnorm(y, yhat, sigma, log = T)
  2*length(coef)-2*sum(loglik)