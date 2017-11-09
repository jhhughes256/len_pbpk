# Fitting to single tissues
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/WINDOWS/system32",
      "C:/Users/hugjh001/Desktop", "C:/Users/hugjh001/Documents/len_pbpk")

    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2] | getwd() == wd[4]) {
      git.dir <- "C:/Users/hugjh001/Documents"
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Load additional libraries
  library(deSolve)
  library(splines)
  library(ggplot2)

# Source functions, data and models
  source(paste("functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# -----------------------------------------------------------------------------
# Ready input for simulation
# Add ID numbers to dataiv.av
  dataiv.av$ID <- factor(dataiv.av$DOSEMGKG)
  levels(dataiv.av$ID) <- 1:4
  dataiv.av$ID <- as.numeric(as.character(dataiv.av$ID))

# Create input data.frame but also...
# Determine coefficients for linear forcing function
  Cart_int <- function(T, sub) {
    Cart <- approxfun(
      dataiv.av$TIME[dataiv.av$ID == sub], 
      dataiv.av$PLA[dataiv.av$ID == sub],
      method = "linear", rule = 2
    )
    Cart(T)
  }
  
# Setup models
  DESflow <- function(T, C, THETAin, sub) {
    dC <- vector(len = 1)
    
    Cart <- Cart_int(T, sub)
    dC[1] <- Q*(Cart -C[1])/V
    
    list(dC, "Cart" = Cart)
  }
  
  DESmemb <- function(T, C, THETAin, sub) {
    dC <- vector(len = 2)
    
    Cart <- Cart_int(T, sub)
    dC[1] <- (Q*(Cart  -C[1]) +PS*(C[2] -C[1]))/V1  
    dC[2] <- PS*(C[1] -C[2])/V2  
    
    list(dC, "Cart" = Cart)
  }

# -----------------------------------------------------------------------------
# Fit values for single tissues
# Fixed parameters are:
# - Q, Vreal
# Interested in fitting:
# - V (flow limited)
# - V1, V2, PS (memb limited)

# Create fitting function
  fitMLE <- function(par, times, sigma, fixed = NULL) {
    browser()
    # Simulate using deSolve
    C_0 <- c(Cven = 0, C2 = 0)
    paramlist <- c(fixed["Q"], par)
    simdata <- lsoda(C_0, data$time, DESmemb, paramlist, sub = 1)
    simdata <- data.frame(simdata)
    Vrat <- par["V"]/fixed["Vreal"]
    
    # browser()
    if(class(simdata) == "try-error") browser()
    # Remove time == 0 from data
    remove <- unique(
      which(is.nan(data$dv)), which(simdata$Ctis == 0)
    )
    dv <- data$dv[-remove]
    pred <- simdata$Ctis[-remove]
    # Determine log-likelihood
    loglik <- dnorm(dv, pred, abs(pred)*sigma, log = T)
    return(-1*sum(loglik))
  }

# Create initial parameter string to feed to function (names of parameters matter!)
  flow_par <- c(V = 0.01)
  memb_par <- c(V1 = 10, V2 = 10, PS = 0.01)

# Run optim
  # flowres <- ddply(fitdata, .(ID), function(input) {
  #   out <- optim(
  #     flow_par,
  #     fitMLE,
  #     method = "L-BFGS-B", hessian = T,
  #     lower = 1e-4, upper = 1e4,
  #     data = input, mod = flowmod, sigma = 0.01,
  #     fixed = c("Q" = 0.1573, "Vreal" = 0.0875)
  #   )
  #   data.frame(
  #     model = "flowlim",
  #     ofv = out$value,
  #     V1 = out$par,
  #     V1se = sqrt(diag(solve(out$hessian)))/out$par*100
  #   )
  # })

  membres <- ddply(dataiv.av, .(ID), function(input) {
    out <- optim(
      memb_par,
      fitMLE,
      method = "L-BFGS-B", hessian = T,
      lower = 1e-4, upper = 1e2,
      data = input, sigma = 0.1,
      fixed = c("Q" = 0.1573, "Vreal" = 0.0875)
    )
    data.frame(
      model = "memblim",
      ofv = out$value,
      V1 = out$par[1],
      V1se = sqrt(diag(solve(out$hessian)))/out$par*100[1],
      V2 = out$par[2],
      V2se = sqrt(diag(solve(out$hessian)))/out$par*100[2],
      PS = out$par[3],
      PSse = sqrt(diag(solve(out$hessian)))/out$par*100[3]
    )
  })
