# Fitting to single tissues
# -----------------------------------------------------------------------------
# Set up directories
if (!exists("git.dir")) {
  rm(list = ls(all = T))
  wd <- c("C:/Users/Jim Hughes/Documents", "C:/WINDOWS/system32",
          "C:/Users/hugjh001/Desktop", "C:/Users/hugjh001/Documents/len_pbpk")
  
  graphics.off()
  if (getwd() == wd[1]) {
    gir.dir <- paste0(getwd(), "/GitRepos")
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
library(mrgsolve)
library(splines)
library(ggplot2)

# Source functions, data and models
source(paste(git.dir, reponame, "functions",
             "utility.R", sep = "/"))
source(paste(git.dir, reponame, "rscripts",
             "data_iv.R", sep = "/"))
# loads dataiv & dataiv.av
source(paste(git.dir, reponame, "models", "single_tissues",
             "flow_limited_linear.R", sep = "/"))
source(paste(git.dir, reponame, "models", "single_tissues",
             "memb_limited_linear.R", sep = "/"))

# -----------------------------------------------------------------------------
# Ready an input data.frame for simulation
# Add ID numbers to dataiv.av
dataiv.av$ID <- factor(dataiv.av$DOSEMGKG)
levels(dataiv.av$ID) <- 1:4
dataiv.av$ID <- as.numeric(as.character(dataiv.av$ID))

# Create input data.frame but also...
# Determine coefficients for linear forcing function
fitdata <- ddply(dataiv.av, .(ID), function(x) {
  time <- c(0, x$TIME)
  pla <- c(0, x$PLA)
  tis <- c(NaN, x$HRT)
  last <- length(x$PLA)
  slope <- diff(pla)/diff(time)
  slope <- c(slope, slope[last])
  int <- pla-slope*time
  cbind(round(time, 2), round(slope, 2), round(int, 2), round(tis, 4))
})
names(fitdata)[-1] <- c("time", "M", "B", "dv")

fitdata <- merge(fitdata, fitdata, all.x = T, all.y = T)
fitdata$cmt <- 1
fitdata$M <- locf(fitdata$M)
fitdata$B <- locf(fitdata$B)

# Order by ID then time (do in reverse naturally)
fitdata <- fitdata[order(fitdata$time), ]
fitdata <- fitdata[order(fitdata$ID), ]

# -----------------------------------------------------------------------------
# Fit values for single tissues
# Fixed parameters are:
# - Q, Vreal
# Interested in fitting:
# - V (flow limited)
# - V1, V2, PS (memb limited)

# Create fitting function
fitMLE <- function(par, data, mod, sigma, fixed = NULL) {
  # Add parameters to mrgsolve input
  for (i in 1:length(names(par))) {
    data[names(par)[i]] <- abs(par[i])
  }
  # Add fixed parameters to mrgsolve input if desired
  if (!is.null(fixed)) {
    for (i in 1:length(names(fixed))) {
      data[names(fixed)[i]] <- fixed[i]
    }
  }
  # Simulate using mrgsolve
  simdata <- try(as.data.frame(mrgsim(
    data_set(mod, data)
  )))  # mrgsim
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
  flow_par <- c("V" = 0.01)
  memb_par <- c("V1" = 10, "V2" = 10, "PS" = 0.01)
  
# Run optim
  flowres <- ddply(fitdata, .(ID), function(input) {
    out <- optim(
      flow_par,
      fitMLE,
      method = "L-BFGS-B", hessian = T,
      lower = 1e-4, upper = 1e4,
      data = input, mod = flowmod, sigma = 0.01,
      fixed = c("Q" = 0.9227, "Vreal" = 0.125)
    )
    data.frame(
      model = "flowlim",
      ofv = out$value,
      V1 = out$par,
      V1se = sqrt(diag(solve(out$hessian)))/out$par*100
    )
  })
  
  indata <- fitdata
  indata$Q <- c(rep(0.9227, 38))
  indata$Vreal <- c(rep(0.125), 38)
  indata$V <- c(rep(flowres$V1[1:2], each = 9), rep(flowres$V1[3:4], each = 10))
  plotdata <- as.data.frame(mrgsim(data_set(flowmod, indata)))
  
  
  p <- NULL
  p <- ggplot(data = plotdata)
  p <- p + geom_line(aes(x = time, y = Cart),
    size = 1, alpha = 0.5, colour = "red")
  p <- p + geom_line(aes(x = time, y = Cven),
    size = 1, alpha = 0.5, colour = "blue")
  p <- p + geom_line(aes(x = time, y = Ctis),
    size = 1, alpha = 0.5, colour = "green4")
  p <- p + scale_y_log10()
  p <- p + geom_point(aes(x = time, y = dv), data = indata)
  p <- p + facet_wrap(~ID, ncol = 2, scales = "free")
  p


membres <- ddply(fitdata, .(ID), function(input) {
  out <- optim(
    memb_par,
    fitMLE,
    method = "L-BFGS-B", hessian = T,
    lower = 1e-2, upper = 1e3,
    data = input, mod = membmod, sigma = 0.1,
    fixed = c("Q" = 0.9227, "Vreal" = 0.125)
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