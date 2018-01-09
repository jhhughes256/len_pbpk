# Fitting to single tissues
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/WINDOWS/system32",
      "C:/Users/hugjh001/Desktop")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2]) {
      git.dir <- "C:/Users/hugjh001/Documents"
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Load additional libraries

# Source functions, data and models
  source(paste(git.dir, reponame, "functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# -----------------------------------------------------------------------------
# Ready an input data.frame for simulation
# Melt average DV data.frame
  meltiv.av <- melt(dataiv.av, id = c("DOSEMGKG","TADNOM", "DOSEMG", "AMT", "WT", "TIME"))
  colnames(meltiv.av) <- c(head(colnames(meltiv.av), 6), "TISSUE", "DV")

# Change all NaN to NA
  meltiv.av$DV[is.nan(meltiv.av$DV)] <- NA

# Add ID numbers for DOSEMGKG
  meltiv.av$ID <- factor(meltiv.av$DOSEMGKG)
  levels(meltiv.av$ID) <- 1:4
  meltiv.av$ID <- as.numeric(as.character(meltiv.av$ID))

# Convert TISSUE column to number format
  meltiv.av$TISSUEf <- factor(meltiv.av$TISSUE)
  levels(meltiv.av$TISSUEf) <- 1:8
  meltiv.av$TISSUE <- as.numeric(as.character(meltiv.av$TISSUEf))
  # 1 - "PLA"
  # 2 - "BRA"
  # 3 - "LVR"
  # 4 - "MSC" 
  # 5 - "HRT"
  # 6 - "SPL"
  # 7 - "LUN"
  # 8 - "KID"

# Create input data.frame but also...
# Determine coefficients for linear forcing function
  lindata <- ddply(meltiv.av[meltiv.av$TISSUE == 1,], .(ID), function(x) {
    time <- c(0, x$TIME)
    pla <- c(0, x$DV)
    last <- length(x$DV)
    slope <- diff(pla)/diff(time)
    slope <- c(slope, slope[last])
    int <- pla-slope*time
    cbind(round(time, 2), round(slope, 2), round(int, 2))
  })
  names(lindata)[-1] <- c("TIME", "M", "B")
  lindata <- data.frame(
    ID = rep(lindata$ID, 8),
    TISSUE = rep(1:8, each = 38),
    TIME = rep(lindata$TIME, 8),
    M = rep(lindata$M, 8),
    B = rep(lindata$B, 8)
  )

  ivdata <- ddply(meltiv.av, .(ID, TISSUE), function(x) {
    data.frame(
      TIME = c(0, x$TIME),
      DV = c(NA, x$DV)
    )
  })

  fitdata <- merge(ivdata, lindata,  all.x = T, all.y = T)

  fitdata <- fitdata[order(fitdata$TIME), ]
  fitdata <- fitdata[order(fitdata$ID), ]
  fitdata <- fitdata[order(fitdata$TISSUE), ]

  fitdata$M <- locf(fitdata$M)
  fitdata$B <- locf(fitdata$B)
  fitdata$MDV <- 0
  fitdata$MDV[is.na(fitdata$DV)] <- 1

  names(fitdata)[1] <- "#ID"

  write.csv(fitdata, "E:/Hughes/Data/PBPK/SINGLE_TISSUES/nmprep.csv",
    quote = FALSE, row.names = FALSE, na = ".")

# Order by ID then time (do in reverse naturally)


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
  # flowres <- ddply(fitdata, .(ID), function(input) {
  #   out <- optim(
  #     flow_par,
  #     fitMLE,
  #     method = "L-BFGS-B", hessian = T,
  #     lower = 1e-4, upper = 1e4,
  #     data = input, mod = flowmod, sigma = 0.01,
  #     fixed = c("Q" = 0.9227, "Vreal" = 0.125)
  #   )
  #   data.frame(
  #     model = "flowlim",
  #     ofv = out$value,
  #     V1 = out$par,
  #     V1se = sqrt(diag(solve(out$hessian)))/out$par*100
  #   )
  # })

  membres <- ddply(fitdata, .(ID), function(input) {
    out <- optim(
      memb_par,
      fitMLE,
      method = "L-BFGS-B", hessian = T,
      lower = 1e-4, upper = 1e2,
      data = input, mod = membmod, sigma = 0.1,
      fixed = c("Q" = 0.9227, "Vreal" = 0.125)
    )
    browser()
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
