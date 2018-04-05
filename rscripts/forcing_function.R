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
    git.dir <- "E:/Hughes/Git"
    # git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    reponame <- "len_pbpk"
  }
  # Load libraries
  library(reshape2)
  library(GA)
  library(plyr)

  # Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# Set up sum of exponential function
  order.sumexp <- function(par, n, a) {
  # Sorts sum of exponential parameters to enable pred.sumexp
  # par = sum of exponential parameters
  # n = number of exponentials
  # a = absorption status
    m <- -abs(head(par, n+a))  # slope parameters (prevents exponential growth)
    b <- tail(par, -(n+a))  # intercept parameters
    m.ord <- order(m, decreasing = T)  # slope order (terminal first, absorption last)
    b.ord <- m.ord  # intercept order (match slopes)
    if (a) b.ord <- order(m[-m.ord[n+a]], decreasing = T)  # if absorption curve remove extra term
    unname(c(m[m.ord], b[b.ord]))  # ordered parameters
  }

  pred.sumexp <- function(par, x) {
  # Provides predictions according to model parameters
  # par = sum of exponential parameters
  # x = independent variable (time)
  # d = order of derivative (uses dth derivative of model)
  # Define objects
    l <- length(par)  # number of parameters
    a <- l %% 2 == 1  # absorption status (odd parameter length == absorption)
    n <- ceiling(l/2)  # number of exponentials
    m <- -abs(par[1:n])  # slope parameters (prevents exponential growth)
    b <- par[(n+1):l]  # intercept parameters
  # Order parameters (allows for flip-flop)
    m.ord <- order(m, decreasing = T)  # slope order (terminal first, absorption last)
    b.ord <- m.ord  # intercept order (match slopes)
    if (a) b.ord <- order(m[-m.ord[n]], decreasing = T)  # if absorption curve remove extra term
    p <- c(m[m.ord], b[b.ord])  # ordered parameters
  # Sum of exponentials
    for (i in 1:n) {  # for each exponential
      if (i == 1) yhat <- exp(p[i]*x + p[n+i])  # first exponential defines yhat
      else if (i != n | !a) yhat <- yhat + exp(p[i]*x + p[n+i])  # following exponentials add to yhat
      else if (a) yhat <- yhat - exp(p[i]*x)*sum(exp(p[(n+1):(2*n-1)]))  # for absorption curve apply final term
    }
    return(yhat)  # predicted dependent variable (drug concentration)
  }

  mle.sumexp <- function(par, x, y, sigma, ga = F) {
  # Determine log likelihood of given model parameters
  # par = sum of exponential parameters
  # x = independent variable (time)
  # y = observed dependent variable (drug concentration)
  # sigma = proportional error
  # ga = genetic algorithm status
    z <- ifelse(ga, 2, -2)  # adjust ofv for minimisation or maximisation
    yhat <- pred.sumexp(par, x)  # sum of exponential model prediction
    loglik <- dnorm(y, yhat, abs(yhat)*sigma, log = T)  # log likelihood
    return(z*sum(loglik))  # objective function value
  }

  pred.lambdaz <- function(dv, t) {
    i <- 2
    j <- 0
    bestR2 <- -1
    bestk <- rep(0, 3)
    terminal <- which(dv == max(dv))[1]:length(dv)
    rem <- matrix(length(dv)+1)
    if (length(terminal) >= i) {
      repeat {
        for (l in 1:ncol(rem)) {
          mod <- suppressWarnings(lm(log(tail(dv[-rem[,l]], i)) ~ tail(unique(t)[-rem[,l]], i)))
          k <- unname(mod$coefficients)
          R2 <- suppressWarnings(as.numeric(summary(mod)["adj.r.squared"]))
          if (!is.finite(k[2])) {
            R2 <- -1
            break
          }
          if (is.nan(R2)) R2 <- suppressWarnings(as.numeric(summary(mod)["r.squared"]))
          if (k[2] < 0) {
            if (R2 > bestR2) {
              if (i > 2) bestR2 <- R2
              bestk <- c(k, sd(residuals(mod)))
            } else {
              break  # break out of for loop
            }
          }
        }
        if (R2 < bestR2) break  # take model that is better than latest model
        if (i == length(terminal)) {  # once all points have been included test for model validity
          if (!is.finite(bestk[2])) { browser() }
          else if (bestk[2] == 0) {  # if all current models are invalid, trial removing random points
            if (length(terminal) != j + 3) {  # as long as there are still points to remove
              j <- j + 1
              i <- j + 2
              rem <- length(dv) - combn(i + 1, j) + 1
            } else {  # if not make a last ditch effort (intended for simulation study only)
              bestk <- c(max(dv), -log(2)/56, 0.01)
              break  # break out of repeat loop
            }
          } else { break }  # if there was a valid model then finish
        }
        i <- i + 1
      }
    }
    bestk
  }

  optim.sumexp.new <- function(data, oral = F, nexp = 4) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Set up objects
    res <- list(par = NULL, sumexp = NULL, value = NULL,
    error = NULL, hessian = NULL, message = NULL)
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
  # Estimate candidate model parameters
    init <- pred.lambdaz(y, x)
    for (i in 1:nexp) {
      gares <- try(ga("real-valued",
        mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
        min = c(rep(init[2]*50, i + oral), rep(init[1]-2, i)),
        max = c(rep(init[2]/50, i + oral), rep(init[1]+2, i)),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      ))
      if (class(gares) == "try-error") browser()
      ga.par <- gares@solution[1, ]
      optres <- try(optim(
        ga.par,
        mle.sumexp,
        method = "BFGS", hessian = T,
        x = x, y = y, sigma = 0.01
      ))
      if (class(optres) == "try-error") {
        optres <- list(
          par = gares@solution[1,],
          value = -gares@fitnessValue,
          counts = NA,
          hessian = matrix(NA, ncol = length(ga.par), nrow = length(ga.par)),
          convergence = NA,
          message = NA
        )
      }
    # Create output
      par.ord <- order.sumexp(optres$par, i, oral)
      res$par[[i]] <- optres$par
      res$sumexp[[i]] <- par.ord
      res$value[[i]] <- c(ofv = optres$value, optres$counts)
      res$error[[i]] <- c(0.01, "fixed")
      res$hessian[[i]] <- optres$hessian
      res$message[[i]] <- c(convergence = optres$convergence,
        message = ifelse(is.null(optres$message), "NULL", optres$message)
      )
    }
    res
  }

# Grab plasma data
# Convert data to units PKSim is expecting
# Units are ng/mL aka ug/L
# convert to umol/L as when multiplied by volume in the arterial blood it gives
# umol which is PKSim's native unit
  dataiv.pla <- with(dataiv[!is.na(dataiv$PLA),], data.frame(
    TADNOM = TADNOM,
    DOSEMGKG = DOSEMGKG,
    TIME = TIME,
    DV = PLA/259.26,
    DVNORM = PLA/(DOSEMG*259.26)
  ))

  dataiv.av <- ddply(dataiv.pla, .(TADNOM, DOSEMGKG), function(x) {
    with(x, data.frame(
      TIME = mean(TIME, na.rm = T),
      DV = mean(DV, na.rm = T),
      DVNORM = mean(DVNORM, na.rm = T)
    ))
  })

# Fit using DVNORM data then check dose transformed predictions against original data
  out05 <- optim.sumexp.new(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 0.5],
    y = dataiv.av$DVNORM[dataiv.av$DOSEMGKG == 0.5]
  ))
  plot(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 0.5],
    y = log(dataiv.av$DV[dataiv.av$DOSEMGKG == 0.5])
  ))
  lines(data.frame(
    x = seq(0, 600, by = 1),
    y = log(0.015*pred.sumexp(out05$sumexp[[3]], seq(0, 600, 1)))
  ))

  out15 <- optim.sumexp.new(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 1.5],
    y = dataiv.av$DVNORM[dataiv.av$DOSEMGKG == 1.5]
  ))
  plot(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 1.5],
    y = log(dataiv.av$DV[dataiv.av$DOSEMGKG == 1.5])
  ))
  lines(data.frame(
    x = seq(0, 600, by = 1),
    y = log(0.045*pred.sumexp(out15$sumexp[[3]], seq(0, 600, 1)))
  ))

  out5 <- optim.sumexp.new(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 5],
    y = dataiv.av$DVNORM[dataiv.av$DOSEMGKG == 5]
  ))
  plot(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 5],
    y = log(dataiv.av$DV[dataiv.av$DOSEMGKG == 5])
  ))
  lines(data.frame(
    x = seq(0, 600, by = 1),
    y = log(0.15*pred.sumexp(out5$sumexp[[4]], seq(0, 600, 1)))
  ))

  out10 <- optim.sumexp.new(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 10],
    y = dataiv.av$DVNORM[dataiv.av$DOSEMGKG == 10]
  ))
  plot(data.frame(
    x = dataiv.av$TIME[dataiv.av$DOSEMGKG == 10],
    y = log(dataiv.av$DV[dataiv.av$DOSEMGKG == 10])
  ))
  lines(data.frame(
    x = seq(0, 600, by = 1),
    y = log(0.3*pred.sumexp(out10$sumexp[[4]], seq(0, 600, 1)))
  ))
