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

  pred.sumexp <- function(par, x, d = 0) {
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
      if (i == 1) yhat <- p[i]^d*exp(p[i]*x + p[n+i])  # first exponential defines yhat
      else if (i != n | !a) yhat <- yhat + p[i]^d*exp(p[i]*x + p[n+i])  # following exponentials add to yhat
      else if (a) yhat <- yhat - p[i]^d*exp(p[i]*x)*sum(exp(p[(n+1):(2*n-1)]))  # for absorption curve apply final term
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

  optim.sumexp <- function(data, oral = F, nexp = 3) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Set up objects
    res <- list(par = NULL, sumexp = NULL, value = NULL,
    error = NULL, hessian = NULL, message = NULL)
    # opt.par <- list(NULL)
    # opt.val <- list(NULL)
    # opt.gra <- list(NULL)
    # opt.con <- list(NULL)
    # opt.mes <- list(NULL)
    # opt.hes <- list(NULL)
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
    lmres <- -pred.lambdaz(y, x)[1]
    lm.sub <- which(y == max(y))[1]:length(y)
    repeat {
      lm.mod <- lm(log(y[lm.sub]) ~ x[lm.sub])
      lmres <- unname(lm.mod$coefficients)
      if (is.na(lmres[2])) {
        lmres <- c(max(y), -pred.lambdaz(y, x))
        break
      }
      if (lmres[2] < 0) break
      else lm.sub <- tail(lm.sub, -1)
    }
  # Estimate candidate model parameters
    for (i in 1:nexp) {
      gares <- try(ga("real-valued",
        mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
        min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
        max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      ))
      if (class(gares) == "try-error") browser()
      optres <- try(optim(
        gares@solution[1, ],
        mle.sumexp,
        method = "BFGS", hessian = T,
        x = x, y = y, sigma = 0.01
      ))
      if (class(optres) == "try-error") {
        optres <- list(
          par = gares@solution[1,],
          value = mle.sumexp(gares@solution[1,], x, y, 0.01),
          counts = c("function" = 501, gradient = NA),
          convergence = 99,
          message = "zero gradient",
          hessian = matrix(NA,
            ncol = length(gares@solution[1,]),
            nrow = length(gares@solution[1,])
          )
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

  best.sumexp.aic <- function(opt) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    k <- unlist(lapply(opt$par, length))
    aic <- ofv + 2*k
    return(sapply(opt, function(x) x[which(aic == min(aic))]))
  }

pred.lambdaz <- function(dv, t) {
    if (t[1] == 0) dv[1] <- 0
    mdv <- which(dv == 0)
    i <- 3
    bestR2 <- -1
    bestk <- 0
    if (length(dv[-mdv]) >= i) {
      repeat {
        mod <- suppressWarnings(lm(log(tail(dv[-mdv], i)) ~ tail(unique(t)[-mdv], i)))
        k <- -1*mod$coefficients["tail(unique(t)[-mdv], i)"]
        R2 <- suppressWarnings(as.numeric(summary(mod)["adj.r.squared"]))
        if (is.na(k)) browser()
        if (is.nan(R2)) R2 <- suppressWarnings(as.numeric(summary(mod)["r.squared"]))
        if (k > 0) {
          if (R2 > bestR2) {
            if (i > 2) bestR2 <- R2
            bestk <- k
          } else {
            break
          }
        }
        if (i == 5 & bestk == 0) {  #
          mdv <- c(mdv, which(dv == max(tail(dv[-mdv], 3))))
          i <- 1
        }
        if (i == length(dv[-mdv])) {  # last ditch effort (intended for simulation study only)
          if (bestk > 0) break
          else {
            mod <- suppressWarnings(lm(log(tail(dv, 2)) ~ tail(unique(t), 2)))
            bestk <- -1*mod$coefficients["tail(unique(t), 2)"]
            if (bestk > 0) break
            else {
              bestk <- log(2)/56
              break
            }
          }
        }
        i <- i + 1
      }
    }
    bestk
  }
