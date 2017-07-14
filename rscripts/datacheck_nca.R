# Data Check Script for PO Data
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
    #git.dir <- "E:/Hughes/Git"
    git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    reponame <- "len_pbpk"
  }

# Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_load.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_po.R", sep = "/"))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load NCA functions
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
  auc.interv <- function(times, fit.par, fn, log = F) {
    C <- do.call(fn, list(x = times, p = fit.par))
    auc <- c(NULL)
    for (i in 2:length(C)) {
      h <- times[i] - times[i-1]
      dC <- C[i-1] - C[i]
      if (log & dC > 0) auc[i-1] <- dC*h/log(C[i-1]/C[i])
      else auc[i-1] <- (C[i-1] + C[i])*h/2
    }
    return(sum(auc))
  }
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load NCA functions
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# https://www.boomer.org/c/p4/c20/c2001.html
  nca_calc <- function(times, dv, dose, method = "linlog") {
    # define data variables
    not_na <- which(!is.na(dv))
    times <- times[not_na]
    dv <- dv[not_na]
    dvt <- dv*times
    # frequently used variable calculations
    l <- length(times)
    h <- times[-1] - times[-l]
    dC <- dv[-l] - dv[-1]
    dCt <- dvt[-l] - dvt[-1]
    # determine the slope of the terminal concentrations
    best <- list("r2" = 0, "k" = 0, "n" = 0)
    diff_ldv_maxdv <- length(dv) - which(dv == max(dv))
    if (diff_ldv_maxdv < 3) {
      stop("Less than three terminal concentrations")
    }
    for (i in 3:diff_ldv_maxdv) {
      fit <- lm(log(tail(dv, i)) ~ tail(times, i))
      k <- -fit$coefficients["tail(times, i)"]
      r2 <- as.numeric(summary(fit)["adj.r.squared"])
      if (r2 > (best$r2 - 0.0001) & k > 0) {
        best$r2 <- r2
        best$k <- unname(k)
        best$n <- i
      } else if (best$k != 0) {
        break
      }
    }
    if (best$k == 0) stop("Terminal slope best fit by positive rate constant")
    if (diff(times[c((l-best$n+1),l)]) < 2*log(2)/best$k) {
      warning("Terminal slope does not cover more than two half-lives")
    }
    # determine final parameters
    vec_auc <- (dv[-1] + dv[-l])*h/2
    vec_aumc <- (dvt[-1] + dvt[-l])*h/2
    if (method == "linlog") {
      neg <- which(dC < 0)
      negt <- which(dCt < 0)
      sum_auc <- sum(c(vec_auc[neg], dC[-neg]*h[-neg]/log((dv[-l]/dv[-1])[-neg])))
      sum_aumc <- sum(c(vec_aumc[negt], dCt[-negt]*h[-negt]/log((dvt[-l]/dvt[-1])[-negt])))
    } else if (method == "linear") {
      sum_auc <- sum(vec_auc)
      sum_aumc <- sum(vec_aumc)
    } else if (method == "log") {
      sum_auc <- sum(dC*h/log((dv[-l]/dv[-1])))
      sum_aumc <- sum(dCt*h/log((dvt[-l]/dvt[-1])))
    } else {
      stop(paste(method, "is not a valid option for method"))
    }
    term_auc <- dv[l]/best$k
    term_aumc <- dvt[l]/best$k + dv[l]/best$k^2
    auc <- sum_auc + term_auc
    aumc <- sum_aumc + term_aumc
    mrt <- aumc/auc
    ke <- 1/mrt
    cl <- dose/auc
    vss <- cl*mrt
    return(c(
      "auc_0_24" = sum_auc,
      "auc_24_inf" = term_auc,
      "auc_0_inf" = auc,
      "aumc_0_24" = sum_aumc,
      "aumc_24_inf" = term_aumc,
      "aumc_0_inf" = aumc,
      "mrt" = mrt,
      "kel" = ke,
      "CL" = cl,
      "Vss" = vss,
      "dose" = dose
    ))  # return
  }
  po_nca_calc <- function(iv_nca, po_nca) {
    mat <- unname(po_nca["mrt"] - iv_nca["mrt"])
    ka <- 1/mat
    f1 <- unname(po_nca["auc_0_inf"]*iv_nca["dose"]/(iv_nca["auc_0_inf"]*po_nca["dose"]))
    names(iv_nca) <- paste0("iv_", names(iv_nca))
    names(po_nca) <- paste0("po_", names(po_nca))
    return(c(
      iv_nca,
      po_nca,
      "po_mat" = mat,
      "po_ka" = ka,
      "po_f1" = f1
    ))
  }
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load NCA functions
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
  ivsub05 <- dataiv.av[dataiv.av$DOSEMGKG == 0.5,]
  ivsub05.nca <- nca_calc(c(0, ivsub05$TIME), c(0, ivsub05$PLA), mean(ivsub05$AMT))
  ivsub15 <- dataiv.av[dataiv.av$DOSEMGKG == 1.5,]
  ivsub15.nca <- nca_calc(c(0, ivsub15$TIME), c(0, ivsub15$PLA), mean(ivsub15$AMT))
  ivsub5 <- dataiv.av[dataiv.av$DOSEMGKG == 5,]
  ivsub5.nca <- nca_calc(c(0, ivsub5$TIME), c(0, ivsub5$PLA), mean(ivsub5$AMT))
  ivsub10 <- dataiv.av[dataiv.av$DOSEMGKG == 10,]
  ivsub10.nca <- nca_calc(c(0, ivsub10$TIME), c(0, ivsub10$PLA), mean(ivsub10$AMT))
  posub05 <- datapo.av[datapo.av$DOSEMGKG == 0.5,]
  posub05.nca <- nca_calc(c(0, posub05$TIME), c(0, posub05$DV), mean(posub05$AMT))
  posub10 <- datapo.av[datapo.av$DOSEMGKG == 10,]
  posub10.nca <- nca_calc(c(0, posub10$TIME), c(0, posub10$DV), mean(posub10$AMT))

  mean(c(
    po_nca_calc(ivsub05.nca, posub05.nca)["po_ka"],
    po_nca_calc(ivsub15.nca, posub05.nca)["po_ka"],
    po_nca_calc(ivsub5.nca, posub05.nca)["po_ka"],
    po_nca_calc(ivsub10.nca, posub05.nca)["po_ka"],
    po_nca_calc(ivsub05.nca, posub10.nca)["po_ka"],
    po_nca_calc(ivsub15.nca, posub10.nca)["po_ka"],
    po_nca_calc(ivsub5.nca, posub10.nca)["po_ka"],
    po_nca_calc(ivsub10.nca, posub10.nca)["po_ka"]
  ))

  c(
    ivsub05.nca["kel"],
    ivsub15.nca["kel"],
    ivsub5.nca["kel"],
    ivsub10.nca["kel"],
    posub05.nca["kel"],
    posub10.nca["kel"]
  )
