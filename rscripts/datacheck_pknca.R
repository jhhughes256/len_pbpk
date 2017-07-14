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

  library(PKNCA)
  library(plyr)
  library(reshape2)
  data_conc <- data.frame(
    DOSEMGKG = c(dataiv.av$DOSEMGKG, 0.5, 10, datapo.av$DOSEMGKG),
    TADNOM = c(dataiv.av$TADNOM, 0, 0, datapo.av$TADNOM),
    TIME = c(dataiv.av$TIME, 0, 0, datapo.av$TIME),
    AMT = c(dataiv.av$AMT, 14000, 280000, datapo.av$AMT),
    DV = c(dataiv.av$PLA, 0, 0, datapo.av$DV),
    WT = c(dataiv.av$WT, 28, 28, datapo.av$WT),
    route = c(
      rep("intravascular", dim(dataiv.av)[1]),
      rep("extravascular", dim(datapo.av)[1]+2)
    )
  )
  data_conc$DV[is.nan(data_conc$DV)] <- NA
  data_conc$TIME[data_conc$TADNOM == 2 & data_conc$route == "intravascular"] <- 0
  nca_conc <- PKNCAconc(data_conc, DV~TIME|DOSEMGKG+route)
  data_dose <- ddply(data_conc, .(DOSEMGKG, route), function(x) {
    out <- x[1,]
    out$AMT <- mean(x$AMT)
    out$DOSEMG <- mean(x$DOSEMG)
    out$WT <- mean(x$WT)
    out
  })
  nca_dose <- PKNCAdose(data_dose, AMT~TIME|DOSEMGKG+route, route = "intravascular")
  ncadata <- PKNCAdata(nca_conc, nca_dose)
  resultnca.melt <- pk.nca(ncadata)$result
  resultnca <- dcast(resultnca.melt, route+DOSEMGKG ~ PPTESTCD, value.var = "PPORRES")
  resultnca$dose <- arrange(data_dose, route, DOSEMGKG)$AMT
  resultnca$wt <- arrange(data_dose, route, DOSEMGKG)$WT
  resultnca$cl <- with(resultnca, dose/aucinf.obs)

http://onlinelibrary.wiley.com/doi/10.1002/jps.2600820718/pdf

  resultncaiv[resultncaiv$PPTESTCD == "half.life",]
  resultncaiv[resultncaiv$PPTESTCD == "aucinf.obs",]
