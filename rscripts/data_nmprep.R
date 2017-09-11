# NONMEM Preparation Script for Observed Data
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
# Load libraries
  library(reshape2)

# Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_po.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_ip.R", sep = "/"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataiv$UID <- 1:170
  dataiv.melt <- melt(dataiv,
    c("UID", "ID", "TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME"),
    variable.name = "TISSUE",
    value.name = "DV"
  )
  cmt <- factor(dataiv.melt$TISSUE,
    levels = c(
      "INT", "PLA", "LUN", "LVR", "BRA", "GUT", "SPLR", "SPL", "KID",
      "FILT", "TUBC", "URINE", "HRT", "MSC", "BOD")
  )
  dataiv.dv <- data.frame(
    CID = dataiv.melt$UID,
    TIME = dataiv.melt$TIME,
    AMT = NA,
    EVID = 0,
    DV = dataiv.melt$DV,
    CMT = as.numeric(cmt),
    MDV = 0,
    WT = dataiv.melt$WT,
    TADNOM = dataiv.melt$TADNOM,
    DOSEMGKG = dataiv.melt$DOSEMGKG,
    DOSE = dataiv.melt$DOSEMG,
    ROUTE = 1
  )
  dataiv.dv$MDV[is.na(dataiv.dv$DV)] <- 1
  dataiv.dose <- data.frame(
    CID = dataiv$UID,
    TIME = 0,
    AMT = dataiv$AMT,
    EVID = 1,
    DV = NA,
    CMT = 2,
    MDV = 1,
    WT = dataiv$WT,
    TADNOM = dataiv$TADNOM,
    DOSEMGKG = dataiv$DOSEMGKG,
    DOSE = dataiv$DOSEMG,
    ROUTE = 1
  )
  datapo.dv <- data.frame(
    CID = datapo$UID+170,
    TIME = datapo$TIME,
    AMT = NA,
    EVID = 0,
    DV = datapo$DV,
    CMT = 2,
    MDV = 0,
    WT = datapo$WT,
    TADNOM = datapo$TADNOM,
    DOSEMGKG = datapo$DOSEMGKG,
    DOSE = datapo$DOSEMG,
    ROUTE = 2
  )
  datapo.dose <- data.frame(
    CID = datapo$UID+170,
    TIME = 0,
    AMT = datapo$AMT,
    EVID = 1,
    DV = NA,
    CMT = 1,
    MDV = 1,
    WT = datapo$WT,
    TADNOM = datapo$TADNOM,
    DOSEMGKG = datapo$DOSEMGKG,
    DOSE = datapo$DOSEMG,
    ROUTE = 2
  )
  dataip.dv <- data.frame(
    CID = dataip$UID+261,
    TIME = dataip$TIME,
    AMT = NA,
    EVID = 0,
    DV = dataip$DV,
    CMT = 2,
    MDV = 0,
    WT = dataip$WT,
    TADNOM = dataip$TADNOM,
    DOSEMGKG = dataip$DOSEMGKG,
    DOSE = dataip$DOSEMG,
    ROUTE = 3
  )
  dataip.dose <- data.frame(
    CID = dataip$UID+261,
    TIME = 0,
    AMT = dataip$AMT,
    EVID = 1,
    DV = NA,
    CMT = 4,
    MDV = 1,
    WT = dataip$WT,
    TADNOM = dataip$TADNOM,
    DOSEMGKG = dataip$DOSEMGKG,
    DOSE = dataip$DOSEMG,
    ROUTE = 3
  )
  data.nm <- arrange(
    rbind(dataiv.dv, dataiv.dose, datapo.dv, datapo.dose, dataip.dv, dataip.dose),
    CID, TIME
  )
  data.nm$AMT[is.na(data.nm$AMT)] <- "."
  data.nm$MDV[is.na(data.nm$DV)] <- 1
  data.nm$DV[is.na(data.nm$DV)] <- "."
  names(data.nm)[1] <- "#ID"
  filename.out <- "produced_data/nmprep.csv"
  write.csv(data.nm, file = filename.out, quote = F, row.names = F)
