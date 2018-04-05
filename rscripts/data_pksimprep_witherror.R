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
# IV Data
  dataiv.melt <- melt(dataiv[, -(1:2)],
    c("TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME"),
    variable.name = "TISSUE",
    value.name = "DV"
  )
  dataiv.melt$DVNORM <- with(dataiv.melt, DV/DOSEMG)
  dataiv.avsd <- ddply(dataiv.melt, .(TADNOM, DOSEMGKG, TISSUE), function(x) {
    data.frame(
      DOSEMG = round(mean(x$DOSEMG, na.rm = T), 3),
      AMT = round(mean(x$AMT, na.rm = T), 3),
      WT = round(mean(x$WT, na.rm = T), 3),
      DV = round(mean(x$DV, na.rm = T), 3),
      SD = round(sd(x$DV, na.rm = T), 3),
      DVNORM = round(mean(x$DVNORM, na.rm = T), 3),
      SDNORM = round(sd(x$DVNORM, na.rm = T), 3)
    )
  })
  cmt <- factor(dataiv.avsd$TISSUE,
    levels = c(
      "INT", "PLA", "LUN", "LVR", "BRA", "GUT", "SPLR", "SPL", "KID",
      "FILT", "TUBC", "URINE", "HRT", "MSC", "BOD")
  )
  dataiv.obs <- data.frame(
    "Time" = dataiv.avsd$TADNOM,
    "Dose_mgkg" = dataiv.avsd$DOSEMGKG,
    "Dose_mg" = dataiv.avsd$DOSEMG,
    "Concentration" = dataiv.avsd$DV,
    "Error" = dataiv.avsd$SD,
    # "Concentration" = dataiv.avsd$DVNORM,
    # "Error" = dataiv.avsd$SDNORM,
    "Organ" = as.numeric(cmt),
    "Weight" = dataiv.avsd$WT,
    "Route" = "Intravenous"
  )
# PO Data
  datapo$DVNORM <- with(datapo, DV/DOSEMG)
  datapo.avsd <- ddply(datapo, .(TADNOM, DOSEMGKG), function(x) {
    data.frame(
      DOSEMG = round(mean(x$DOSEMG, na.rm = T), 3),
      AMT = round(mean(x$AMT, na.rm = T), 3),
      WT = round(mean(x$WT, na.rm = T), 3),
      DV = round(mean(x$DV, na.rm = T), 3),
      SD = round(sd(x$DV, na.rm = T), 3),
      DVNORM = round(mean(x$DVNORM, na.rm = T), 3),
      SDNORM = round(sd(x$DVNORM, na.rm = T), 3)
    )
  })
  datapo.obs <- data.frame(
    "Time" = datapo.avsd$TADNOM,
    "Dose_mgkg" = datapo.avsd$DOSEMGKG,
    "Dose_mg" = datapo.avsd$DOSEMG,
    "Concentration" = datapo.avsd$DV,
    "Error" = datapo.avsd$SD,
    # "Concentration" = datapo.avsd$DVNORM,
    # "Error" = datapo.avsd$SDNORM,
    "Organ" = 2,
    "Weight" = datapo.avsd$WT,
    "Route" = "Oral Gavage"
  )
# IP Data
  dataip$DVNORM <- with(dataip, DV/DOSEMG)
  dataip.avsd <- ddply(dataip, .(TADNOM, DOSEMGKG), function(x) {
    data.frame(
      DOSEMG = round(mean(x$DOSEMG, na.rm = T), 3),
      AMT = round(mean(x$AMT, na.rm = T), 3),
      WT = round(mean(x$WT, na.rm = T), 3),
      DV = round(mean(x$DV, na.rm = T), 3),
      SD = round(sd(x$DV, na.rm = T), 3),
      DVNORM = round(mean(x$DVNORM, na.rm = T), 3),
      SDNORM = round(sd(x$DVNORM, na.rm = T), 3)
    )
  })
  dataip.obs <- data.frame(
    "Time" = dataip.avsd$TADNOM,
    "Dose_mgkg" = dataip.avsd$DOSEMGKG,
    "Dose_mg" = dataip.avsd$DOSEMG,
    "Concentration" = dataip.avsd$DV,
    "Error" = dataip.avsd$SD,
    # "Concentration" = dataip.avsd$DVNORM,
    # "Error" = dataip.avsd$SDNORM,
    "Organ" = 2,
    "Weight" = dataip.avsd$WT,
    "Route" = "Intraperitoneal"
  )
  units <- c("min", "mg/kg", "mg", "ng/ml", "ng/ml", "", "g", "", "", "", "")
  ave.pksim <- arrange(
    rbind(dataiv.obs, datapo.obs, dataip.obs),
    Time
  )
  lloq.conc <- 0.259
  lloq.sd <- 0.259
  # lloq.conc <- 0.259/ave.pksim$Dose_mg[ave.pksim$Concentration == "NaN"]
  # lloq.sd <- 0.259/ave.pksim$Dose_mg[ave.pksim$Error %in% c(NA, "NaN")]
  ave.pksim$Concentration[ave.pksim$Concentration == "NaN"] <- paste0("<", lloq.conc)
  ave.pksim$Error[ave.pksim$Error %in% c(NA, "NaN")] <- paste0("<", lloq.sd)
  ave.pksim$Molecule <- "Lenalidomide"
  ave.pksim$Species <- "Mouse"
  ave.pksim$Organ <- factor(ave.pksim$Organ)
  levels(ave.pksim$Organ) <- c("Venous Blood", "Lung", "Liver", "Brain", "Spleen",
    "Kidney", "Heart", "Muscle")
  ave.pksim$Organ <- as.character(ave.pksim$Organ)
  ave.pksim$Route <- as.character(ave.pksim$Route)
  ave.pksim$Compartment <- "Tissue"
  ave.pksim$Compartment[ave.pksim$Organ == "Venous Blood"] <- "Plasma"

  out.pksim <- rbind(units, ave.pksim)
  # filename.out <- "produced_data/pksimprep_dvnorm.csv"
  filename.out <- "produced_data/pksimprep_dv.csv"
  write.csv(out.pksim, file = filename.out, row.names = F)
