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
  dataiv.obs <- data.frame(
    "Time" = dataiv.melt$TIME,
    "Concentration" = dataiv.melt$DV,
    "Organ" = as.numeric(cmt),
    "Weight" = dataiv.melt$WT,
    "Nominal_Time" = dataiv.melt$TADNOM,
    "Dose_mgkg" = dataiv.melt$DOSEMGKG,
    "Dose_mg" = dataiv.melt$DOSEMG,
    "Route" = "Intravenous"
  )
  datapo.obs <- data.frame(
    "Time" = datapo$TIME,
    "Concentration" = datapo$DV,
    "Organ" = 2,
    "Weight" = datapo$WT,
    "Nominal_Time" = datapo$TADNOM,
    "Dose_mgkg" = datapo$DOSEMGKG,
    "Dose_mg" = datapo$DOSEMG,
    "Route" = "Oral Gavage"
  )
  dataip.obs <- data.frame(
    "Time" = dataip$TIME,
    "Concentration" = dataip$DV,
    "Organ" = 2,
    "Weight" = dataip$WT,
    "Nominal_Time" = dataip$TADNOM,
    "Dose_mgkg" = dataip$DOSEMGKG,
    "Dose_mg" = dataip$DOSEMG,
    "Route" = "Intraperitoneal"
  )
  units <- c("min", "ng/ml", "", "kg", "min", "mg/kg", "mg", "", "ng/ml", "", "", "")
  all.pksim <- arrange(
    rbind(dataiv.obs, datapo.obs, dataip.obs),
    Nominal_Time
  )
  all.pksim$DoseNormConc <- all.pksim$Concentration/all.pksim$Dose_mg

  all.pksim$Concentration[is.na(all.pksim$Concentration)] <- "<0.259"
  loq <- paste0("<", 0.259/all.pksim$Dose_mg[is.na(all.pksim)])
  all.pksim$DoseNormConc[is.nan(all.pksim$DoseNormConc)] <- loq

  all.pksim$Molecule <- "Lenalidomide Basic"
  all.pksim$Species <- "Mouse"
  all.pksim$Organ <- factor(all.pksim$Organ)
  levels(all.pksim$Organ) <- c("Venous Blood", "Lung", "Liver", "Brain", "Spleen",
    "Kidney", "Heart", "Muscle")
  all.pksim$Organ <- as.character(all.pksim$Organ)
  all.pksim$Route <- as.character(all.pksim$Route)
  all.pksim$Compartment <- "Tissue"
  all.pksim$Compartment[all.pksim$Organ == "Venous Blood"] <- "Plasma"

  all.pksim <- all.pksim[order(all.pksim$Time), ]
  out.pksim <- rbind(units, all.pksim)
  filename.out <- "produced_data/pksimprep_all.csv"
  write.csv(out.pksim, file = filename.out, row.names = F)
