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
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32", 
      "C:/Users/hugjh001/Documents/len_pbpk")
    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2] | getwd() == wd[5]) {
      git.dir <- wd[2]
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
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
  dataiv.av$UID <- 1:170
  dataiv.melt <- melt(dataiv.av,
    c("TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME"),
    variable.name = "TISSUE",
    value.name = "DV"
  )
  dataiv.sdmelt <- melt(dataiv.av,
    c("TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME"),
    variable.name = "TISSUE",
    value.name = "SD"
  )
  cmt <- factor(dataiv.melt$TISSUE,
    levels = c(
      "INT", "PLA", "LUN", "LVR", "BRA", "GUT", "SPLR", "SPL", "KID",
      "FILT", "TUBC", "URINE", "HRT", "MSC", "BOD")
  )
  dataiv.obs <- data.frame(
    "Time" = dataiv.melt$TIME,
    "Concentration" = dataiv.melt$DV,
    "Error" = dataiv.sdmelt$SD,
    "Organ" = as.numeric(cmt),
    "Weight" = dataiv.melt$WT,
    "Nominal_Time" = dataiv.melt$TADNOM,
    "Dose_mgkg" = dataiv.melt$DOSEMGKG,
    "Dose_mg" = dataiv.melt$DOSEMG,
    "Route" = "Intravenous"
  )
  datapo.obs <- data.frame(
    "Time" = datapo.av$TIME,
    "Concentration" = datapo.av$DV,
    "Error" = datapo.sd$DV,
    "Organ" = 2,
    "Weight" = datapo.av$WT,
    "Nominal_Time" = datapo.av$TADNOM,
    "Dose_mgkg" = datapo.av$DOSEMGKG,
    "Dose_mg" = datapo.av$DOSEMG,
    "Route" = "Oral Gavage"
  )
  dataip.obs <- data.frame(
    "Time" = dataip.av$TIME,
    "Concentration" = dataip.av$DV,
    "Error" = dataip.sd$DV,
    "Organ" = 2,
    "Weight" = dataip.av$WT,
    "Nominal_Time" = dataip.av$TADNOM,
    "Dose_mgkg" = dataip.av$DOSEMGKG,
    "Dose_mg" = dataip.av$DOSEMG,
    "Route" = "Intraperitoneal"
  )
  units <- c("min", "ng/ml", "ng/ml", "", "kg", "min", "mg/kg", "mg", "", "", "", "")
  ave.pksim <- arrange(
    rbind(dataiv.obs, datapo.obs, dataip.obs),
    Nominal_Time
  )
  ave.pksim$Concentration[ave.pksim$Concentration == "NaN"] <- "<0.259"
  ave.pksim$Molecule <- "Lenalidomide Basic"
  ave.pksim$Species <- "Mouse"
  ave.pksim$Organ <- factor(ave.pksim$Organ)
  levels(ave.pksim$Organ) <- c("Venous Blood", "Lung", "Liver", "Brain", "Spleen",
    "Kidney", "Heart", "Muscle")
  ave.pksim$Organ <- as.character(ave.pksim$Organ)
  ave.pksim$Route <- as.character(ave.pksim$Route)
  ave.pksim$Compartment <- "Tissue"
  ave.pksim$Compartment[ave.pksim$Organ == "Venous Blood"] <- "Plasma"

  out.pksim <- rbind(units, ave.pksim)
  filename.out <- "produced_data/pksimprep.csv"
  write.csv(out.pksim, file = filename.out, row.names = F)
