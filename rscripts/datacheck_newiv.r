# Data Check Script for IV Data
# -----------------------------------------------------------------------------
# Data already published:
# Rozewski DM, Herman SEM, Towns WH, Mahoney E, Stefanovski MR, Shin JD, et al.
# Pharmacokinetics and Tissue Disposition of Lenalidomide in Mice. AAPS J.
# 2012;14(4):872-82.
# -----------------------------------------------------------------------------
# Ready workspace
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list=ls(all=TRUE))
    graphics.off()
    git.dir <- "E:/Hughes/Git"
    #git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    # git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "len_pbpk"
  }

# Load libraries
  library(readxl)
  library(plyr)
  library(reshape2)
  library(ggplot2)
  library(grid)
  library(stringr)

# Set the working directory
  master.dir <- paste(git.dir, reponame, sep = "/")
  setwd(master.dir)

# Organise working and output directories
  plot.dir <- paste(master.dir, "plot", sep = "/")
  data.dir <- paste(master.dir, "produced_data", sep = "/")
  scriptname <- "datacheck_iv"
  plot.out <- paste(plot.dir, scriptname, sep = "/")
  data.out <- paste(data.dir, scriptname, sep = "/")
  if (!file.exists(plot.out)) {
    dir.create(plot.out)
  }
  if (!file.exists(data.out)) {
    dir.create(data.out)
  }

# Source functions to be used with the script
  source(paste(git.dir, reponame, "functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "functions",
    "endsplitter.R", sep = "/"))
  source(paste(git.dir, reponame, "functions",
    "datacheck_plotfn.R", sep = "/"))

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# -----------------------------------------------------------------------------
# Load in the data
# Using read_excel due to .xls format
  file.name.in <- "E:/Hughes/Data/RAW_NonClinical/All Tissue PK data Summary - copy May4.xls"
  rawnew <- suppressWarnings(read_excel(file.name.in,
    # col_types = c(rep("text", 3), rep("numeric", 65)),
    sheet = "Data"
  ))  # read_excel
# Read in previous data (already processed)
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# -----------------------------------------------------------------------------
# Explore Data for Cleaning and Data Extraction
# View data names and structure
  names(rawnew) <- str_replace_all(names(rawnew),"[ ()#]",".")
  names(rawnew)
  str(rawnew)
  # Important columns for each tissue are weight and DV
  # Careful of column structures, many are characters

# Remove first row and store it for later use
# Contains volumes of water that were added to tissue aliquots
  firstrow <- rawnew[1, ]
  rawnew <- rawnew[-1, ]
  rawnew$DVID[is.na(rawnew$DVID)] <- 0
  rawnew <- rawnew[rawnew$DVID == 0, ]
  rawnew <- rawnew[!str_detect(rawnew$Sample.ID, "Repeat"), ]
  # While end.splitter does a good job, still need to fix repeated entry before
  # use; only need the non-repeat rows

# Subset important columns and rename columns
  rawdv <- rawnew[, str_detect(names(rawnew), "DV") & !str_detect(names(rawnew), "Norm")]
  names(rawdv) <- c("DVID", "PLA", "BRA", "LVR", "MSC", "HRT", "SPL", "LUN", "KID", rep(".", 4))
  rawwt <- rawnew[, str_detect(names(rawnew), "Wt")]
  names(rawwt) <- paste(c("", "BRA", "LVR", "MSC", "HRT", "SPL", "LUN", "KID"), "WT", sep = ".")

# Run end.splitter to get columns with which to merge with dataiv
  splitID <- end.splitter(rawnew$Sample.ID)

# Check that data is safe to merge
  length(splitID$UID) == length(dataiv$UID)
  rawwt$.WT %in% dataiv$WT
  # should be two falses at rows 113 & 114

# Bind the merging dataset together
  datanew <- cbind(splitID[1], rawwt[, -c(1, 9:10)], rawdv[, -c(1:2, 10:13)])
  # UID is the only part of splitID required for the merge
  # remove DVID, blank DV and blank weight columns
  # remove plasma DV and total weight columns as they will be provided by dataiv

# Merge with dataiv, excluding tissue DVs (columns 1:9)
  alldata <- merge(dataiv[1:9], datanew)

# -----------------------------------------------------------------------------
# Correct tissue concentrations
# Concentrations are in nano-moles per litre
# Tissue had "equal" volume of water added to the tissue for incubation before processing
# Volume was equal to the target weight of the tissue, not actual
# Collate volumes of water for incubation from text in spreadsheet (litres)
  vol_bra <- 250e-6
  vol_lvr <- 150e-6
  vol_msc <- 150e-6
  vol_hrt <- 150e-6
  vol_spl <- 70e-6
  vol_lun <- 50e-6
  vol_kid <- 150e-6

# Additionally DVs required transformation due to insufficient tissue for
# separate standard curves. Only the plasma standard curve was used.
# Ratios between known slope of neat and slope of tissue were used to transform
# Collate conversion ratios from text in spreadsheet
  rat_bra <- 0.040698/0.0130808
  rat_lvr <- 0.0171569/0.00547052  # ratio between slope of plasma (not neat) awaiting correspondence with Dolly
  rat_msc <- 0.040698/0.0632832
  rat_hrt <- 0.040698/0.0434212
  rat_spl <- 0.040698/0.0263323
  rat_lun <- 0.040698/0.0372036
  rat_kid <- 0.040698/0.0211495

# Calculations required for dv conversion
# Transform DV to account for plasma standard curve - multiply by conversion ratio
# Determine nano-moles in sample - multiply by volume of liquid added to tissue (litres)
# Determine concentration of drug in tissue - divide by mass of tissue (milligrams)
# Convert from nano-moles to nanograms - multiply by 259.26
# Desired final concentration is in ng/mg (ug/g or mg/kg)
  dv_conversion <- function(col) {
    # col - column name for conversion
    # vol - volume of liquid added to tissue
    # rat - conversion ratio to adjust for use of plasma standard curve
    dv <- suppressWarnings(as.numeric(alldata[[col]]))
    wt <- suppressWarnings(as.numeric(alldata[[paste0(col, ".WT")]]))
    vol <- get(paste0("vol_", tolower(col)))
    rat <- get(paste0("rat_", tolower(col)))
    out <- dv*rat*vol*259.26/wt
  }
  subdata <- alldata[1:8]
  subdata$PLA <- alldata$PLA*259.26/1e6  # convert from nM to ng/L to ng/uL (ug/mL or mg/L)
  subdata$BRA <- dv_conversion("BRA")
  subdata$LVR <- dv_conversion("LVR")
  subdata$MSC <- dv_conversion("MSC")
  subdata$HRT <- dv_conversion("HRT")
  subdata$SPL <- dv_conversion("SPL")
  subdata$LUN <- dv_conversion("LUN")
  subdata$KID <- dv_conversion("KID")

# -----------------------------------------------------------------------------
