# Data Check Script for IP Data
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

  scriptname <- "datacheck_ip"
  plot.out <- paste(plot.dir, scriptname, sep = "/")
  data.out <- paste(data.dir, scriptname, sep = "/")
  if (!file.exists(plot.out)) {
    dir.create(plot.out)
  }
  if (!file.exists(data.out)) {
    dir.create(data.out)
  }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explore Data for Cleaning and Data Extraction
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# View data names and structure
  names(rawip)
  sort(names(rawip))
  str(rawip)

# Check distribution of doses and weights
  with(rawip, table(Dose..mg.kg., useNA = "always"))
  with(rawip, table(Mouse.Wt..g., useNA = "always"))
  # 39 NA weights

# Check NA's for non-binned data
  any(is.na(rawip$Sample.ID))
  any(is.na(rawip$Time..min.))

# Isolate NA's
  rawip[which(is.na(rawip$Mouse.Wt..g.)), ]
  # NA weights are from 0.5 mg/kg dose level

# Number of samples
  length(unique(rawip$Sample.ID))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Clean & Extract Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Convert data for all data
  dataip <- data.frame("ID" = rawip$Sample.ID, "DOSEMGKG" = rawip$Dose..mg.kg.)
  dataip$DOSEMG <- dataip$DOSEMGKG*rawip$Mouse.Wt..g./1000
  dataip$AMT <- dataip$DOSEMG*10^6	#dose in ng
  dataip$WT <- rawip$Mouse.Wt..g.
  dataip$TIME <- rawip$Time..min.
  dataip$DV <- rawip$Plasma.DV..ng.mL.

# Clean sample IDs
  IDiv <- end.splitter(dataip$ID)
  dataip <- cbind(IDiv, dataip[,-1])
  # Remove the ID column and bind the new ID columns in its place

# Replace missing WT information with median value (28 g)
  dataip$WT[dataip$DOSEMGKG == 0.5] <- 28
  dataip$DOSEMG[dataip$DOSEMGKG == 0.5] <- 0.5*28/1000
  dataip$AMT[dataip$DOSEMGKG == 0.5] <- 0.5*28*10^3

# Create data.frame with average values for each timeslot
  dataip.av <- ddply(dataip, .(DOSEMGKG, TADNOM), function(x) {
    data.frame(
      "DOSEMG" = mean(x$DOSEMG, na.rm = T),
      "AMT" = mean(x$AMT, na.rm = T),
      "WT" = mean(x$WT, na.rm = T),
      "TIME" = mean(x$TIME, na.rm = T),
      "DV" = mean(x$DV, na.rm = T)
    )
  })
