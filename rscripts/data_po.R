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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explore Data for Cleaning and Data Extraction
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# View data names and structure
  names(rawpo)
  sort(names(rawpo))
  str(rawpo)

# Check distribution of doses and weights
  with(rawpo, table(Dose..mg.kg., useNA = "always"))
  with(rawpo, table(Mouse.Wt..g., useNA = "always"))
  # no weight data!
  # will have to impute, mean or median of all IV mice

# Check NA's for non-binned data
  any(is.na(rawpo$Time..min.))
  with(rawpo, table(Sample.ID, useNA = "always"))[1]
  # 5 cases of Sample.ID of ""

# Isolate ""'s in Sample.ID
  rawpo[rawpo$Sample.ID == "", ]

# Number of samples (no repeat samples)
  dim(rawpo)[1]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Clean & Extract Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Convert data for all data
  datapo <- data.frame("ID" = rawpo$Sample.ID, "DOSEMGKG" = rawpo$Dose..mg.kg.)
  datapo$DOSEMG <- datapo$DOSEMGKG*28/1000
  datapo$AMT <- datapo$DOSEMG*10^6  # ng
  datapo$WT <- 28  # median weight of IV bolus mice (g)
  datapo$TIME <- rawpo$Time..min.
  datapo$DV <- rawpo$Plasma.DV..ng.mL.

# Clean sample IDs
  IDpo <- end.splitter(datapo$ID)
  datapo <- cbind(IDpo, datapo[, -1])

# Create data.frame with average values for each timeslot
  datapo.av <- ddply(datapo, .(DOSEMGKG, TADNOM), function(x) {
    data.frame(
      "DOSEMG" = mean(x$DOSEMG, na.rm = T),
      "AMT" = mean(x$AMT, na.rm = T),
      "WT" = mean(x$WT, na.rm = T),
      "TIME" = mean(x$TIME, na.rm = T),
      "DV" = mean(x$DV, na.rm = T)
    )
  })
