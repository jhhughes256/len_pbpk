# Data Check Script for IV Data
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
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_load.R", sep = "/"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explore Data for Cleaning and Data Extraction
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# View data names and structure
  names(rawiv)
  sort(names(rawiv))
  str(rawiv)

# Check distribution of doses and weights
  with(rawiv, table(Dose..mg.kg., useNA = "always"))
  with(rawiv, table(Mouse.Wt..g., useNA = "always"))

# Check NA's for non-binned data
  any(is.na(rawiv$Sample.ID))
  any(is.na(rawiv$Time..min.))

# Isolate NA's
  rawiv[which(is.na(rawiv$Dose..mg.kg.)), ]
  rawiv[which(is.na(rawiv$Mouse.Wt..g.)), ]
  rawiv[which(is.na(rawiv$Time..min.)), ]
  # NA's appear to be as a result of a repeat sample, it was ascertained from
  # Mitch that these are in fact for the same mouse. Should use the values from
  # each entry to inform the other to form one row.

# Number of samples (includes repeats)
  length(unique(rawiv$Sample.ID))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Clean & Extract Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Convert data for all data
  dataiv <- data.frame("ID" = rawiv$Sample.ID, "DOSEMGKG" = rawiv$Dose..mg.kg.)
  dataiv$DOSEMG <- dataiv$DOSEMGKG*rawiv$Mouse.Wt..g./1000
  dataiv$AMT <- dataiv$DOSEMG*10^6	#dose in ng
  dataiv$WT <- rawiv$Mouse.Wt..g.
  dataiv$TIME <- rawiv$Time..min.
  dataiv$PLA <- rawiv$Plasma.DV..ng.mL.
  dataiv$BRA <- rawiv$Brain..ng.mL.
  dataiv$LVR <- rawiv$Liver..ng.mL.
  dataiv$MSC <- rawiv$Mscl..ng.mL.
  dataiv$HRT <- rawiv$Hrt..ng.mL.
  dataiv$SPL <- rawiv$Spln..ng.mL.
  dataiv$LUN <- rawiv$Lung..ng.mL.
  dataiv$KID <- rawiv$Kidney..ng.mL.

# Clean sample IDs
  IDiv <- end.splitter(dataiv$ID)
  dataiv <- cbind(IDiv, dataiv[,-1])
  # Remove the ID column and bind the new ID columns in its place

# Check the TADNOM against DOSEMGKG
  with(dataiv, table(DOSEMGKG, TADNOM))
  dataiv[which(dataiv$TADNOM == "25"), ]
  # There are one set of values stated to have TADNOM as 25, but values were
  # measured at 20mins, one at 13 which is odd
  # This TADNOM is likely meant to be 20
  dataiv$TADNOM[which(dataiv$TADNOM == "25")] <- "20"

# Fix repeated UID
# The original sample has NA WT
# The repeated sample has NA DOSEMGKG
  subiv <- dataiv[which(with(dataiv, TADNOM == "20" & DOSEMGKG == 5 & !is.na(dataiv$WT))), ]
  IDori <- dataiv[which(is.na(dataiv$WT)), ]
  IDrep <- arrange(dataiv[which(is.na(dataiv$DOSEMGKG)), ], UID)
  tissue.mean <- function(x) {
    colwise(mean)(x[c("PLA", "BRA", "LVR", "MSC", "HRT", "SPL", "LUN", "KID")])
  }
  rbind(tissue.mean(subiv), tissue.mean(IDori), tissue.mean(IDrep))
  # it looks like the repeat row is no good, and was done due to the plasma
  # concs not working in the original test.

# Combine into one row, use DV and DOSEMGKG from orig, use WT & TIME from repeat
  IDori$TIME <- IDrep$TIME
  IDori$WT <- IDrep$WT
  IDori$DOSEMG <- IDori$DOSEMGKG*IDori$WT/1000
  IDori$AMT <- IDori$DOSEMG*10^6  # units: ng
  IDori$PLA <- NA

# Then replace both samples with this combined sample and order the data.frame
  IDrem1 <- which(is.na(dataiv$WT))
  IDrem2 <- which(is.na(dataiv$DOSEMGKG))
  dataiv <- arrange(
    rbind(
      dataiv[-c(IDrem1, IDrem2), ], IDori
    ),  # rbind
    DOSEMGKG, TIME
  )  # arrange

# Make TADNOM numeric
  dataiv$TADNOM <- as.numeric(dataiv$TADNOM)

# Create data.frame with average values for each timeslot
  dataiv.av <- ddply(dataiv, .(DOSEMGKG, TADNOM), function(x) {
    data.frame(
      "DOSEMG" = mean(x$DOSEMG, na.rm = T),
      "AMT" = mean(x$AMT, na.rm = T),
      "WT" = mean(x$WT, na.rm = T),
      "TIME" = mean(x$TIME, na.rm = T),
      "PLA" = mean(x$PLA, na.rm = T),
      "BRA" = mean(x$BRA, na.rm = T),
      "LVR" = mean(x$LVR, na.rm = T),
      "MSC" = mean(x$MSC, na.rm = T),
      "HRT" = mean(x$HRT, na.rm = T),
      "SPL" = mean(x$SPL, na.rm = T),
      "LUN" = mean(x$LUN, na.rm = T),
      "KID" = mean(x$KID, na.rm = T)
    )
  })
  
  dataiv.sd <- ddply(dataiv, .(DOSEMGKG, TADNOM), function(x) {
    data.frame(
      "DOSEMG" = sd(x$DOSEMG, na.rm = T),
      "AMT" = sd(x$AMT, na.rm = T),
      "WT" = sd(x$WT, na.rm = T),
      "TIME" = sd(x$TIME, na.rm = T),
      "PLA" = sd(x$PLA, na.rm = T),
      "BRA" = sd(x$BRA, na.rm = T),
      "LVR" = sd(x$LVR, na.rm = T),
      "MSC" = sd(x$MSC, na.rm = T),
      "HRT" = sd(x$HRT, na.rm = T),
      "SPL" = sd(x$SPL, na.rm = T),
      "LUN" = sd(x$LUN, na.rm = T),
      "KID" = sd(x$KID, na.rm = T)
    )
  })
