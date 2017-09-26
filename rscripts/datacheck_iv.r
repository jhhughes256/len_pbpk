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
    rm(list=ls(all=TRUE))
    graphics.off()
    git.dir <- "E:/Hughes/Git"
    #git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    reponame <- "len_pbpk"
  }

# Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_load.R", sep = "/"))

  scriptname <- "datacheck_iv"
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

  dataiv.med <- ddply(dataiv, .(DOSEMGKG, TADNOM), function(x) {
    data.frame(
      "DOSEMG" = median(x$DOSEMG, na.rm = T),
      "AMT" = median(x$AMT, na.rm = T),
      "WT" = median(x$WT, na.rm = T),
      "TIME" = median(x$TIME, na.rm = T),
      "PLA" = median(x$PLA, na.rm = T),
      "BRA" = median(x$BRA, na.rm = T),
      "LVR" = median(x$LVR, na.rm = T),
      "MSC" = median(x$MSC, na.rm = T),
      "HRT" = median(x$HRT, na.rm = T),
      "SPL" = median(x$SPL, na.rm = T),
      "LUN" = median(x$LUN, na.rm = T),
      "KID" = median(x$KID, na.rm = T)
    )
  })

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data Check
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Check subject numbers
  with(dataiv, table(ID))
  with(dataiv, table(UID))
  any(with(dataiv, table(UID)) > 2)
  # Successfully removed repeat sample

# Check PK dose data
  with(dataiv, table(DOSEMGKG, TADNOM))

# Check the dose columns
  with(dataiv, table(DOSEMG))
  hist(dataiv$DOSEMG)

# Check distribution of DV
  meltiv <- melt(dataiv, id = c("UID", "ID", "TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME"))
  colnames(meltiv) <- c(head(colnames(meltiv), 8), "TISSUE", "DV")
  iv.distplot(meltiv, "alldata", plot.out)

  meltiv.av <- melt(dataiv.av, id = c("DOSEMGKG","TADNOM", "DOSEMG", "AMT", "WT", "TIME"))
  colnames(meltiv.av) <- c(head(colnames(meltiv.av), 6), "TISSUE", "DV")
  iv.distplot(meltiv.av, "meandata", plot.out)

  meltiv.med <- melt(dataiv.med, id = c("DOSEMGKG","TADNOM", "DOSEMG", "AMT", "WT", "TIME"))
  colnames(meltiv.med) <- c(head(colnames(meltiv.med), 6), "TISSUE", "DV")
  iv.distplot(meltiv.med, "mediandata", plot.out)

# Calculate dose normalized concentrations and mark missing DV
# Units are ng/ml per mg
  meltiv$DVNORM <- meltiv$DV/meltiv$DOSEMG
  meltiv$MDV <- ifelse(is.na(meltiv$DV), 1, 0)
  meltiv.av$DVNORM <- meltiv.av$DV/meltiv.av$DOSEMG
  meltiv.av$MDV <- ifelse(is.na(meltiv.av$DV), 1, 0)

# Check maximum concentrations per dose and per tissue
  ddply(meltiv.av, .(TISSUE, DOSEMGKG), function(x) {
    max(x$DV, na.rm = T)
  })

  ddply(meltiv.av, .(TISSUE, DOSEMGKG), function(x) {
    max(x$DVNORM, na.rm = T)/1000
  })

# Plot PK data
  iv.CvTplot(meltiv[meltiv$DOSEMGKG == 0.5,], "alldata", plot.out)
  iv.CvTplot(meltiv[meltiv$DOSEMGKG == 1.5,], "alldata", plot.out)
  iv.CvTplot(meltiv[meltiv$DOSEMGKG == 5,], "alldata", plot.out)
  iv.CvTplot(meltiv[meltiv$DOSEMGKG == 10,], "alldata", plot.out)
  iv.CvTplot(meltiv.av[meltiv.av$DOSEMGKG == 0.5,], "meandata", plot.out)
  iv.CvTplot(meltiv.av[meltiv.av$DOSEMGKG == 1.5,], "meandata", plot.out)
  iv.CvTplot(meltiv.av[meltiv.av$DOSEMGKG == 5,], "meandata", plot.out)
  iv.CvTplot(meltiv.av[meltiv.av$DOSEMGKG == 10,], "meandata", plot.out)

  iv.CvTplot(meltiv, "alldata", plot.out, dosenorm = T)
  iv.CvTplot(meltiv.av, "meandata", plot.out, dosenorm = T)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Covariate Data Check
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Count missing by dose group
  missingbydose <- ddply(dataiv, .(DOSEMGKG), colwise(calculate.percent.missing))

  filename.out <- paste(data.out, "Missing_by_dosegroup.csv", sep = "/")
  write.csv(missingbydose, file = filename.out, row.names = F)

# Missing by TADNOM
  missingbytad <- ddply(dataiv, .(TADNOM), colwise(calculate.percent.missing))

  filename.out <- paste(output.dir, "Missing_by_nominaltime.csv", sep = "/")
  write.csv(missingbytad, file = filename.out, row.names = F)

# Missing summary
  missingsummary <- ddply(dataiv, .(DOSEMGKG, TADNOM), colwise(calculate.percent.missing))

  filename.out <- paste(output.dir, "Missing_summary.csv", sep = "/")
  write.csv(missingbytad, file = filename.out, row.names = F)

# DV count
  DVsum <- function(data, group, uid, col) {
    ind <- unique(data[, group])
    out <- data.frame(
      matrix(nrow = length(ind), ncol = length(col)*2 + 1)
    )
    out[1] <- ind
    for (i in 1:length(ind)) {
      sub <- data[data[group] == ind[i], ]
      out[i, 2] <- length(unique(sub[, uid]))
      for (j in 1:length(col)) {
        out[i, j*2 + 1] <- length(which(!is.na(sub[, col[j]])))
        out[i, j*2 + 2] <- out[i, j*2 + 1]/out[i, 2]
      }
    }
    names(out) <- c(
      group, "SUBcount", paste0(rep(col, each = 2), c("count", "perSUB"))
    )
    return(out)
  }
  DVsummary <- DVsum(dataiv, "DOSEMGKG", "UID",
    c("AMT", "PLA", "BRA", "LVR", "MSC", "HRT", "SPL", "LUN", "KID")
  )

  filename.out <- paste(data.out, "DV_sum.csv",sep="/")
  write.csv(DVsummary, file = filename.out)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Covariate Data Check
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
  plotdata <- meltiv
  BINnumber <- 3

  plotdata$DOSEMGKGf <- as.factor(plotdata$DOSEMGKG)
  plotdata$TISSUEf <- as.factor(plotdata$TISSUE)
  plotdata$WT_bin <- as.factor(ave(plotdata$WT, cut(plotdata$WT, BINnumber), FUN = median))
  plotdata$DOSE_bin <- as.factor(ave(plotdata$DOSEMG, cut(plotdata$DOSEMG, BINnumber), FUN = median))

  plotByFactor("DOSE_bin", "Binned Dose (ug)", plotdata, plot.out)
  plotByFactor("WT_bin", "Binned Weight (g)", plotdata, plot.out)
