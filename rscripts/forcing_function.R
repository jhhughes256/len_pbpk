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
    # git.dir <- "E:/Hughes/Git"
    # git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "len_pbpk"
  }
  # Load libraries
  library(reshape2)
  library(GA)
  library(plyr)

  # Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# Grab plasma data
# Convert data to units PKSim is expecting
# Units are ng/mL aka ug/L
# convert to umol/L as when multiplied by volume in the arterial blood it gives
# umol which is PKSim's native unit
  dataiv.pla <- with(dataiv[!is.na(dataiv$PLA),], data.frame(
    TADNOM = TADNOM,
    DOSEMGKG = DOSEMGKG,
    TIME = TIME,
    DV = PLA/259.26,
    DVNORM = PLA/(DOSEMG*259.26)
  ))

  dataiv.av <- ddply(dataiv.pla, .(TADNOM, DOSEMGKG), function(x) {
    with(x, data.frame(
      TIME = mean(TIME, na.rm = T),
      DV = mean(DV, na.rm = T),
      DVNORM = mean(DVNORM, na.rm = T)
    ))
  })
  
# Define linear interpolation function
# This can easily be used with ddply
  interpLinear <- function(x, y) {
    m <- diff(y)/diff(x)
    b <- head(y, -1) - head(x, -1)*m
    data.frame(
      TADNOM = x,
      DV = y,
      m = c(m, tail(m, 1)),
      c = c(b, tail(b, 1))
    )
  }
  
  ffdata <- ddply(dataiv.av, .(DOSEMGKG), function(x) {
    interpLinear(c(0, x$TADNOM), c(0, x$DV))
  })
  
# Check to see if this is working as expected
# First set up linear interpolation data
  ffpred <- data.frame(
    DOSEMGKG = rep(c(0.5, 1.5, 5, 10), each = 601),
    TIME = rep(seq(0, 600, by = 1), times = 4)
  )
  roundDownToValue <- function(x, value) {
    #debug
    # x <- ffpred$TIME
    # value <- unique(ffdata$TIME)
    out <- double(length(x))
    for (i in 1:length(value)) {
      out[x >= value[i]] <- value[i]
    }
    out
  }
  
  ffpred$TADNOM <- roundDownToValue(ffpred$TIME, unique(ffdata$TADNOM))
  plotpre <- merge(ffpred, ffdata, by = c("DOSEMGKG", "TADNOM"))
  plotpre$DV <- with(plotpre, TIME*m + c)
    
# Set up factors for plotting
  plotobs <- dataiv.av
  plotobs$DOSEMGKGf <- factor(plotobs$DOSEMGKG)
  plotpre$DOSEMGKGf <- factor(plotpre$DOSEMGKG)
  
# Then plot
  p <- NULL
  p <- ggplot()
  p <- p + geom_point(aes(x = TADNOM, y = DV), data = plotobs)
  p <- p + geom_line(aes(x = TIME, y = DV), data = plotpre)
  p <- p + scale_y_log10()
  p <- p + facet_wrap(~DOSEMGKGf)
  p 
  