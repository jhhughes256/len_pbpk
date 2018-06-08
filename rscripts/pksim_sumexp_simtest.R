# PKSim/MoBi Forcing Function Potential Fix? - Testing Script
#
# Fitting sum of exponential equation to PKSim simulation data for use in a
# MoBi simulation that will have it's results checked back against the PKSim
# simulation data for tissues to see if all of this mess has been sorted out!
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32", "C:/Users/hugjh001/Documents/len_pbpk")

    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2] | getwd() == wd[5]) {
      git.dir <- "C:/Users/hugjh001/Documents"
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Load libraries
  library(readxl)
  library(stringr)
  library(plyr)
  library(reshape2)
  library(ggplot2)

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit sum of exponentials to PKSim arterial data
# Source sum of exponentials functions
  library(GA)
  source(paste(git.dir, "optinterval", "paper_functions.R", sep = "/"))
  
# Read in PKSim data
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "PKSim_withArterial.xlsx"
  pksim.raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Rename columns
  pksim <- pksim.raw
  names(pksim) <- c(
    "TIME", 
    paste0(
      rep(c("pls", "bc"), times = 2), 
      rep(c("Venous", "Arterial"), each = 2)
    ),
    paste0(
      rep(c("pls", "bc", "is", "ic", "tis"), times = 4), 
      rep(c("Brain", "Heart", "Kidney", "Lung"), each = 5)
    )
  )
  
# Convert time to minutes
  pksim$TIME <- pksim$TIME*60
  
# Prepare data.frame for using optim.sumexp
  fit.pksim <- data.frame(
    Time = pksim$TIME,
    DV = pksim$plsArterial
  )
  # fit.pksim <- data.frame(
  #   Time = pksim$TIME,
  #   DV = pksim$plsVenous
  # )
  
# Fit sum of exponentials
  sumexp.out <- optim.sumexp(fit.pksim[-c(1, fit.pksim$Time > 25),], nexp = 2)
  sumexp.mod <- best.sumexp.aic(sumexp.out)
  sumexp.mod
  
# Visual goodness of fit plot
  # debug
  # venous
  # best <- c(-0.003988262, -0.02121703,  -0.08636453, -0.9093104, -5.594721923, -0.25598858, 0.89099724, 2.6222481)
  # arterial
  best <- c(-0.003954853, -0.02089863, -0.140071, -4.155630, -5.604366331, -0.28593179, 1.851636, 4.388731)
  # fit.pksim$PRED <- pred.sumexp(sumexp.mod$sumexp, fit.pksim$Time)
  fit.pksim$PRED <- pred.sumexp(best, fit.pksim$Time)
  
  p <- NULL
  p <- ggplot(data = fit.pksim)
  p <- p + geom_point(aes(x = Time, y = DV), colour = "red")
  p <- p + geom_line(aes(x = Time, y = PRED), linetype = "dashed")
  p <- p + scale_y_log10()  # scale used as zero is used
  # p <- p + coord_cartesian(xlim = c(0, 50))
  p
  
  p <- NULL
  p <- ggplot(data = fit.pksim)
  p <- p + geom_point(aes(x = Time, y = (PRED - DV)/PRED))
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + coord_cartesian(ylim = c(-0.5, 1))
  p
  