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
  
# Read in MoBi data for normal Kidney simulation
  model.name <- "PKSimKidney2.xlsx"
  kidlo.raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Rename columns
  kidlo <- kidlo.raw
  names(kidlo) <- c(
    "TIME", 
    paste0(
      c(rep(c("pls", "bc"), times = 3), "is", "ic", "ur", "tis"),
      c(rep(c("Venous", "Arterial"), each = 2), rep("Kidney", times = 6))
    )
  )
  
# Convert time to minutes
  kidlo$TIME <- kidlo$TIME*60
  kidlo$RES <- kidlo$tisKidney - pksim$tisKidney
  kidlo$PROPRES <- with(kidlo, RES/tisKidney)
  
# Read in MoBi data for fixed Kidney simulation
  model.name <- "PKSimKidneyFix2.xlsx"
  kidhi.raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Rename columns
  kidhi <- kidhi.raw
  names(kidhi) <- c(
    "TIME", 
    paste0(
      c(rep(c("pls", "bc"), times = 3), "is", "ic", "ur", "tis"),
      c(rep(c("Venous", "Arterial"), each = 2), rep("Kidney", times = 6))
    )
  )
  
# Convert time to minutes
  kidhi$TIME <- kidhi$TIME*60
  kidhi$RES <- kidhi$tisKidney - pksim$tisKidney
  kidhi$PROPRES <- with(kidhi, RES/tisKidney)
  
# Read in MoBi data when proportional gain uses half the elim. rate const.
  model.name <- "PKSimKidneyHalf2.xlsx" #"PKSimKidneyEst.xlsx"
  kidmi.raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Rename columns
  kidmi <- kidmi.raw
  names(kidmi) <- c(
    "TIME", 
    paste0(
      c(rep(c("pls", "bc"), times = 3), "is", "ic", "ur", "tis"),
      c(rep(c("Venous", "Arterial"), each = 2), rep("Kidney", times = 6))
    )
  )
  
# Convert time to minutes
  kidmi$TIME <- kidmi$TIME*60
  kidmi$RES <- kidmi$tisKidney - pksim$tisKidney
  kidmi$PROPRES <- with(kidmi, RES/tisKidney)
  
# Plot kidney concentrations against each other
  # p <- NULL
  # p <- ggplot()
  # p <- p + geom_line(aes(x = TIME, y = tisKidney), data = pksim, colour = "red")
  # p <- p + geom_line(aes(x = TIME, y = tisKidney), data = kidlo, linetype = "dashed")
  # p <- p + geom_line(aes(x = TIME, y = tisKidney), data = kidhi, linetype = "dashed")
  # p <- p + geom_line(aes(x = TIME, y = tisKidney), data = kidmi, linetype = "dotted", colour = "blue")
  # p <- p + scale_y_log10()
  # p
  # 
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = TIME, y = PROPRES), data = kidlo, colour = "red")
  p <- p + geom_line(aes(x = TIME, y = PROPRES), data = kidhi, colour = "blue")
  p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "black")
  p <- p + geom_hline(yintercept = mean(kidlo$PROPRES, na.rm = T), 
    linetype = "dashed", colour = "red")
  p <- p + geom_hline(yintercept = mean(kidhi$PROPRES, na.rm = T), 
    linetype = "dashed", colour = "blue")
  p <- p + geom_line(aes(x = TIME, y = PROPRES), data = kidmi, colour = "green4")
  p <- p + geom_hline(yintercept = mean(kidmi$PROPRES, na.rm = T),
    linetype = "dashed", colour = "green4")
  p
  