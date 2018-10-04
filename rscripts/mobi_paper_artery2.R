# -----------------------------------------------------------------------------
# Ready workspace
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
  library(scales)

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Simulation Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation data is output from MoBi
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
  
# Tidying script consists of reading in arterial data for all five single-tissue
#   models for both mouse and humans (10 simulations total). These are then 
#   collated into one dataset with columns for TISSUE and SPECIES. Plots will 
#   use TISSUE for colour and SPECIES for facetting.
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Mouse Data
  
# Define simulation file names
  file.dir <- "raw_data/MoBi_paper/"
  mouse.model.art <- "PKSimArterial_Base.xlsx"
  mouse.model.bra <- "PKSimBrain_Base.xlsx"
  mouse.model.hrt <- "PKSimHeart_Base.xlsx"
  mouse.model.kid <- "PKSimKidney_Base.xlsx"
  mouse.model.lng <- "PKSimLung_Base.xlsx"
  
# Define forcing function
  mouseFF <- function(x) {
    exp(-0.004047586*x + -5.591089002) + 
    exp(-0.02130095*x + -0.29374929) + 
    exp(-0.170071*x + 1.851636) +
    exp(-4.15563*x + 3.528731)
  }
  
# Read in simulation data and rename/subset columns
  raw.art <- read_excel(paste0(file.dir, mouse.model.art))  # read
  artdata <- as.data.frame(raw.art)
  names(artdata) <- c("TIME", "affPL", "affBC")  # rename
  
  raw.bra <- read_excel(paste0(file.dir, mouse.model.bra))  # read
  bradata <- as.data.frame(raw.bra)
  names(bradata) <- c("TIME", "effPL", "effBC", "affPL", "affBC")  # rename
  bradata <- bradata[, c("TIME", "affPL", "affBC")]  # subset
  
  raw.hrt <- read_excel(paste0(file.dir, mouse.model.hrt))  # read
  hrtdata <- as.data.frame(raw.hrt)
  names(hrtdata) <- c("TIME", "effPL", "effBC", "affPL", "affBC")  # rename
  hrtdata <- hrtdata[, c("TIME", "affPL", "affBC")]  # subset
  
  raw.kid <- read_excel(paste0(file.dir, mouse.model.kid))  # read
  kiddata <- as.data.frame(raw.kid)
  names(kiddata) <- c("TIME", "effPL", "effBC", "affPL", "affBC")  # rename
  kiddata <- kiddata[, c("TIME", "affPL", "affBC")]  # subset
  
# For the lung the venous concentrations are named as afferent concentrations
  raw.lng <- read_excel(paste0(file.dir, mouse.model.lng))  # read
  lngdata <- as.data.frame(raw.lng)
  names(lngdata) <- c("TIME", "affPL", "affBC", "effPL", "effBC")  # rename
  lngdata <- lngdata[, c("TIME", "affPL", "affBC")]  # subset
  
# Define tissue column and bind together before defining species column
  artdata$TISSUE <- "Artery"
  bradata$TISSUE <- "Brain"
  hrtdata$TISSUE <- "Heart"
  kiddata$TISSUE <- "Kidney"
  lngdata$TISSUE <- "Lung"
  mdata <- rbind(artdata, bradata, hrtdata, kiddata, lngdata)
  mdata$SPECIES <- "Mouse"
  
# Convert time to minutes and define forcing function concentrations
  mdata$TIME <- mdata$TIME*60
  mdata$PRED <- mouseFF(mdata$TIME)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Human Data
  
# Define simulation file names
  file.dir <- "raw_data/MoBi_paper/"
  human.model.art <- "PKSimArterial_BaseHuman.xlsx"
  human.model.bra <- "PKSimBrain_BaseHuman.xlsx"
  human.model.hrt <- "PKSimHeart_BaseHuman.xlsx"
  human.model.kid <- "PKSimKidney_BaseHuman.xlsx"
  human.model.lng <- "PKSimLung_BaseHuman.xlsx"
  
# Define forcing function
  humanFF <- function(x) {
    exp(-0.01466504*x + 0.70681102) + 
    exp(-0.06788317*x + 2.23756615) + 
    exp(-2.35813742*x + 4.89291408)
  }
  
# Read in simulation data and rename/subset columns
  raw.art <- read_excel(paste0(file.dir, human.model.art))  # read
  artdata <- as.data.frame(raw.art)
  names(artdata) <- c("TIME", "affPL", "affBC")  # rename
  
  raw.bra <- read_excel(paste0(file.dir, human.model.bra))  # read
  bradata <- as.data.frame(raw.bra)
  names(bradata) <- c("TIME", "effPL", "effBC", "affPL", "affBC")  # rename
  bradata <- bradata[, c("TIME", "affPL", "affBC")]  # subset
  
  raw.hrt <- read_excel(paste0(file.dir, human.model.hrt))  # read
  hrtdata <- as.data.frame(raw.hrt)
  names(hrtdata) <- c("TIME", "effPL", "effBC", "affPL", "affBC")  # rename
  hrtdata <- hrtdata[, c("TIME", "affPL", "affBC")]  # subset
  
  raw.kid <- read_excel(paste0(file.dir, human.model.kid))  # read
  kiddata <- as.data.frame(raw.kid)
  names(kiddata) <- c("TIME", "effPL", "effBC", "affPL", "affBC")  # rename
  kiddata <- kiddata[, c("TIME", "affPL", "affBC")]  # subset
  
  raw.lng <- read_excel(paste0(file.dir, human.model.lng))  # read
  lngdata <- as.data.frame(raw.lng)
  names(lngdata) <- c("TIME", "affPL", "affBC", "effPL", "effBC")  # rename
  lngdata <- lngdata[, c("TIME", "affPL", "affBC")]  # subset
  
# Define tissue column and bind together before defining species column
  artdata$TISSUE <- "Artery"
  bradata$TISSUE <- "Brain"
  hrtdata$TISSUE <- "Heart"
  kiddata$TISSUE <- "Kidney"
  lngdata$TISSUE <- "Lung"
  hdata <- rbind(artdata, bradata, hrtdata, kiddata, lngdata)
  hdata$SPECIES <- "Human"
  
# Convert time to minutes and define forcing function concentrations
  hdata$TIME <- hdata$TIME*60
  hdata$PRED <- humanFF(hdata$TIME)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Collated Data
  data <- rbind(mdata, hdata)
  
# Make plotdata
  plotdata <- data[data$TIME > 0, ]
  plotdata$RES <- with(plotdata, affPL - PRED)
  plotdata$PROPRES <- with(plotdata, RES/PRED*100)
  
# Define colourblind-friendly palette
  cbPalette <- c("#D55E00", "#E69F00", "#CC79A7", "#009E73", "#0072B2")
  
# Make plots
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = TIME/60, y = PROPRES, colour = TISSUE), 
    data = plotdata, size = 1.5)
  p1 <- p1 + geom_hline(yintercept = 0, linetype = "dashed")
  p1 <- p1 + scale_colour_manual("Model", values = cbPalette)
  p1 <- p1 + scale_x_continuous("Time (hour)", breaks = c(0, 4, 8, 12, 16))
  p1 <- p1 + scale_y_continuous("Proportional Error (%)", labels = comma)
  p1 <- p1 + facet_wrap(~SPECIES)
  p1
  
# Produce Figure 2 for the MoBi methods paper
  ggsave("produced_data/mobi_paper_artery2.png", width = 17.4, height = 10.4, 
    units = c("cm"))
  # ggsave("produced_data/Figure2.eps", width = 17.4, height = 23.4, 
  #   units = c("cm"), dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  