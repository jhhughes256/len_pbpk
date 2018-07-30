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
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Simulation data is output from MoBi
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
  
# Read in simulation data and rename columns
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "PKSimArterialAve.xlsx"
  data_raw <- read_excel(paste0(file.dir, model.name))
  data <- as.data.frame(data_raw)
  names(data) <- c("TIME", "Plasma", "Blood Cells")
  
# Convert time column to minutes
  data$TIME <- data$TIME*60
  
# Create forcing function and 
  forcingFunction <- function(x) {
    exp(-0.003954853*x + -5.604366331) + 
    exp(-0.02089863*x + -0.28593179) + 
    exp(-0.140071*x + 1.851636) +
    exp(-4.15563*x + 4.388731)
  }

# Make plotdata
  plotdata <- data
  plotdata$PRED <- forcingFunction(plotdata$TIME)
  plotdata$RES <- with(plotdata, PRED - Plasma)
  plotdata$PROPRES <- with(plotdata, RES/PRED*100)
  
# Make plots
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_point(aes(x = TIME, y = PROPRES), data = plotdata[-1,])
  p1 <- p1 + geom_hline(yintercept = 0, linetype = "dashed", colour = "green4")
  p1 <- p1 + xlab("Time (min)")
  p1 <- p1 + scale_y_continuous("Proportional Error (%)", labels = comma)
  p1
  
  