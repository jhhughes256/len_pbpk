# global.R for len_pbpk model exploration
# -----------------------------------------------------------------------------
# Load package libraries
  library(shinydashboard)  # Dashboard ui
  library(reshape2)  # Conversion of data format from wide to long
  library(ggplot2)	# Plotting
  library(mrgsolve)	 # Metrum differential equation solver for pharmacometrics
  library(markdown)
# Define a custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 20))
# Set colour palette
  cPalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF")
# Set times for simulation
  TIME <- seq(from = 0, to = 100, by = 0.2)
# Source the models
  source("model.R")
  NRmelt <- read.csv("www/data.csv")
