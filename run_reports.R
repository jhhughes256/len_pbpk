# If you wish to run this you will require one of two types of installations
#   - R + pandoc
#   - R + RStudio
# This will allow the .Rmd files to render using pandoc

# Set working directory
#  work.dir <- "C:/Users/hugjh001/Documents/len_pbpk"
  work.dir <- "E:/Hughes/Git/len_pbpk"
  setwd(work.dir)

# Functions required if initialising directory
#  chooseCRANmirror()
#  packrat::init()

# Source packages from packrat and bespoke functions
  source("packrat/init.R")
  source("functions/make_report.R")

# Run reports
