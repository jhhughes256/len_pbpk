# Using PKSim simulation data and plotting in the R environment
# -----------------------------------------------------------------------------
# Done in RStudio using projects, so getwd() returns git.dir
  rm(list = ls(all = T))
  reponame <- "len_pbpk"
  git.dir <- substr(getwd(), 1, nchar(getwd()) - nchar(reponame) - 1)

  
# Load libraries
  library(ggplot2)

# Source data
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))
  ref_data <- read.csv("raw_data/PKSim_data/AllTis.Ref-Results.csv", 
    sep = ";", quote = "")
  base_data <- read.csv("raw_data/PKSim_data/AllTis.Base-Results.csv", 
    sep = ";", quote = "")
  obs_data <- dataiv.av
  
# Format data
  col_names <- c("id", "time", "peripheral.venous.blood", "PLA", 
    "brain", "heart", "kidney", "lumen", "liver", "lung", "muscle", "spleen")
  names(ref_data) <- col_names
  names(base_data) <- col_names
  
# -----------------------------------------------------------------------------
  p1 <- 