# Data Loading Script for Lenalidomide Mouse Data
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
    #git.dir <- "E:/Hughes/Git"
    git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    reponame <- "len_pbpk"
  }

# Set the working directory
  master.dir <- paste(git.dir, reponame, sep = "/")
  setwd(master.dir)

# Load libraries
  library(readxl)
  library(plyr)
  library(reshape2)
  library(ggplot2)
  library(grid)
  library(stringr)

# Source functions to be used with the script
  source(paste(git.dir, reponame, "functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "functions",
    "endsplitter.R", sep = "/"))
  source(paste(git.dir, reponame, "functions",
    "datacheck_plotfn.R", sep = "/"))

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Organise working and output directories
  plot.dir <- paste(master.dir, "plot", sep = "/")
  data.dir <- paste(master.dir, "produced_data", sep = "/")

# -----------------------------------------------------------------------------
# Load in the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Using read_excel due to .xls format
  file.name.in <- "raw_data/All_Tissue_PK_Data_Summary.xls"
  dataraw <- suppressWarnings(read_excel(file.name.in,
    col_types = c(rep("text", 2), rep("numeric", 19)),
    sheet = "Data"
  ))  # read_excel

# Remove problematic characters from column names
  colnames(dataraw) <- gsub.all(
    c(" ", "(", ")", "/"),
    rep(".", 4),
    colnames(dataraw)
  )

# General cleaning of data
  datanew <- dataraw[!(is.na(dataraw$Dose..mg.kg.) & is.na(dataraw$Sample.ID)), ]
  datanew$Dose..mg.kg. <- as.numeric(datanew$Dose..mg.kg.)
  datanew[is.na(str_detect(datanew$Sample.ID, "_")), 1] <- ""
  # Final bit of cleaning here is for endsplitter to work

# Split the data into 3 datasets based on the titles in the spreadsheet
# +1 & -1 used to avoid using the rows that split up the data
# po and pr data have no tissue data, only first 7 columns needed
  rowsplit <- c(
    which(str_detect(datanew$Sample.ID, "Data")),
    length(datanew$Sample.ID)
  )

  rawiv <- data.frame(datanew[c((rowsplit[1] + 1):(rowsplit[2] - 1)), ])
  rawiv.av <- data.frame(rawiv[!is.na(rawiv$Avg.time..min.),])
  rawpo <- data.frame(datanew[c((rowsplit[2] + 1):(rowsplit[3] - 1)), ][1:7])
  rawpo.av <- data.frame(rawpo[!is.na(rawpo$Plasma.Avg..ng.mL.), ])
  rawip <- data.frame(datanew[c((rowsplit[3] + 1):(rowsplit[4])), ][1:7])
  rawip.av <- data.frame(rawip[!is.na(rawip$Plasma.Avg..ng.mL.), ])
