# Calculation of protein binding from equilibrium dialysis results
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
      reponame <- "len_pbpk"
      data_folder <- paste(git.dir, reponame, "raw_data", sep = "/")
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
      data_folder <- "E:/Hughes/Data/RAW_NonClinical/protein_binding/20171009"
    }
    rm("wd")
  }

# Load libraries
  library(readxl)
  library(stringr)
  library(plyr)

# -----------------------------------------------------------------------------
# Load in data

  filename_in <- paste0(data_folder, "/mouse plasma protein binding results_9282017_std passing_Short.xls")
  mouse_plas <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)
  mouse_plas$Spreadsheet <- "mouse"

  filename_in <- paste0(data_folder, "/Plasma protein binding_plasma_human_rerun_finalresults_9282017_Short.xls")
  human_plas <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)
  human_plas$Spreadsheet <- "human"

  filename_in <- paste0(data_folder, "/Plasmaproteinbinding_PBS_results_10022017_Short_171006133211.xls")
  both_pbs <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)
  both_pbs$Spreadsheet <- "pbs"

# Combine data
  any(names(mouse_plas) != names(human_plas) & names(mouse_plas) != names(both_pbs))
  alldata <- rbind(mouse_plas, human_plas, both_pbs)

# Fix column names
  name_correction <- c("Specified.Amount", "Calculated.Amount")
  names(alldata) <- str_replace_all(names(alldata),"[ ()#]",".")
  names(alldata)[names(alldata) == "Amount"] <- name_correction

# Take subset of sample data and add covariate info
  subdata <- alldata[(alldata$Sample.Type == "Std Bracket Sample" | alldata$Sample.Type == "QC Sample") & !is.na(alldata$Sample.Type) & alldata$Filename != "zero",]
  subdata$Level <- as.numeric(subdata$Level)

# Linear regression
  with(subdata[subdata$Spreadsheet == "mouse",], plot(Area.Ratio ~ Level))
  with(subdata[subdata$Spreadsheet == "mouse",], 1/Level**2)

  lm(Area.Ratio ~ Level, data = subdata, subset = Spreadsheet == "mouse", weight = 1/Level**2)
