# Global script for application for choosing standard curves
# -----------------------------------------------------------------------------
# Set directory strings
  data_dir <- paste(dirname(dirname(getwd())), "raw_data", sep = "/")

# Load libraries
  library(readxl)
  library(stringr)
  library(shiny)
  library(ggplot2)
  library(plyr)

# Customise ggplot2 theme
  theme <- theme_set(theme_bw(base_size = 18))
  
# Load in data
  filename_in1 <- paste0(data_dir, "/mouse plasma protein binding results_9282017_std passing_Short.xls")
  mouse_plas <- read_excel(filename_in1, sheet = "Lenalidomide", skip = 4)
  mouse_plas$Spreadsheet <- 1

  filename_in2 <- paste0(data_dir, "/Plasma protein binding_plasma_human_rerun_finalresults_9282017_Short.xls")
  human_plas <- read_excel(filename_in2, sheet = "Lenalidomide", skip = 4)
  human_plas$Spreadsheet <- 2

  filename_in3 <- paste0(data_dir, "/Plasmaproteinbinding_PBS_results_10022017_Short_171006133211.xls")
  both_pbs <- read_excel(filename_in3, sheet = "Lenalidomide", skip = 4)
  both_pbs$Spreadsheet <- 3
  
  alldata <- rbind(mouse_plas, human_plas, both_pbs)
  
# Fix and subset data
  name_correction <- c("Specified.Amount", "Calculated.Amount")
  names(alldata) <- str_replace_all(names(alldata),"[ ()#]",".")
  names(alldata)[names(alldata) == "Amount"] <- name_correction
  subdata <- alldata[!is.na(alldata$Sample.Type) & 
    (alldata$Sample.Type == "Std Bracket Sample" | alldata$Sample.Type == "QC Sample") &
    alldata$Filename != "zero",
  ]
  subdata$Level <- as.numeric(subdata$Level)
  subdata$Area.Ratio <- as.numeric(subdata$Area.Ratio)
  subdata$Specified.Amount <- as.numeric(subdata$Specified.Amount)
  
  subdata$batch <- 1
  subdata$batch[str_detect(subdata$Filename, "R")] <- 2
  subdata$batch[subdata$Filename == "QC Sample"] <- 0
  
# Create a concentration vector for reference
  concVec <- sort(unique(
    subdata$Specified.Amount[subdata$Sample.Type == "Std Bracket Sample"]
  ))
  
# Create custom row column ui chunk
  fluidRowCol <- function(x) {
    fluidRow( 
      column(6, style = "padding-top: 7px",
        radioButtons(paste0("set", x),
          NULL,
          choices = list(1, 2),
          inline = TRUE
        )  # radioButtons
      ),  # column
      column(6,
        checkboxInput(paste0("status", x),
          NULL,
          value = TRUE
        )  # checkboxInput
      )  # column
    )  # fluidRow
  }