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
  mouse_plas1 <- read_excel(filename_in1, sheet = "Lenalidomide", skip = 4)
  mouse_plas1$Run <- 1
  mouse_plas1$Spreadsheet <- 1
  
  filename_in2 <- paste0(data_dir, "/mouse plasma 4h_results_10182017_Short.xls")
  mouse_plas2 <- read_excel(filename_in2, sheet = "Lenalidomide", skip = 4)
  mouse_plas2$Run <- 2
  mouse_plas2$Spreadsheet <- 1

  filename_in3 <- paste0(data_dir, "/Plasma protein binding_plasma_human_rerun_finalresults_9282017_Short.xls")
  human_plas1 <- read_excel(filename_in3, sheet = "Lenalidomide", skip = 4)
  human_plas1$Run <- 1
  human_plas1$Spreadsheet <- 2
  
  filename_in4 <- paste0(data_dir, "/human plasma protein binding_4h_LenaAPCI_results_10182017_Short.xls")
  human_plas2 <- read_excel(filename_in4, sheet = "Lenalidomide", skip = 4)
  human_plas2$Run <- 2
  human_plas2$Spreadsheet <- 2

  filename_in5 <- paste0(data_dir, "/Plasmaproteinbinding_PBS_results_10022017_Short_171006133211.xls")
  both_pbs1 <- read_excel(filename_in5, sheet = "Lenalidomide", skip = 4)
  both_pbs1$Run <- 1
  both_pbs1$Spreadsheet <- 3
  
  filename_in6 <- paste0(data_dir, "/Plasmaproteinbinding_4h_PBS_results_10182017_Short.xls")
  both_pbs2 <- read_excel(filename_in6, sheet = "Lenalidomide", skip = 4)
  both_pbs2$Run <- 2
  both_pbs2$Spreadsheet <- 3
  
  alldata <- rbind(mouse_plas1, mouse_plas2, human_plas1, human_plas2, both_pbs1, both_pbs2)
  
# Fix and subset data
  name_correction <- c("Specified.Amount", "Calculated.Amount")
  names(alldata) <- str_replace_all(names(alldata),"[ ()#]",".")
  names(alldata)[names(alldata) == "Amount"] <- name_correction
  subdata <- alldata[!is.na(alldata$Sample.Type) & 
    alldata$Sample.Type != "Blank Sample" &
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
  fluidRowCol <- function(x, init) {
    fluidRow( 
      column(6, style = "padding-top: 7px",
        radioButtons(paste0("set", x),
          NULL,
          choices = list(1, 2),
          inline = TRUE,
          selected = init$set[x]
        )  # radioButtons
      ),  # column
      column(6,
        checkboxInput(paste0("status", x),
          NULL,
          value = init$status[x]
        )  # checkboxInput
      )  # column
    )  # fluidRow
  }
  
# Create preset initial values
  initVals <- list(
    Oct2 = list(
      mouse = list(
        set = c(1, 2, rep(1, 7)),
        status = c(F, T, F, rep(T, 6))
      ),
      human = list(
        set = rep(1, 9),
        status = c(F, T, F, rep(T, 6))
      ),
      pbs = list(
        set = rep(2, 9),
        status = c(F, T, F, rep(T, 5), F)
      )
    ),
    Oct18 = list(
      mouse = list(
        set = c(rep(1, 5), rep(2, 4)),
        status = c(F, F, rep(T, 7))
      ),
      human = list(
        set = c(rep(2, 7), 1, 1),
        status = c(T, F, rep(T, 4), F, T, T)
      ),
      pbs = list(
        set = c(rep(1, 5), 2, 1, 1, 1),
        status = c(F, rep(T, 5), F, T, T)
      )
    )
  )