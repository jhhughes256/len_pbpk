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
  
# Set up environment with observed data present
  source(paste(git.dir, reponame, "rscripts", "data_iv.R", sep = "/"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy simulation data

# Simulation data is output from MoBi
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol

# Define simulation file names
  file.dir <- "raw_data/MoBi_paper/"
  nomusc.model.5 <- "Fig1Data_NoMusclePGP5.xlsx"
  musc.model.5 <- "Fig1Data_MusclePGP5.xlsx"

# Read in simulation data and subset/rename columns
  raw.nomusc <- read_excel(paste0(file.dir, nomusc.model.5))  # read
  nomusc.data <- as.data.frame(raw.nomusc)[, c(1:2, 26)]  # subset
  names(nomusc.data) <- c("TIME", "Plasma", "Tissue")  # rename
  
  raw.musc <- read_excel(paste0(file.dir, musc.model.5))  # read
  musc.data <- as.data.frame(raw.musc)[, c(1:2, 26)]  # subset
  names(musc.data) <- c("TIME", "Plasma", "Tissue")  # rename
  
# Differentiate the datasets and combine
  nomusc.data$PGP <- F
  musc.data$PGP <- T
  simdata <- melt(rbind(nomusc.data, musc.data), id.vars = c("TIME", "PGP"),
    variable.name = "Tissue", value.name = "Concentration")
  
# Convert time to mins and concentrations to ng/mL
  simdata$TIME <- simdata$TIME*60
  simdata$Concentration <- simdata$Concentration*259.26
  
# Define obsdata and subset for dosage
  obsdata <- melt(dataiv[c(4, 8:9, 12)], id.vars = c("TIME", "DOSEMGKG"),
    variable.name = "Tissue", value.name = "Concentration")
  plotdata <- obsdata[obsdata$DOSEMGKG == 5, ]
  levels(plotdata$Tissue) <- c("Plasma", "Tissue")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Produce Figure1 for the MoBi methods paper
# Define colourblind-friendly palette
  cbpalette <- c("#D55E00", "#0072B2", "#009E73")  # red, blue, green
  
# Define plot
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = TIME, y = Concentration, colour = Tissue), 
    data = simdata, size = 0.8)
  p1 <- p1 + geom_point(aes(x = TIME, y = Concentration, colour = Tissue), 
    data = plotdata)
  p1 <- p1 + scale_colour_manual("Compartment", values = cbpalette)
  p1 <- p1 + xlab("Time (min)")
  p1 <- p1 + scale_y_log10("Concentration (ng/mL)")
  p1 <- p1 + facet_wrap(~PGP, ncol = 2)
  p1
  ggsave("produced_data/mobi_paper_intro4.png", width = 17.4, height = 10.4,
    units = c("cm"))
  # ggsave("produced_data/Figure1.eps", width = 17.4, height = 10.4,
  #   units = c("cm"), dpi = 1200, device = cairo_ps, fallback_resolution = 1200)

