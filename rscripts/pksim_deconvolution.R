# Testing numerical deconvolution from PKSim afferent blood concentrations to
# forcing function concentrations
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

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare PKSim output and corresponding forcing function input
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Read in file
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "BrainSumexp52.xlsx"
  braindata_raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Remove hydrolysis metabolite columns (if they exist)
  braindata <- braindata_raw[!str_detect(names(braindata_raw), "Hydrolysis.Metabolite")]
  
# Rename columns
  names(braindata) <- c(
    "TIME", "Venous Blood Plasma", "Venous Blood Blood Cells", 
    "Arterial Blood Plasma", "Arterial Blood Blood Cells", "Tissue Plasma", 
    "Tissue Blood Cells", "Tissue Interstitial", "Tissue Intracellular", "Whole Tissue"
  )
  
# Correct times from hours to minutes
  braindata$TIME <- braindata$TIME*60
  
# Subset for number of samples desired
  simdata <- braindata
  
# Create forcing function data and bind with PKSim data and times
  forcingFunction <- function(x) {
    0.015*(exp(-3.81E-7*x - 1.34) + exp(-0.0145*x + 3.02) + exp(-0.101*x + 5.99))
  }
  # nolagdf <- data.frame(
  #   time = simdata$TIME,
  #   fin = forcingFunction(simdata$TIME),
  #   out = simdata$`Arterial Blood Plasma`
  # )
  lagdf <- data.frame(
    time = tail(simdata$TIME, -1),
    fin = forcingFunction(head(simdata$TIME, -1)),
    out = tail(simdata$`Arterial Blood Plasma`, -1)
  )

# Quick look at the data
  p <- NULL
  p <- ggplot(data = lagdf)
  p <- p + geom_point(aes(x = time, y = fin), colour = "red")
  p <- p + geom_line(aes(x = time, y = out), linetype = "dashed")
  p <- p + coord_trans(y = "log10")
  p 
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Deconvolution
# First add zero padding to the data
  indv <- c(0, 0, lagdf$fin, 0, 0)
  intime <- c(-2, -1, lagdf$time, max(lagdf$time) + 1, max(lagdf$time) + 2)
  oudv <- c(0, 0, lagdf$out, 0, 0)
  outime <- intime
  npoints <- length(lagdf$fin)
  
# fft requires equal numbers of evenly spaced data points
  
# Deconvolution in the time domain is division in the frequency domain
  # unit = output/input
  
# Change to frequency domain
  fftfin <- fft(lagdf$fin)
  fftout <- fft(lagdf$out)
  
# Deconvolve in frequency domain to find the unit response function
  unitfft <- (fftout/fftfin)/npoints
  
# The unit response function affects concentrations at each unit of time
  # fine for imaginary parts to be discarded at this step
  unit <- as.double(fft(unitfft, inverse = T))
  
  plot(unit[-1] ~ lagdf$time[-1], col = "red", type = "b", main = "Unit Response Function")
  
  