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
# Read in simulation data
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "HumanBrainSumexp52.xlsx"
  braindata_raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Remove hydrolysis metabolite columns
  braindata <- braindata_raw[!str_detect(names(braindata_raw), "Hydrolysis.Metabolite")]
  
# Rename columns
  names(braindata) <- c(
    "TIME", "Venous Blood Plasma", "Venous Blood Blood Cells", 
    "Arterial Blood Plasma", "Arterial Blood Blood Cells", "Tissue Plasma", 
    "Tissue Blood Cells", "Tissue Interstitial", "Tissue Intracellular", "Whole Tissue"
  )
  
# Read in simulation data
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "HumanHeartSumexp52.xlsx"
  heartdata_raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Remove hydrolysis metabolite columns
  heartdata <- heartdata_raw[!str_detect(names(heartdata_raw), "Hydrolysis.Metabolite")]
  
# Rename columns
  names(heartdata) <- c(
    "TIME", "Venous Blood Plasma", "Venous Blood Blood Cells", 
    "Arterial Blood Plasma", "Arterial Blood Blood Cells", "Tissue Plasma", 
    "Tissue Blood Cells", "Tissue Interstitial", "Tissue Intracellular", "Whole Tissue"
  )
  
# Read in simulation data
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "HumanKidneySumexp52.xlsx"
  kidneydata_raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Remove hydrolysis metabolite columns
  kidneydata <- kidneydata_raw[!str_detect(names(kidneydata_raw), "Hydrolysis.Metabolite")]
  
# Rename columns
  names(kidneydata) <- c(
    "TIME", "Venous Blood Plasma", "Venous Blood Blood Cells", 
    "Arterial Blood Plasma", "Arterial Blood Blood Cells", "Tissue Plasma", 
    "Tissue Blood Cells", "Tissue Interstitial", "Tissue Intracellular", "Whole Tissue"
  )
  
# Read in simulation data
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "HumanLungSumexp52.xlsx"
  lungdata_raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Remove hydrolysis metabolite columns
  lungdata <- lungdata_raw[!str_detect(names(lungdata_raw), "Hydrolysis.Metabolite")]
  
# Rename columns
  names(lungdata) <- c(
    "TIME", "Venous Blood Plasma", "Venous Blood Blood Cells", 
    "Arterial Blood Plasma", "Arterial Blood Blood Cells", "Tissue Plasma", 
    "Tissue Blood Cells", "Tissue Interstitial", "Tissue Intracellular", "Whole Tissue"
  )
  
# Add tissue column and bind together
  braindata$TISSUE <- "Brain"
  braindata$K <- 0.70331/0.41015
  heartdata$TISSUE <- "Heart"
  heartdata$K <- 0.29700/0.41015
  kidneydata$TISSUE <- "Kidney"
  kidneydata$K <- 0.98668/0.41015
  lungdata$TISSUE <- "Lung"
  lungdata$K <- 5.5327/0.94206
  simdata1 <- rbind(braindata[-1,], heartdata[-1,], kidneydata[-1,])
  simdata1$OBS <- simdata1$`Arterial Blood Plasma`
  simdata2 <- lungdata[-1,]
  simdata2$OBS <- simdata2$`Venous Blood Plasma`
  simdata <- rbind(simdata1, simdata2)
  simdata$TIME <- simdata$TIME*60
  simdata$TISSUEf <- factor(simdata$TISSUE)
  
# Update time units to minutes and concentration units to ng/mL (ug/L)
  # umol2ug <- function(x) {
  #   x*259.26
  # }
  # simdata[-1] <- colwise(umol2ug)(simdata[-1])
  rbcK <- 0.47787
  tickRate <- 1/diff(simdata$TIME)[1]
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validate Forcing Function
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Blood flows of each tissue are
  forcingFunction <- function(x) {
    0.015*(exp(-3.81E-7*x - 1.34) + exp(-0.0145*x + 3.02) + exp(-0.101*x + 5.99))
  }
  simdata$PREDp <- forcingFunction(simdata$TIME)
  simdata$PREDb <- simdata$PREDp*rbcK
  simdata$RESp <- with(simdata, PREDp - OBS)
  simdata$PROPRESp <- with(simdata, RESp/PREDp)
    simdata <- ddply(simdata, .(TISSUE), function(x) {
    x$PREDLAGp <- forcingFunction(c(0, head(x$TIME, -1)))
    x$PREDLAGb <- x$PREDLAGp*rbcK
    x
  })
  # simdata$Cp2t <- with(simdata, PREDLAGp*K/tickRate)
  # simdata$Cb2t <- with(simdata, PREDLAGb*K/tickRate)
  simdata$Cp2t <- with(simdata, PREDLAGp*(1 - exp(-K/tickRate)))
  simdata$Cb2t <- with(simdata, PREDLAGb*(1 - exp(-K/tickRate)))

  simdata$DV <- with(simdata, PREDLAGp - Cp2t)
  simdata$DVRES <- with(simdata, DV - OBS)
  simdata$DVRESPROP <- with(simdata, abs(DVRES)/DV)
  
  p0 <- NULL
  p0 <- ggplot()
  p0 <- p0 + geom_line(aes(x = TIME, y = OBS, colour = TISSUE), data = simdata[-1,])
  p0 <- p0 + coord_trans(y = "log10")
  p0

  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = PREDp, y = RESp, colour = TISSUE), data = simdata[-1,])
  p1

  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_line(aes(x = TIME, y = PROPRESp, colour = TISSUE), data = simdata[-1,])
  p2 <- p2 + geom_hline(aes(yintercept = K/tickRate, colour = TISSUE), linetype = "dashed", data = simdata[-1])
  p2 <- p2 + coord_trans(y = "log10")
  p2
  
  p3 <- NULL
  p3 <- ggplot()
  p3 <- p3 + geom_point(aes(x = OBS, y = DV, colour = TISSUE), size = 1, data = simdata)
  p3 <- p3 + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  p3
  
  p4 <- NULL
  p4 <- ggplot()
  p4 <- p4 + geom_line(aes(x = PREDp, y = PROPRESp, colour = TISSUE), data = simdata[-1,])
  p4

# Check out the slope of the forcing function vs the residuals
  brainmod <- lm(RESp ~ PREDp, data = simdata[simdata$TISSUE == "Brain",])
  kidneymod <- lm(RESp ~ PREDp, data = simdata[simdata$TISSUE == "Kidney",])
  brainmod$coefficients["PREDp"]
  kidneymod$coefficients["PREDp"]
  