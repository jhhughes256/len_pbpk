# Validation of forcing functions used in MoBi
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
  model.name <- "BrainSumexp52.xlsx"
  braindata_raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Remove hydrolysis metabolite columns
  braindata <- braindata_raw[!str_detect(names(braindata_raw), "Hydrolysis.Metabolite")]
  
# Rename columns
  names(braindata) <- c(
    "TIME", "Venous Blood Plasma", "Venous Blood Blood Cells", 
    "Arterial Blood Plasma", "Arterial Blood Blood Cells", "Tissue Plasma", 
    "Tissue Blood Cells", "Tissue Interstitial", "Tissue Intracellular", "Whole Tissue"
  )
  
# Clean up
  V <- 0.00031271
  Q <- 0.00017816
  rbcK <- 0.47787
  HCT <- 0.45
  braindata$TIME <- braindata$TIME*60
  simdata <- braindata
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
  simdata$Cp2t <- with(simdata, PREDp*Q/(tickRate*V))
  simdata$Cb2t <- with(simdata, PREDb*Q/(tickRate*V))
  simdata$Ap2t <- with(simdata, Cp2t*V*(1-HCT))
  simdata$Ab2t <- with(simdata, Cb2t*V*HCT)
  simdata$PREDLAGp <- forcingFunction(c(0, head(simdata$TIME, -1)))
  simdata$DV <- with(simdata, PREDp - Cp2t)
  simdata$RESp <- with(simdata, PREDp - `Arterial Blood Plasma`)
  simdata$PROPRESp <- with(simdata, RESp/PREDp)

  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = TIME, y = PROPRESp), data = simdata[-1,])
  p1 <- p1 + geom_hline(yintercept = (tissueQ/bloodV)/tickRate, linetype = "dashed")
  p1
  
  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_line(aes(x = TIME, y = DV), data = simdata)
  p2

  p3 <- NULL
  p3 <- ggplot()
  p3 <- p3 + geom_point(aes(x = `Arterial Blood Plasma`, y = DV), data = simdata)
  p3 <- p3 + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  p3
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hasn't quite worked out, so lets continue this idea, but with the assumption
# that there are DES ticks and then there are reporting ticks!
# Decide number of DES ticks between reporting ticks
  oldTIME <- simdata$TIME
  repTick <- diff(oldTIME)[1]
  nTick <- 20
  newTick <- repTick/nTick
  newTIME <- cumsum(rep(newTick, nTick*(length(oldTIME)-1)))
  newdata <- data.frame(TIME = newTIME)
  
# Set tolerances for DES
  
# Set constants
  Qt <- 1.7816e-4     # tissue blood flow
  fu <- 0.7           # fraction unbound
  HCT <- 0.45         # hematocrit
  Vart <- 3.1271e-4   # volume of arterial compartment
  Vtis <- 2.3297e-4   # volume of tissue compartment
  Vven <- 7.1923e-4   # volume of venous compartment
  ftv <- 0.037        # fraction of vascular - tissue
  fti <- 0.004        # fraction of interstitial - tissue
  Kpb <- 0.47787      # partition coefficient - plasma/blood
  Kpi <- 0.78467      # partition coefficient - plasma/interstitial
  Kiw <- 0.89210      # partition coefficient - interstitial/water
  Kcw <- 1.1552       # partition coefficient - intracellular/water
  SA_apb <- 1.6886    # surface area - arterial plasma/blood
  SA_vpb <- 3.8838    # surface area - venous plasma/blood
  SA_tpb <- 0.046548    # surface area - tissue plasma/blood
  SA_tpi <- 0.0081890    # surface area - tissue plasma/interstitial
  P_endo <- 1.8863e-5 # velocity - tissue endothelium
  PSA_ic <- 5.1648e-8 # velocity*surface area - tissue interstitial/intracellular

  phvec <- c(
    TIME = 0, Cinp = forcingFunction(0), Cinb = forcingFunction(0)*Kpb,
    Cap = 0, Cab = 0, Ctp = 0, Ctb = 0, Cti = 0, Ctc = 0, Cvp = 0, Cvb = 0, 
    Aap2b = 0, Ap2t = 0, Ab2t = 0, Atp2b = 0,  Atp2i = 0, 
    Ati2c = 0, At2p = 0, At2b = 0, Avp2b = 0
  )
  phdata <- as.data.frame(matrix(phvec, nrow = 1))
  names(phdata) <- names(phvec)
  for (i in 1:6) {
    for (j in 1:nTick) {
    ## Determine mass balance
    # But first set time and apply forcing function if required
      phvec["TIME"] <- newTIME[nTick*(i-1)+j]
      if (j == 1) {
        phvec["Cap"] <- phvec["Cinp"]
        phvec["Cab"] <- phvec["Cinb"]
      }
    # Afferent Vessel Equilibrium
      phvec["Aap2b"] <- newTick*SA_apb*fu*(1*phvec["Cap"] - 1*phvec["Cab"]/Kpb)
    # Afferent Vessel to Tissue Transfer
      phvec["Ap2t"] <- newTick*phvec["Cap"]*Qt*(1-HCT)
      phvec["Ab2t"] <- newTick*phvec["Cab"]*Qt*HCT
    # Tissue Equilibrium
      phvec["Atp2b"] <- newTick*SA_tpb*fu*(1*phvec["Ctp"] - 1*phvec["Ctb"]/Kpb)
      phvec["Atp2i"] <- newTick*fu*P_endo*SA_tpi*(phvec["Ctp"] - phvec["Cti"]/Kpi)
      phvec["Ati2c"] <- newTick*PSA_ic*(phvec["Cti"]*Kiw - phvec["Ctc"]*Kcw)
    # Tissue to Efferent Vessel Transfer
      phvec["At2p"] <- newTick*phvec["Ctp"]*Q*(1-HCT)
      phvec["At2b"] <- newTick*phvec["Ctb"]*Q*HCT
    # Efferent Vessel Equilibrium
      phvec["Avp2b"] <- newTick*SA_vpb*fu*(1*phvec["Cvp"] - 1*phvec["Cvb"]/Kpb)
      
    ## Determine new concentrations
      # browser()
    # Arterial Plasma
      phvec["Cap"] <- phvec["Cap"] + sum(
          -phvec["Aap2b"],          # Plasma Blood Cell Equilibrium
          -phvec["Ap2t"]            # Arterial Plasma to Tissue Plasma
        )/(Vart*(1-HCT))            # Arterial Plasma Volume
    # Arterial Blood Cells
      phvec["Cab"] <- phvec["Cab"] + sum(   
          phvec["Aap2b"],           # Plasma Blood Cell Equilibrium
          -phvec["Ab2t"]            # Arterial Blood Cells to Tissue Blood Cells
        )/(Vart*HCT)                # Arterial Blood Cells Volume
    # Tissue Plasma
      phvec["Ctp"] <- phvec["Ctp"] + sum(
          phvec["Ap2t"],            # Arterial Plasma to Tissue Plasma
          -phvec["Atp2b"],          # Plasma Blood Cell Equilibrium
          -phvec["Atp2i"],          # Plasma Interstitial Equilibrium
          -phvec["At2p"]            # Tissue Plasma to Venous Plasma
        )/(Vtis*ftv*(1-HCT))        # Tissue Plasma Volume
    # Tissue Blood Cells
      phvec["Ctb"] <- phvec["Ctb"] + sum(
          phvec["Ab2t"],            # Arterial Blood Cells to Tissue Blood Cells
          phvec["Atp2b"],           # Plasma Blood Cell Equilibrium
          -phvec["At2b"]            # Tissue Blood Cells to Venous Blood Cells
        )/(Vtis*ftv*HCT)            # Tissue Plasma Volume
    # Tissue Interstitial
      phvec["Cti"] <- phvec["Cti"] + sum(
        phvec["Atp2i"],             # Plasma Interstitial Equilibrium
        -phvec["Ati2c"]             # Interstitial Intracellular Equilibrium
      )/(Vtis*fti)                  # Tissue Interstitial Volume
    # Tissue Intracellular
      phvec["Ctc"] <- phvec["Ctc"] + sum(
        phvec["Ati2c"]              # Interstitial Intracellular Equilibrium
      )/(Vtis*(1 - sum(fti, ftv)))  # Tissue Intracellular Volume
    # venous Plasma and Blood Cells
      phvec["Cvp"] <- phvec["Cvp"] + phvec["At2p"]/(Vven*(1-HCT)) 
      phvec["Cvb"] <- phvec["Cvb"] + phvec["At2b"]/(Vven*(1-HCT)) 
  
    ## Ready for next loop
    # Update forcing function if required and placeholder data.frame
      if (j == nTick) {
        phvec["Cinp"] <- forcingFunction(phvec["TIME"])
        phvec["Cinb"] <- phvec["Cinp"]*Kpb
      }
      phdata <- rbind(phdata, phvec)
    }
  }
  # out <- phdata
c
phvec