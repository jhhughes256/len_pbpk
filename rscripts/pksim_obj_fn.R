# Simulation processing script for Lenalidomide PKSim IV Model
# -----------------------------------------------------------------------------
# Incorporates observed data which was previously published:
# Rozewski DM, Herman SEM, Towns WH, Mahoney E, Stefanovski MR, Shin JD, et al.
# Pharmacokinetics and Tissue Disposition of Lenalidomide in Mice. AAPS J.
# 2012;14(4):872-82.
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

# Set up environment with observed data present
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# Load additional libraries
  library(plyr)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Simulation Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Simulation data is output from PKSim
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Read in simulation data

  filename <- list(
    ZeroBBB.HydroBrain = list(
      Dose0_5 = "raw_data/PKSim_data/M.HydroBrain.PGP.0.5.NoGall.xlsx",
      Dose1_5 = "raw_data/PKSim_data/M.HydroBrain.PGP.1.5.NoGall.xlsx",
      Dose5 = "raw_data/PKSim_data/M.HydroBrain.PGP.5.NoGall.xlsx",
      Dose10 = "raw_data/PKSim_data/M.HydroBrain.PGP.10.NoGall.xlsx"
    ),
    NoBBB.NoHydro = list(
      Dose0_5 = "raw_data/PKSim_data/M.SecrPGP.NoBBB.0.5.xlsx",
      Dose1_5 = "raw_data/PKSim_data/M.SecrPGP.NoBBB.1.5.xlsx",
      Dose5 = "raw_data/PKSim_data/M.SecrPGP.NoBBB.5.xlsx",
      Dose10 = "raw_data/PKSim_data/M.SecrPGP.NoBBB.10.xlsx"
    ),
    NoBBB.HydroAll = list(
      Dose0_5 = "raw_data/PKSim_data/M.HydroLit.0.5.xlsx",
      Dose1_5 = "raw_data/PKSim_data/M.HydroLit.1.5.xlsx",
      Dose5 = "raw_data/PKSim_data/M.HydroLit.5.xlsx",
      Dose10 = "raw_data/PKSim_data/M.HydroLit.10.xlsx"
    ),
    zeroBBB.HydroBAll = list(
      Dose0_5 = "raw_data/PKSim_data/M.HydroBrainALL.PGP.0.5.xlsx",
      Dose1_5 = "raw_data/PKSim_data/M.HydroBrainALL.PGP.1.5.xlsx",
      Dose5 = "raw_data/PKSim_data/M.HydroBrainALL.PGP.5.xlsx",
      Dose10 = "raw_data/PKSim_data/M.HydroBrainALL.PGP.10.xlsx"
    )
  )
  
  llply(filename, function(x) {
  # check to see if NoGall is in the title (needs to be removed if so)
    gall.check <- is.na(as.numeric(substr(x[[1]], nchar(x[[1]]) - 5, nchar(x[[1]]) - 5)))
    
  # Start up read excel loop
    read.in <- lapply(x, function(y) {
      df <- as.data.frame(read_excel(y))  # read in excel as data.frame
      if (gall.check) {  # if NoGall then remove NoGall for next step
        y <- paste0(substr(y, 1, nchar(y) - 11), "xlsx")
      }
      mgkg <- as.numeric(substr(y, nchar(y) - 7, nchar(y) - 5))  # determine dose
      if (is.na(mgkg)) mgkg <- 5  # if dose is na then it is 5
      df$DOSEMGKG <- mgkg  # fill column in the data with this info
    # change names so rbind can work before binding the data.frames together in a list
      names(df) <- str_replace(names(df), substr(y, 21, nchar(y) - 5), "")
      df
    })
    
  # List of data.frames is easily bound together with do.call
    simdata_raw <- do.call("rbind", read.in)
    
  # Remove hydrolysis metabolite columns
    simdata <- simdata_raw[!str_detect(names(simdata_raw), "Hydrolysis.Metabolite")]
    
  # Rename columns
    names(simdata) <- c(
      "TIME", "Venous Blood Plasma", "bone", "brain_pl", "brain_bc", "brain_is", 
      "brain_ic", "Brain Tissue", "Heart Tissue", "urine", "Kidney Tissue", 
      "Liver Tissue", "Lung Tissue", "Muscle Tissue", "Spleen Tissue", "DOSEMGKG" 
    )
    names(dataiv) <- c(
      "UID", "ID", "TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME", 
      "Venous Blood Plasma", "Brain Tissue", "Liver Tissue", "Muscle Tissue", 
      "Heart Tissue", "Spleen Tissue", "Lung Tissue", "Kidney Tissue"
    )
    
  # Update time units to minutes and concentration units to ng/mL (ug/L)
    simdata$TIME <- simdata$TIME*60
    umol2ug <- function(x) {
      x*259.26
    }
    simdata[-c(1, 16)] <- colwise(umol2ug)(simdata[-c(1, 16)])
    
  # Melt the data for relevant concentration columns
    plotdata <- melt(simdata[c(1,2,8,9,11:16)], id.vars = c("TIME", "DOSEMGKG"),
      variable.name = "Tissue", value.name = "Concentration")
    plotdata$TIME <- round(plotdata$TIME, 1)
    obsdata <- melt(dataiv[c(3:4, 9:16)], id.vars = c("TADNOM", "DOSEMGKG"),
      variable.name = "Tissue", value.name = "Concentration")
    names(obsdata) <- c("TIME", "DOSEMGKG", "TISSUE", "OBS")
    
    predata <- plotdata[plotdata$TIME %in% unique(obsdata$TIME), ]
    names(predata) <- c("TIME", "DOSEMGKG", "TISSUE", "PRED")
    
    alldata <- merge(obsdata, predata)
    cleandata <- alldata[!is.na(alldata$OBS),]
    
  # Subset to see fits of individual tissues if interested
    cleandata <- cleandata[cleandata$TISSUE == "Heart Tissue",]
    
    obj <- -2*sum(with(cleandata, dnorm(OBS, PRED, abs(PRED)*0.3, log = T)))
  })
