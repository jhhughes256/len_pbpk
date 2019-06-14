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
  library(scales)
  library(cowplot)

# Update ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Simulation Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Simulation data is output from PKSim
# Model names are: 
# M.BBB  M.BBB.HydroAll  M.BBB.HydroBrain  M.LitPGP  M.NoBBB  
# M.NoBBB.HydroAll  M.NoBBB.HydroPlas  M.BBB.HydroComb
# Remember to choose the model name with suffix "diag"
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Read in simulation data
  file.dir <- "raw_data/PKSim_paper/"
  model.name <- "M.NoBBB.HydroPlas"
  
  simdata0_5 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".Diag.0.5.xlsx")))
  simdata1_5 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".Diag.1.5.xlsx")))
  simdata5 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".Diag.5.xlsx")))
  simdata10 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".Diag.10.xlsx")))
  
# Specify dosage for each dataset
  simdata0_5$DOSEMGKG <- 0.5
  simdata1_5$DOSEMGKG <- 1.5
  simdata5$DOSEMGKG <- 5
  simdata10$DOSEMGKG <- 10
  
# Ensure names of columns are the same to allow rbind
  names(simdata0_5) <- str_replace(names(simdata0_5), paste0(model.name, ".0.5"), "")
  names(simdata1_5) <- str_replace(names(simdata1_5), "sim", "")
  names(simdata5) <- str_replace(names(simdata5), "sim", "")
  names(simdata10) <- str_replace(names(simdata10), "sim", "")
  simdata_raw <- rbind(simdata0_5, simdata1_5, simdata5, simdata10)
  
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
  obsdata <- melt(dataiv[c(2:4, 9:16)], id.vars = c("TADNOM", "DOSEMGKG", "ID"),
    variable.name = "Tissue", value.name = "Concentration")

# LLOQ data
# Regular lloq is 1 nmol/L == 259.26 ng/L == 0.25926 ng/mL
  lloqdata <- data.frame(
    Tissue = c(
      "Venous Blood Plasma", "Brain Tissue", "Liver Tissue", "Muscle Tissue", 
      "Heart Tissue", "Spleen Tissue", "Lung Tissue", "Kidney Tissue"),
    LLOQ = 0.25926*c(1, 3.11, 3.14, 0.643, 0.937, 1.55, 1.09, 1.92)
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Merge plotdata and obsdata
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Subset plotdata and rename obsdata time column
# Round times, as simulation times are only correct to 6 decimal places
  plotdata$TIME <- round(plotdata$TIME)
  subdata <- plotdata[plotdata$TIME %in% unique(dataiv$TADNOM),]
  names(obsdata)[1] <- "TIME"
  
# Duplicate simulation values and give ID column
# Will enable the casting of data later
  subdata <- rbind(
    data.frame(subdata, ID = 1),
    data.frame(subdata, ID = 2),
    data.frame(subdata, ID = 3),
    data.frame(subdata, ID = 4),
    data.frame(subdata, ID = 5)
  )
  
# Add column to define data type and then bind together
  subdata$DATA <- "PRED"
  obsdata$DATA <- "DV"
  alldata <- rbind(subdata, obsdata)
  
# Cast the data.frame to give two columns, one for DV and one for PRED
  castdata <- dcast(alldata, TIME + DOSEMGKG + Tissue + ID ~ DATA, 
    value.var = "Concentration")
  
# Remove all rows with no observed values and calculate residuals
  castdata <- castdata[!is.na(castdata$DV),]
  castdata$EPS <- with(castdata, log(DV) - log(PRED))
  popvar <- with(castdata, 1/1165*sum(EPS^2))
  castdata$W <- with(castdata, sqrt(log(PRED)^2*popvar^2))
  castdata$WRES <- with(castdata, EPS/W)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Diagnostic plots using cowplot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define colourblind palette and custom palette
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  myPalette <- with(cbPalette, 
    c(blue, red, green, orange, pink, grey, yellow, skyblue)
  )

# Plot 1 - OBS vs PRED (log scale)
  xy.lim2 <- c(min(c(castdata$PRED, castdata$DV), na.rm = T), 
    max(c(castdata$PRED, castdata$DV), na.rm = T))

  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_point(aes(x = PRED, y = DV), colour = cbPalette$blue,
	  shape = 1, data = castdata[castdata$Tissue == "Venous Blood Plasma",])
  p1 <- p1 + geom_abline(intercept = 0, slope = 1, 
    colour = "black")  # add line of identity
  p1 <- p1 + geom_smooth(aes(x = PRED, y = DV), data = castdata,
	  method = loess, se = T, colour = "red")  # add loess smoothing line
	p1 <- p1 + scale_x_log10(NULL, lim = xy.lim2)
	p1 <- p1 + scale_y_log10(NULL, lim = xy.lim2)
  p1 <- p1 + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )  # remove legend
  p1