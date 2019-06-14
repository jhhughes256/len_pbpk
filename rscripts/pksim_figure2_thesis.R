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

# Set thesis ggplot theme
 	theme_bw2 <- theme_set(theme_bw(14))
	theme_bw2 <- theme_update(
		plot.margin = unit(c(0.5, 0.5, 1, 0.5), "lines")
	)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Simulation Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Simulation data is output from PKSim
# Model names are: 
# M.BBB  M.BBB.HydroAll  M.BBB.HydroBrain  M.LitPGP  M.NoBBB  
# M.NoBBB.HydroAll  M.NoBBB.HydroPlas  M.BBB.HydroComb
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Create function to read in simulation data
  MakePlotData <- function(file.dir, model.name) {
    simdata0_5 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".0.5.xlsx")))
    simdata1_5 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".1.5.xlsx")))
    simdata5 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".5.xlsx")))
    simdata10 <- as.data.frame(read_excel(paste0(file.dir, model.name, ".10.xlsx")))
    
  # Specify dosage for each dataset
    simdata0_5$DOSEMGKG <- 0.5
    simdata1_5$DOSEMGKG <- 1.5
    simdata5$DOSEMGKG <- 5
    simdata10$DOSEMGKG <- 10
    
  # Ensure names of columns are the same to allow rbind
    names(simdata0_5) <- str_replace(names(simdata0_5), paste0(model.name, ".0.5"), "")
    names(simdata1_5) <- str_replace(names(simdata1_5), paste0(model.name, ".1.5"), "")
    names(simdata5) <- str_replace(names(simdata5), paste0(model.name, ".5"), "")
    names(simdata10) <- str_replace(names(simdata10), paste0(model.name, ".10"), "")
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
    plotdata$SIM <- model.name
    plotdata
  }
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Make plots
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Use function to read in plotdata
  plotPaper <- MakePlotData("raw_data/PKSim_paper/", "M.Paper")
  plotLitPGP <- MakePlotData("raw_data/PKSim_paper/", "M.LitPGP")
  plotBBB <- MakePlotData("raw_data/PKSim_paper/", "M.BBB.HydroBrain")
  plotNoBBB <- MakePlotData("raw_data/PKSim_paper/", "M.NoBBB")
  alldata <- rbind(plotPaper, plotLitPGP, plotBBB, plotNoBBB)
  brainpred <- alldata[alldata$Tissue == "Brain Tissue",]
  
# Quick fix to script to make brainobs subset work
  names(dataiv) <- c(
    "UID", "ID", "TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME", 
    "Venous Blood Plasma", "Brain Tissue", "Liver Tissue", "Muscle Tissue", 
    "Heart Tissue", "Spleen Tissue", "Lung Tissue", "Kidney Tissue"
  )
  
# Create melted dataiv dataset to represent observed values, subset for brain
  obsdata <- melt(dataiv[c(4, 8:16)], id.vars = c("TIME", "DOSEMGKG"),
    variable.name = "Tissue", value.name = "Concentration")
  brainobs <- obsdata[obsdata$Tissue == "Brain Tissue",]
  brainobs.old <- brainobs
  brainobs <- brainobs[brainobs$DOSEMGKG %in% c(5, 10),]
  brainobs$DOSEMGKGf <- factor(brainobs$DOSEMGKG)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
  brainpred.old <- brainpred
  brainpred <- brainpred[brainpred$DOSEMGKG %in% c(5, 10),]
  brainpred$DOSEMGKGf <- factor(brainpred$DOSEMGKG)
  brainpred$SIMf <- factor(brainpred$SIM, 
    levels = levels(factor(brainpred$SIM))[c(4, 2, 1, 3)])
  levels(brainpred$SIMf) <- c("A - No P-gp", "B - P-gp Expression", 
    "C - Brain Metabolism", "D - Intracellular Transport")

  brainobs$DOSEMGKGf <- factor(brainobs$DOSEMGKG)
  
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
  myPalette <- with(cbPalette, c(red, blue))
  
# Define plot
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = TIME, y = Concentration, colour = DOSEMGKGf), 
    data = brainpred)
  p <- p + geom_point(aes(x = TIME, y = Concentration, colour = DOSEMGKGf), 
    data = brainobs, alpha = 0.3)
  p <- p + geom_hline(yintercept = 0.25926*3.11, linetype = "dashed", 
    colour = cbPalette$pink)
  p <- p + facet_wrap(~SIMf, ncol = 2)
  p <- p + xlab("\nTime (min)")
  p <- p + scale_y_log10("Concentration (ng/mL)\n", labels = comma, 
    breaks = c(1e4, 1e2, 1, 1e-2))
  p <- p + scale_colour_manual(name = "Dose (mg/kg)", values = myPalette)
  p <- p + coord_cartesian(xlim = c(0, 300), ylim = c(0.001, 1000))
  p
  
  ggsave(paste0("produced_data/Figure2_thesis.png"), 
    width = 17.4, height = 16.5, units = "cm")