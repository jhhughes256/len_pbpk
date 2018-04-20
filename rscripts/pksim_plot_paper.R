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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Simulation Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Simulation data is output from PKSim
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Read in simulation data
  # filename0_5 <- "raw_data/PKSim_data/M.HydroBrain.PGP.0.5.NoGall.xlsx"
  # filename1_5 <- "raw_data/PKSim_data/M.HydroBrain.PGP.1.5.NoGall.xlsx"
  # filename5 <- "raw_data/PKSim_data/M.HydroBrain.PGP.5.NoGall.xlsx"
  # filename10 <- "raw_data/PKSim_data/M.HydroBrain.PGP.10.NoGall.xlsx"
  # filename0_5 <- "raw_data/PKSim_data/M.SecrPGP.NoBBB.0.5.xlsx"
  # filename1_5 <- "raw_data/PKSim_data/M.SecrPGP.NoBBB.1.5.xlsx"
  # filename5 <- "raw_data/PKSim_data/M.SecrPGP.NoBBB.5.xlsx"
  # filename10 <- "raw_data/PKSim_data/M.SecrPGP.NoBBB.10.xlsx"
  # filename0_5 <- "raw_data/PKSim_data/M.HydroLit.0.5.xlsx"
  # filename1_5 <- "raw_data/PKSim_data/M.HydroLit.1.5.xlsx"
  # filename5 <- "raw_data/PKSim_data/M.HydroLit.5.xlsx"
  # filename10 <- "raw_data/PKSim_data/M.HydroLit.10.xlsx"
  filename0_5 <- "raw_data/PKSim_data/M.HydroBrainALL.PGP.0.5.xlsx"
  filename1_5 <- "raw_data/PKSim_data/M.HydroBrainALL.PGP.1.5.xlsx"
  filename5 <- "raw_data/PKSim_data/M.HydroBrainALL.PGP.5.xlsx"
  filename10 <- "raw_data/PKSim_data/M.HydroBrainALL.PGP.10.xlsx"
  
  simdata0_5 <- as.data.frame(read_excel(filename0_5))
  simdata1_5 <- as.data.frame(read_excel(filename1_5))
  simdata5 <- as.data.frame(read_excel(filename5))
  simdata10 <- as.data.frame(read_excel(filename10))
  
# Specify dosage for each dataset
  simdata0_5$DOSEMGKG <- 0.5
  simdata1_5$DOSEMGKG <- 1.5
  simdata5$DOSEMGKG <- 5
  simdata10$DOSEMGKG <- 10
  
# Ensure names of columns are the same to allow rbind
  # names(simdata0_5) <- str_replace(names(simdata0_5), "M.HydroBrain.PGP.0.5", "")
  # names(simdata1_5) <- str_replace(names(simdata1_5), "M.HydroBrain.PGP.1.5", "")
  # names(simdata5) <- str_replace(names(simdata5), "M.HydroBrain.PGP.5", "")
  # names(simdata10) <- str_replace(names(simdata10), "M.HydroBrain.PGP.10", "")
  # names(simdata0_5) <- str_replace(names(simdata0_5), "M.SecrPGP.NoBBB.0.5", "")
  # names(simdata1_5) <- str_replace(names(simdata1_5), "M.SecrPGP.NoBBB.1.5", "")
  # names(simdata5) <- str_replace(names(simdata5), "M.SecrPGP.NoBBB.5", "")
  # names(simdata10) <- str_replace(names(simdata10), "M.SecrPGP.NoBBB.10", "")
  # names(simdata0_5) <- str_replace(names(simdata0_5), "M.HydroLit.0.5", "")
  # names(simdata1_5) <- str_replace(names(simdata1_5), "M.HydroLit.1.5", "")
  # names(simdata5) <- str_replace(names(simdata5), "M.HydroLit.5", "")
  # names(simdata10) <- str_replace(names(simdata10), "M.HydroLit.10", "")
  names(simdata0_5) <- str_replace(names(simdata0_5), "M.HydroBrainALL.PGP.0.5", "")
  names(simdata1_5) <- str_replace(names(simdata1_5), "M.HydroBrainALL.PGP.1.5", "")
  names(simdata5) <- str_replace(names(simdata5), "M.HydroBrainALL.PGP.5", "")
  names(simdata10) <- str_replace(names(simdata10), "M.HydroBrainALL.PGP.10", "")
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
  obsdata <- melt(dataiv[c(4, 8:16)], id.vars = c("TIME", "DOSEMGKG"),
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
# Plot data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
  plotdata$DOSEMGKGf <- factor(plotdata$DOSEMGKG)
  obsdata$DOSEMGKGf <- factor(obsdata$DOSEMGKG)+
  
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = TIME, y = Concentration, colour = DOSEMGKGf), 
    data = plotdata)
  p <- p + geom_point(aes(x = TIME, y = Concentration, colour = DOSEMGKGf), 
    data = obsdata, alpha = 0.3)
  p <- p + geom_hline(aes(yintercept = LLOQ), 
    data = lloqdata, linetype = "dashed", colour = "magenta")
  p <- p + facet_wrap(~Tissue, ncol = 2, scales = "free")
  p <- p + xlab("\nTime (mins)")
  p <- p + scale_y_log10("Concentration (ng/mL)\n", labels = comma, breaks = c(1e4, 1e2, 1, 1e-2))
  p <- p + scale_colour_manual(name = "Dose (mg/kg)", 
    values = c("red", "green4", "blue", "purple"))
  p <- p + coord_cartesian(xlim = (c(0, 800)))
  p
  
  setwd("produced_data")
  ggsave("Figure 2_HydroAll.png", width = 23.2, height = 17.2, units = "cm")
  