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
# Model names are: 
# M.BBB  M.BBB.HydroAll  M.BBB.HydroBrain  M.LitPGP  M.NoBBB  
# M.NoBBB.HydroAll  M.NoBBB.HydroPlas  M.BBB.HydroComb
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Read in simulation data
  file.dir <- "raw_data/PKSim_paper/"
  model.name <- "M.NoBBB.HydroPlas"
  
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

# Remove brain tissue predictions for dosages with no observed data
  plotdata <- plotdata[!(
    plotdata$Tissue == "Brain Tissue" & plotdata$DOSEMGKG %in% c(0.5, 1.5)
  ),]
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
  plotdata$DOSEMGKGf <- factor(plotdata$DOSEMGKG)
  obsdata$DOSEMGKGf <- factor(obsdata$DOSEMGKG)
  
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
  myPalette <- with(cbPalette, c(orange, green, red, blue))
  
# Define plot
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = TIME, y = Concentration, colour = DOSEMGKGf), 
    data = plotdata)
  p <- p + geom_point(aes(x = TIME, y = Concentration, colour = DOSEMGKGf), 
    data = obsdata, alpha = 0.3)
  p <- p + geom_hline(aes(yintercept = LLOQ), 
    data = lloqdata, linetype = "dashed", colour = cbPalette$pink)
  p <- p + facet_wrap(~Tissue, ncol = 2, scales = "free")
  p <- p + theme(legend.direction = "horizontal", legend.position = c(-0.02, -0.1),
    legend.background = element_rect(fill = "white", colour = NA))
  p <- p + xlab("\nTime (mins)")
  p <- p + scale_y_log10("Concentration (ng/mL)\n", labels = comma, 
    breaks = c(1e4, 1e2, 1, 1e-2))
  p <- p + scale_colour_manual(name = "       Dose (mg/kg)", values = myPalette,
    guide = guide_legend(title.position = "top"))
  p <- p + coord_cartesian(xlim = (c(0, 500)))
  p
  
  ggsave(paste0("produced_data/Figure3_poster.png"), 
    width = 17.4, height = 17.2, units = "cm")
  
  # ggsave(paste0("produced_data/Figure3_poster.eps"), 
  #   width = 17.4, height = 17.2, units = "cm",
  #   dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  