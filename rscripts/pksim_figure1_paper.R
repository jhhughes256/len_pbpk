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
  p1 <- p1 + geom_point(aes(x = PRED, y = DV, colour = Tissue),
	  shape = 1, data = castdata)
  p1 <- p1 + geom_abline(intercept = 0, slope = 1, 
    colour = "black", data = castdata)  # add line of identity
  p1 <- p1 + geom_hline(aes(yintercept = LLOQ, group = Tissue), 
    data = lloqdata, colour = cbPalette$pink, linetype = "dashed")
  p1 <- p1 + geom_smooth(aes(x = PRED, y = DV), data = castdata,
	  method = loess, se = T, colour = "red")  # add loess smoothing line
	p1 <- p1 + scale_x_log10("Population Predicted Concentration (ng/mL)", lim = xy.lim2)
	p1 <- p1 + scale_y_log10("Observed Concentration (ng/mL)", lim = xy.lim2)
  p1 <- p1 + scale_colour_manual(name = "Tissue", values = myPalette)
  p1 <- p1 + facet_wrap(~Tissue, ncol = 2)
  p1 <- p1 + theme(legend.position = "none")  # remove legend
  p1
  
# Write plot to file
  ggsave("produced_data/Figure1_paper.png", width = 18.2, height = 24, units = "cm")
  ggsave("produced_data/Figure1_paper.eps", width = 18.2, height = 24, 
    units = "cm", dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  
# Plot 2 - WRES vs TAD - facet_wrap
  y.lim <- max(abs(castdata$WRES), na.rm = T)
  if (y.lim < 0.1) y.lim <- 0.1
      
  p2 <- NULL
  p2 <- ggplot(castdata)
  p2 <- p2 + geom_point(aes(x = TIME, y = WRES, colour = Tissue),
	  shape = 1)
  p2 <- p2 + geom_abline(aes(x = TIME, y = WRES),
	  intercept = 0, slope = 0, colour = "black")  # add line of identity
  p2 <- p2 + geom_smooth(aes(x = TIME, y = WRES),
	  method = loess, se = T, colour = "red")  # add loess smoothing line
  p2 <- p2 + scale_x_continuous("Time After Dose (hours)", 
      breaks = 0:16*60)
  p2 <- p2 + scale_y_continuous("Weighted Residual") # , breaks = -2:2*3)
  p2 <- p2 + coord_cartesian(ylim = c(-y.lim, y.lim))
  p2 <- p2 + scale_colour_manual(values = myPalette)
  p2 <- p2 + facet_wrap(~Tissue, ncol = 2)
  p2 <- p2 + theme(legend.position = "none")  # remove legend
  p2
  
# Plot 3 - Prop RES vs PRED
  p4 <- NULL
  p4 <- ggplot(castdata)
  p4 <- p4 + geom_point(aes(x = PRED, y = WRES, colour = Tissue),
	  shape = 1)  #, alpha = 0.5, size = 1)
  p4 <- p4 + geom_abline(aes(x = PRED, y = WRES),
	  intercept = 0, slope = 0, colour = "black")  # add zero line
  p4 <- p4 + geom_smooth(aes(x = PRED, y = WRES),
	  method = loess, se = T, colour = "red")  # add loess smoothing line
  p4 <- p4 + scale_x_log10("Population Predicted (ng/mL)")
  p4 <- p4 + scale_y_continuous("Weighted Residual") # , breaks = -2:2*3)
  p4 <- p4 + coord_cartesian(ylim = c(-5, 5))
  p4 <- p4 + scale_colour_manual(values = myPalette)
  p4 <- p4 + facet_wrap(~Tissue, ncol = 2, scales = "free_y")
  p4 <- p4 + theme(legend.position = "none")  # remove legend
  p4
  
# Plot 4 - WRES vs TAD - cowplot
# Using cowplot instead of facet_wrap to give each plot its own fixed limits
# Create plot function
  WRESvTIMEplot <- function(data, subtissue, remx, funPalette) {
  # Subset data for tissue of interest
    subdata <- data[data$Tissue == subtissue,]
    
  # Create tissue specific limits
    y.lim <- max(abs(subdata$WRES), na.rm = T)
    if (y.lim < 0.1) y.lim <- 0.1
    
  # Define plot
    p3 <- NULL
    p3 <- ggplot(subdata)
    p3 <- p3 + geom_point(aes(x = TIME, y = WRES, colour = Tissue),
  	  shape = 1)  #, alpha = 0.5, size = 1)
    p3 <- p3 + geom_abline(aes(x = TIME, y = WRES),
  	  intercept = 0, slope = 0, colour = "black")  # add zero line
    p3 <- p3 + geom_smooth(aes(x = TIME, y = WRES),
  	  method = loess, se = T, colour = "red")  # add loess smoothing line
    p3 <- p3 + scale_x_log10(NULL, limits = c(1.5, 500),
      breaks = sort(c(10^(0:2), 3*10^(0:2))))
    p3 <- p3 + ylab(NULL)
    p3 <- p3 + coord_cartesian(ylim = c(-y.lim, y.lim))
    p3 <- p3 + scale_colour_manual(values = funPalette)
    p3 <- p3 + facet_wrap(~Tissue)
    p3 <- p3 + theme(
      legend.position = "none",  # remove legend
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "lines")
    )
    if (remx) {
      p3 <- p3 + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    p3
  }
  
# Create plots and use cowplot to create grid
  p31 <- WRESvTIMEplot(castdata, "Venous Blood Plasma", T, cbPalette$blue)
  p32 <- WRESvTIMEplot(castdata, "Brain Tissue", T, cbPalette$red)
  p33 <- WRESvTIMEplot(castdata, "Heart Tissue", T, cbPalette$green)
  p34 <- WRESvTIMEplot(castdata, "Kidney Tissue", T, cbPalette$orange)
  p35 <- WRESvTIMEplot(castdata, "Liver Tissue", T, cbPalette$pink)
  p36 <- WRESvTIMEplot(castdata, "Lung Tissue", T, cbPalette$grey)
  p37 <- WRESvTIMEplot(castdata, "Muscle Tissue", F, cbPalette$yellow)
  p38 <- WRESvTIMEplot(castdata, "Spleen Tissue", F, cbPalette$skyblue)
  pgrid <- plot_grid(p31, p32, p33, p34, p35, p36, p37, p38, align = "hv", ncol = 2)
  
# Create text grobs for common y and x axis labels
  y.grob <- textGrob("Weighted Residual", vjust = 0.7,
    gp = gpar(fontface = "plain", col = "black", fontsize = 14), rot = 90)
  x.grob <- textGrob("Time After Dose (mins)", hjust = 0.6, vjust = 0.7,
    gp = gpar(fontface = "plain", col = "black", fontsize = 14))

# Produce final figure
  p3 <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(pgrid, left = y.grob, bottom = x.grob)
  )
  plot_grid(p3)
  
# Write plot to file
  ggsave("produced_data/ESMFig1_paper.png", width = 18.2, height = 24, units = "cm")
  # ggsave("produced_data/ESMFig1_paper.eps", width = 18.2, height = 24, 
  #   units = "cm", dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  
# Plot 5 - WRES vs TAD - cowplot
# Using cowplot instead of facet_wrap to give each plot its own fixed limits
# Create plot function
  x.lim <- c(min(castdata$PRED, na.rm = T), max(castdata$PRED, na.rm = T))
  WRESvPREDplot <- function(data, subtissue, remx, funPalette) {
  # Subset data for tissue of interest
    subdata <- data[data$Tissue == subtissue,]
    
  # Create tissue specific limits
    y.lim <- max(abs(subdata$WRES), na.rm = T)
    if (y.lim < 0.1) y.lim <- 0.1
    
  # Define plot
    p5 <- NULL
    p5 <- ggplot(subdata)
    
    p5 <- p5 + geom_point(aes(x = PRED, y = WRES, colour = Tissue),
  	  shape = 1)  #, alpha = 0.5, size = 1)
    p5 <- p5 + geom_abline(aes(x = PRED, y = WRES),
  	  intercept = 0, slope = 0, colour = "black")  # add zero line
    p5 <- p5 + geom_smooth(aes(x = PRED, y = WRES),
  	  method = loess, se = T, colour = "red")  # add loess smoothing line
    p5 <- p5 + scale_x_log10(NULL)
    p5 <- p5 + ylab(NULL)
    p5 <- p5 + coord_cartesian(xlim = x.lim, ylim = c(-y.lim, y.lim))
    p5 <- p5 + scale_colour_manual(values = funPalette)
    p5 <- p5 + facet_wrap(~Tissue)
    p5 <- p5 + theme(
      legend.position = "none",  # remove legend
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "lines")
    )
    if (remx) {
      p5 <- p5 + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    p5
  }
  
# Create plots and use cowplot to create grid
  p51 <- WRESvPREDplot(castdata, "Venous Blood Plasma", T, cbPalette$blue)
  p52 <- WRESvPREDplot(castdata, "Brain Tissue", T, cbPalette$red)
  p53 <- WRESvPREDplot(castdata, "Heart Tissue", T, cbPalette$green)
  p54 <- WRESvPREDplot(castdata, "Kidney Tissue", T, cbPalette$orange)
  p55 <- WRESvPREDplot(castdata, "Liver Tissue", T, cbPalette$pink)
  p56 <- WRESvPREDplot(castdata, "Lung Tissue", T, cbPalette$grey)
  p57 <- WRESvPREDplot(castdata, "Muscle Tissue", F, cbPalette$yellow)
  p58 <- WRESvPREDplot(castdata, "Spleen Tissue", F, cbPalette$skyblue)
  pgrid <- plot_grid(p51, p52, p53, p54, p55, p56, p57, p58, align = "hv", ncol = 2)
  
# Create text grobs for common y and x axis labels
  y.grob <- textGrob("Weighted Residual", vjust = 0.7,
    gp = gpar(fontface = "plain", col = "black", fontsize = 14), rot = 90)
  x.grob <- textGrob("Population Predicted (ng/mL)", hjust = 0.6, vjust = 0.7,
    gp = gpar(fontface = "plain", col = "black", fontsize = 14))

# Produce final figure
  p5 <- gridExtra::grid.arrange(
    gridExtra::arrangeGrob(pgrid, left = y.grob, bottom = x.grob)
  )
  plot_grid(p5)
  
# Write plot to file
  ggsave("produced_data/ESMFig2_paper.png", width = 18.2, height = 24, units = "cm")
  # ggsave("produced_data/ESMFig2_paper.eps", width = 18.2, height = 24, 
  #   units = "cm", dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  