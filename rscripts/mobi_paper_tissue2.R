# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  library(plyr)
  library(ggplot2)
  library(scales)

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

# Source data then redefine objects to allow them to be row bound
  
# The workspace will contain:
#   obsdata - PKSim data melted and cast into afferent, efferent and tissue
#     columns for each of the four PKSim compartments. The tissue type is in 
#     long format for use with facet_wrap().
#   simdata - MoBi data rbound into columns that match that of obsdata. Again,
#     tissue type is in long format.
  
# Units for all data are:
#   time - minutes
#   concentration - nmol/mL

# First redefine mouse data
  source("rscripts/mobi_dataprep_base.R")
  mobsdata <- obsdata
  msimdata <- simdata
  
# Then the human data
  source("rscripts/mobi_dataprep_basehuman.R")
  hobsdata <- obsdata
  hsimdata <- simdata
  
# Give the datasets a species column before binding together
  mobsdata$SPECIES <- "Mouse"
  msimdata$SPECIES <- "Mouse"
  hobsdata$SPECIES <- "Human"
  hsimdata$SPECIES <- "Human"
  obsdata <- rbind(mobsdata, hobsdata)
  simdata <- rbind(msimdata, hsimdata)
  
# Create a joint column describing both the tissue and species
  obsdata$MODEL <- with(obsdata, paste(SPECIES, TISSUE))
  simdata$MODEL <- with(simdata, paste(SPECIES, TISSUE))

# First create our plot datasets (pdata)
  desired.columns <- c("TIME", "plAfferent", "bcAfferent", "tiTissue", "MODEL", 
    "TISSUE", "SPECIES")
  obspdata <- obsdata[, desired.columns]
  simpdata <- simdata[, desired.columns]
  simpdata$RES <- simpdata$tiTissue - obspdata$tiTissue
  simpdata$PROPRES <- simpdata$RES/obspdata$tiTissue*100
  
# Create data.frame to dictate the geom_text arguments
# Want human geom_text to be in the lower half of the plot, mouse in the upper
  textdata <- ddply(simpdata, .(MODEL, SPECIES), function(df) {
    meanerr <- mean(df$PROPRES, na.rm = T)
    errtext <- paste0(signif(meanerr, 3), "%")
    out <- data.frame(meanERR = errtext)
    if (unique(df$SPECIES) == "Human") {
      out$yaxis <- -25
    } else {  # if unique(df$SPECIES == "Mouse")
      out$yaxis <- 25
    }
    out
  })
  
# Define colourblind palette
  cbPalette <- c("#D55E00", "#E69F00", "#CC79A7", "#009E73", "#0072B2")
    
# Plot proportional residuals against time
  p <- NULL
  p <- ggplot()
  p <- p + geom_point(aes(x = TIME/60, y = PROPRES, colour = TISSUE),
    shape = 1, data = simpdata)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "green4")
  p <- p + geom_text(aes(label = meanERR, y = yaxis), x = 850/60, data = textdata)
  p <- p + scale_colour_manual("Model", values = cbPalette)
  p <- p + scale_x_continuous("Time (h)", breaks = c(0, 4, 8, 12, 16))
  p <- p + scale_y_continuous("Percent Error (%)", lim = c(-50, 50))
  p <- p + facet_wrap(~MODEL, ncol = 2, dir = "v")
  p <- p + theme(legend.position = "none")  # remove legend
  p
  
# Produce Figure 3 for the MoBi methods paper
  ggsave("produced_data/mobi_paper_tissue2.png", width = 17.4, height = 17.4, 
    units = c("cm"))
  ggsave("produced_data/Figure3.eps", width = 17.4, height = 23.4,
    units = c("cm"), dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  
  ggsave("produced_data/mobi_paganz_tissue2.png", width = 17.4, height = 12.4, 
    units = c("cm"))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extra Plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # p1 - Plot tissue concentrations against time
  # p1 <- NULL
  # p1 <- ggplot()
  # p1 <- p1 + geom_point(aes(x = TIME, y = tiTissue),
  #   data = obspdata, colour = "blue", shape = 1)
  # p1 <- p1 + geom_line(aes(x = TIME, y = tiTissue),
  #   data = simpdata, colour = "red")
  # p1 <- p1 + xlab("Time (min)")
  # p1 <- p1 + scale_y_log10("Concentration (nmol/mL)")
  # p1 <- p1 + facet_wrap(~MODEL, ncol = 4)
  # p1
   
# p2 - Plot residuals against time
  # p2 <- NULL
  # p2 <- ggplot()
  # p2 <- p2 + geom_point(aes(x = TIME, y = RES), data = simpdata)
  # p2 <- p2 + xlab("Time (min)")
  # p2 <- p2 + scale_y_continuous("Residuals (nmol/mL)")
  # p2 <- p2 + facet_wrap(~MODEL, ncol = 4, scales = "free_y")
  # p2
   
# p3 - Plot residuals against concentration
  # p3 <- NULL
  # p3 <- ggplot()
  # p3 <- p3 + geom_point(aes(x = tiTissue, y = RES), data = simpdata)
  # p3 <- p3 + xlab("Concentration (nmol/mL)")
  # p3 <- p3 + scale_y_continuous("Residuals (nmol/mL)")
  # p3 <- p3 + facet_wrap(~MODEL, ncol = 4, scales = "free")
  # p3
  