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
  library(ggplot2)
  library(scales)

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Source data
  source("rscripts/mobi_dataprep.R")
  
# The workspace now contains:
#   obsdata - PKSim data melted and cast into afferent, efferent and tissue
#     columns for each of the four PKSim compartments. The tissue type is in 
#     long format for use with facet_wrap().
#   simdata - MoBi data rbound into columns that match that of obsdata. Again,
#     tissue type is in long format.
  
# Units for all data are:
#   time - minutes
#   concentration - nmol/mL
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# First create our plot datasets (pdata)
  desired.columns <- c("TIME", "plAfferent", "bcAfferent", "tiTissue", "TISSUE")
  obspdata <- obsdata[, desired.columns]
  simpdata <- simdata[, desired.columns]
  simpdata$RES <- simpdata$tiTissue - obspdata$tiTissue
  simpdata$PROPRES <- simpdata$RES/obspdata$tiTissue*100
    
# p1 - Plot tissue concentrations against time
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_point(aes(x = TIME, y = tiTissue),
    data = obspdata, colour = "blue", shape = 2)
  p1 <- p1 + geom_line(aes(x = TIME, y = tiTissue),
    data = simpdata, colour = "red")
  p1 <- p1 + xlab("Time (min)")
  # p1 <- p1 + scale_y_log10("Concentration (nmol/mL)")
  p1 <- p1 + facet_wrap(~TISSUE)
  p1
  
# p2 - Plot residuals against time
  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_point(aes(x = TIME, y = RES), data = simpdata)
  p2 <- p2 + xlab("Time (min)")
  p2 <- p2 + scale_y_continuous("Residuals (nmol/mL)")
  p2 <- p2 + facet_wrap(~TISSUE, scales = "free_y")
  p2
  
# p3 - Plot residuals against concentration
  p3 <- NULL
  p3 <- ggplot()
  p3 <- p3 + geom_point(aes(x = tiTissue, y = RES), data = simpdata)
  p3 <- p3 + xlab("Concentration (nmol/mL)")
  p3 <- p3 + scale_y_continuous("Residuals (nmol/mL)")
  p3 <- p3 + facet_wrap(~TISSUE, scales = "free")
  p3
  
# p4 - Plot proportional residuals against time
  p4 <- NULL
  p4 <- ggplot()
  p4 <- p4 + geom_point(aes(x = TIME, y = PROPRES), data = simpdata)
  p4 <- p4 + geom_hline(yintercept = 0, linetype = "dashed", colour = "green4")
  p4 <- p4 + xlab("Time (min)")
  p4 <- p4 + scale_y_continuous("Percent Error (%)")
  p4 <- p4 + facet_wrap(~TISSUE, scales = "free_y")
  p4
  