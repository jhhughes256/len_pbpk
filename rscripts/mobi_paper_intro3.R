# Testing idea about mean concentrations
  rm(list = ls(all = T))
  
# First load libraries
  library(plyr)
  library(ggplot2)
  library(stringr)
  library(cowplot)
  
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Define forcing functions
  tissue1obs <- function(t) {  # based
    exp(-3.95E-3*t - 5.6) + exp(-2.09E-2*t - 0.285) + 
    exp(-0.14*t + 1.85) + exp(-4.16*t + 4.39)
  }
  tissue1est1 <- function(t) {
    exp(-2.95E-3*t - 5.9) + exp(-2.09E-2*t - 1.48) + 
    exp(-0.08*t - 0.285)
  }
  tissue1est2 <- function(t) {
    exp(-3.95E-3*t - 5.6) + exp(-2.09E-2*t - 0.285) + 
    exp(-0.14*t + 1.85) + exp(-4.16*t + 4.39)
  }
  
  tissue2obs <- function(t) {
    exp(-3.95E-4*t - 5.60) + exp(-2.09E-3*t - 2.286) +
    exp(-0.02*t + 0.25) - exp(-0.05*t)*sum(exp(c(-5.6, -2.286, 0.25)))
  }
  tissue2est1 <- function(t) {
    exp(-8.95E-4*t - 5.60) + exp(-7.09E-2*t - 2.286) +
    exp(-0.02*t - 0.95) - exp(-0.05*t)*sum(exp(c(-5.6, -2.286, -0.95)))
  }
  tissue2est2 <- function(t) {
    exp(-3.95E-4*t - 5.60) + exp(-2.09E-3*t - 2.286) +
    exp(-0.015*t + 0.25) - exp(-0.07*t)*sum(exp(c(-5.6, -2.286, 0.25)))
  }
  
# Create data for observed data
  set.seed(1337)
  obstime <- c(0.03, 0.17, 0.33, 0.75, 1, 1.5, 3, 5, 8, 12)
  obsdata <- data.frame(TIME = rep(obstime, 2))
  ERR <- (1+rnorm(n = length(obstime)*2, mean = 0, sd = 0.2))
  obsdata$OBS <- c(tissue1obs(obstime*60), tissue2obs(obstime*60))*ERR
  obsdata$TIS <- rep(c("Afferent Blood", "Tissue"), each = length(obstime))
  
# Create data for model predictions
  esttime <- seq(0, 16, by = 0.1)
  estdata <- data.frame(TIME = rep(esttime, 2))
  estdata$EST1 <- c(tissue1est1(esttime*60), tissue1est2(esttime*60))
  estdata$EST2 <- c(tissue2est1(esttime*60), tissue2est2(esttime*60))
  estdata$MOD <- rep(c("Circulatory", "Single-Tissue"), each = length(esttime))
  estdata$TIS1 <- "Afferent Blood"
  estdata$TIS2 <- "Tissue"
  
# Combine data for plotting
  plotdata <- data.frame(TIME = rep(estdata$TIME, 2))
  plotdata$EST <- c(estdata$EST1, estdata$EST2)
  plotdata$MOD <- rep(estdata$MOD, 2)
  plotdata$TIS <- c(estdata$TIS1, estdata$TIS2)
  
# Define colourblind-friendly palette
  cbpalette <- c("#D55E00", "#0072B2", "#009E73")  # red, blue, green
  
# Create plot
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = TIME, y = EST, colour = MOD), 
    data = plotdata, size = 0.8)
  p <- p + geom_point(aes(x = TIME, y = OBS), data = obsdata)
  p <- p + scale_colour_manual("Model", values = cbpalette)
  p <- p + scale_y_log10("\nConcentration (mg/L)\n", lim = c(8.5E-5, 2.5))
  p <- p + scale_x_continuous("\nTime (hours)\n", breaks = 0:4*4)
  p <- p + facet_wrap(~TIS)
  p
  
# Produce Figure1 for the MoBi methods paper
  ggsave("produced_data/mobi_paper_intro3.png", width = 17.4, height = 10.4,
    units = c("cm"))
  ggsave("produced_data/Figure1.eps", width = 17.4, height = 10.4,
    units = c("cm"), dpi = 1200, device = cairo_ps, fallback_resolution = 1200)

