# Testing idea about mean concentrations
  rm(list = ls(all = T))
  
# First load libraries
  library(plyr)
  library(ggplot2)
  library(stringr)
  library(cowplot)
  
# Define forcing functions
  tissue1obs <- function(t) {  # based
    exp(-3.95E-3*t - 5.60) + exp(-2.09E-2*t - 0.286) +
    exp(-0.050*t + 1.85) - exp(-0.1*t)*sum(exp(c(-5.6, -0.286, 1.85)))
  }
  tissue1est1 <- function(t) {
    exp(-1.95E-3*t - 7) + exp(-1.59E-2*t - 0.186) +
    exp(-0.04*t + 2.16) - exp(-0.05*t)*sum(exp(c(-7, -0.186, 2.16)))
  }
  tissue1est2 <- function(t) {
    exp(-3.95E-3*t - 5.60) + exp(-2.09E-2*t - 0.286) +
    exp(-0.050*t + 1.85) - exp(-0.075*t)*sum(exp(c(-5.6, -0.286, 1.85)))
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
  obsdata <- data.frame(TIME = c(0.03, 0.17, 0.33, 0.75, 1, 1.5, 3, 5, 8, 12))
  ERR <- (1+rnorm(n = length(obsdata$TIME), mean = 0, sd = 0.2))
  obsdata$OBS1 <- tissue1obs(obsdata$TIME*60)*ERR
  obsdata$OBS2 <- tissue2obs(obsdata$TIME*60)*ERR
  
# Create data for model predictions
  estdata <- data.frame(TIME = seq(0, 16, by = 0.1))
  estdata$EST1P1 <- tissue1est1(estdata$TIME*60)
  estdata$EST1P2 <- tissue1est2(estdata$TIME*60)
  estdata$EST2P1 <- tissue2est1(estdata$TIME*60)
  estdata$EST2P2 <- tissue2est2(estdata$TIME*60)
  
# Define colourblind-friendly palette
# Red = "#D55E00"
# Blue = "#0072B2"
# Green = "#009E73"
  
# Create plot
  fig1plot <- function(obs, est) {
    p <- NULL
    p <- ggplot()
    p <- p + geom_line(aes(x = TIME, y = get(est)), data = estdata,
      size = 0.8, colour = "#D55E00")
    p <- p + geom_point(aes(x = TIME, y = get(obs)), data = obsdata,
      shape = 4, stroke = 1.2, colour = "#0072B2")
    p <- p + scale_y_log10("\nConcentration (mg/L)\n", lim = c(8.5E-5, 2.5))
    p <- p + scale_x_continuous("\nTime (hours)\n", breaks = 0:4*4)
    p
  }
  
  p1 <- fig1plot("OBS1", "EST1P1")
  p2 <- fig1plot("OBS2", "EST2P1")
  p3 <- fig1plot("OBS1", "EST1P2")
  p4 <- fig1plot("OBS2", "EST2P2")
  
  plot_grid(p1, p2, p3, p4)