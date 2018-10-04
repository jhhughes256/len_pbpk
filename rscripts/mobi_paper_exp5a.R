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
  
# Define dataset function
  fig4b <- function(k) {
  # Define simulation resolution and forcing functions
    tal <- 5.2
    forcingFunction <- function(t) {
      0.015*(exp(-3.81E-7*t - 1.34) + exp(-0.0145*t + 3.02) + exp(-0.101*t + 5.99))
    }
  # Create base dataset times
    df <- data.frame(TIME = seq(0, 2496/tal, by = 1/tal))
      
  # Define starting concentrations
  # For original forcing function
    df$FF <- forcingFunction(df$TIME)
    df$FFlag <- forcingFunction(df$TIME - 1/tal)
    
  # For all-corrected forcing function
    df$FFhi <- df$FF
    df$FFhi[1] <- 0
    
  # For simres-corrected simulated forcing function
    df$FFlo <- df$FFlag*exp(-k/tal)
    df$FFlo[1] <- 0
    
  # For Input Function PKSim simulation
    df$IF <- (df$FFhi + df$FFlo)/2
    
    df
  }
  
# Run function on different elimination rate constants
  plotbrain <- fig4b(0.00017816/0.00031271)
  plotbrain$TISSUE <- "Brain"
  # plotheart <- fig4c(0.00038372/0.00031271)
  # plotheart$TISSUE <- "Heart"
  # plotkidney <- fig4c(0.0017816/0.00031271)
  # plotkidney$TISSUE <- "Kidney"
  plotlung <- fig4b(0.0075004/0.00071923)
  plotlung$TISSUE <- "Lung"

# Rbind the plots together
  plotdf <- rbind(plotbrain, plotlung)
  
# Define colourblind-friendly palette
# Red = "#D55E00"
# Blue = "#0072B2"
# Green = "#009E73"
  
# Plot data for first plot row
  plotr1fn <- function(cond) {
    p <- NULL
    p <- ggplot()
  
    p <- p + geom_line(aes(x = TIME, y = FFlo), 
      data = plotdf[plotdf$TISSUE == cond,], colour = "#D55E00", size = 0.8)  #, linetype = "dashed")
  
    p <- p + geom_line(aes(x = TIME, y = FF),
      data = plotdf[plotdf$TISSUE == cond,], colour = "black", size = 0.8, linetype = "dashed")
    p <- p + xlab("\nTime (min)")
    if (cond == "Lung") {
      p <- p + scale_y_log10(NULL, breaks = c(0.01, 0.1, 1))
    } else {
      p <- p + scale_y_log10("Concentration (nmol/mL)\n", breaks = c(0.01, 0.1, 1))
    }
    p <- p + coord_cartesian(xlim = c(0, 6*60), ylim = c(0.0005, 7))
    p <- p + facet_wrap(~TISSUE, ncol = 2, scale = "free_y")
    p
  }
  p1 <- plotr1fn("Brain")
  p2 <- plotr1fn("Lung") 
  
# Plot data for second plot
  plotr2fn <- function(cond, ylim) {
    p <- NULL
    p <- ggplot()
    
    p <- p + geom_point(aes(x = TIME, y = FFlo), 
      data = plotdf[plotdf$TISSUE == cond,], colour = "#D55E00", shape = 4, stroke = 1, size = 1)
    p <- p + geom_point(aes(x = TIME, y = FFhi), 
      data = plotdf[plotdf$TISSUE == cond,], colour = "#0072B2", shape = 4, stroke = 1, size = 1)
    p <- p + geom_point(aes(x = TIME, y = IF), 
      data = plotdf[plotdf$TISSUE == cond,], colour = "#009E73", shape = 4, stroke = 1, size = 1)
    
    p <- p + geom_line(aes(x = TIME, y = FF),
      data = plotdf[plotdf$TISSUE == cond,], colour = "black", size = 0.8, linetype = "dashed")
    p <- p + xlab("\nTime (min)")
    if (cond == "Lung") {
      p <- p + ylab(NULL)
    } else {
      p <- p + ylab("Concentration (nmol/mL)\n")
    }
    p <- p + coord_cartesian(xlim = c(0, 1*5), ylim = ylim)
    p <- p + facet_wrap(~TISSUE, ncol = 2, scale = "free_y")
    p
  }
  p3 <- plotr2fn("Brain", c(3, 6.5))
  # p3 <- plotr2fn(plotdf$TISSUE == "Brain", NULL)
  p4 <- plotr2fn("Lung", NULL) 
  
# Source data for third plot
  source("rscripts/mobi_paper_exp3b.R")
  
# Plot data for third plot
  # Plot data
  plotr3fn <- function(cond1, cond2, ylim) {
    p <- NULL
    p <- ggplot()
    
    p <- p + geom_line(aes(x = TIME, y = PREDnd), 
      data = plotdf2[cond2,], colour = "#D55E00", size = 0.8)  #, linetype = "dashed")
    p <- p + geom_line(aes(x = TIME, y = PREDhi), 
      data = plotdf2[cond2,], colour = "#0072B2", size = 0.8)  #, linetype = "dashed")
    p <- p + geom_line(aes(x = TIME, y = PREDav),
      data = plotdf2[cond2,], colour = "#009E73", size = 0.8)  #, linetype = "dashed")
    
    p <- p + geom_point(aes(x = TIME, y = FFlo), 
      data = plotdf[plotdf$TISSUE == cond1,], colour = "#D55E00", shape = 4, stroke = 0.8, size = 1)
    p <- p + geom_point(aes(x = TIME, y = FFhi), 
      data = plotdf[plotdf$TISSUE == cond1,], colour = "#0072B2", shape = 4, stroke = 0.8, size = 1)
    p <- p + geom_point(aes(x = TIME, y = IF), 
      data = plotdf[plotdf$TISSUE == cond1,], colour = "#009E73", shape = 4, stroke = 0.8, size = 1)
  
    p <- p + geom_line(aes(x = TIME, y = FF),
      data = plotdf[plotdf$TISSUE == cond1,], colour = "black", size = 0.8, linetype = "dashed")
    p <- p + xlab("\nTime (min)")
    if (cond1 == "Lung") {
      p <- p + ylab(NULL)
    } else {
      p <- p + ylab("Concentration (nmol/mL)\n")
    }
    p <- p + coord_cartesian(xlim = c(0, 0.6), ylim = ylim)
    p <- p + facet_wrap(~TISSUE, ncol = 2, scale = "free_y")
    p
  }
  p5 <- plotr3fn("Brain", plotdf2$TISSUE == "Brain", c(4.9, 7.4))
  p6 <- plotr3fn("Lung", plotdf2$TISSUE == "Lung", NULL) 
  
# Plot together using cowplot
  plot_grid(p1, p2, p3, p4, p5, p6, align = "hv", ncol = 2, 
    labels = c("a)", "", "b)", "", "c)", ""), 
    hjust = c(-3.5, 0, -3.4, 0, -3.5, 0)
  )
  
# Produce Figure4 for the MoBi methods paper
  ggsave("produced_data/mobi_paper_exp5.png", width = 17.4, height = 19.4,
    units = c("cm"))
  ggsave("produced_data/Figure5.eps", width = 17.4, height = 19.4,
    units = c("cm"), dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
  