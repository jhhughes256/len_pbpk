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
  dfbrain <- fig4b(0.00017816/0.00031271)
  dfbrain$TISSUE <- "Brain"
  dflung <- fig4b(0.0075004/0.00071923)
  dflung$TISSUE <- "Lung"
  
# Rbind the data together
  dfall <- rbind(dfbrain, dflung)

# Create plot data
  plotdf <- melt(dfall, id.vars = c("TIME", "TISSUE"), 
    variable.name = "Type", value.name = "DV")
  levels(plotdf$Type) <- c("FF", "FFlag", "AC", "NC", "TC")
  
# Define colourblind-friendly palette
  cbpalette <- c("#0072B2", "#D55E00", "#009E73")  # red, blue, green
  
# Plot data for first plot row
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = TIME, y = DV), colour = "#D55E00", size = 0.8,
    data = plotdf[plotdf$Type == "NC",])
  p1 <- p1 + geom_line(aes(x = TIME, y = DV), colour = "black", size = 0.8, 
    data = plotdf[plotdf$Type == "FF",], linetype = "dashed")
  p1 <- p1 + xlab("\nTime (min)")
  p1 <- p1 + scale_y_log10("Concentration (nmol/mL)\n", breaks = c(0.01, 0.1, 1))
  p1 <- p1 + coord_cartesian(xlim = c(0, 6*60), ylim = c(0.0005, 7))
  p1 <- p1 + facet_wrap(~TISSUE, ncol = 2)
  p1
  
# Plot data for second plot
  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_point(aes(x = TIME, y = DV, colour = Type), size = 1, 
    data = plotdf[plotdf$Type %in% c("NC", "AC", "TC"),], shape = 4, stroke = 1)
  p2 <- p2 + geom_line(aes(x = TIME, y = DV), colour = "black", size = 0.8, 
    data = plotdf[plotdf$Type == "FF",], linetype = "dashed")
  p2 <- p2 + scale_colour_manual("", values = c(cbpalette))
  p2 <- p2 + guides(colour = F)
  p2 <- p2 + xlab("\nTime (min)")
  p2 <- p2 + ylab("Concentration (nmol/mL)\n")
  p2 <- p2 + coord_cartesian(xlim = c(0, 5))
  p2 <- p2 + facet_wrap(~TISSUE, ncol = 2)
  p2
  
# Source data for third plot
  source("rscripts/mobi_paper_exp3b.R")
  
# Plot data for third plot
  # Plot data
  plotdf <- dfall
  plotr3fn <- function(cond1, cond2, ylim) {
    p <- NULL
    p <- ggplot()
    
    p <- p + geom_line(aes(x = TIME, y = PREDnd), 
      data = plotdf2[cond2,], colour = "#D55E00", size = 0.8)  #, linetype = "dashed")
    p <- p + geom_line(aes(x = TIME, y = PREDhi), 
      data = plotdf2[cond2,], colour = "#0072B2", size = 0.8)  #, linetype = "dashed")
    p <- p + geom_line(aes(x = TIME, y = PREDmx),
      data = plotdf2[cond2,], colour = "#009E73", size = 0.8)  #, linetype = "dashed")
    
    p <- p + geom_point(aes(x = TIME, y = FFlo), 
      data = plotdf[cond1,], colour = "#D55E00", shape = 4, stroke = 1, size = 1)
    p <- p + geom_point(aes(x = TIME, y = FFhi), 
      data = plotdf[cond1,], colour = "#0072B2", shape = 4, stroke = 1, size = 1)
    p <- p + geom_point(aes(x = TIME, y = IF), 
      data = plotdf[cond1,], colour = "#009E73", shape = 4, stroke = 1, size = 1)
  
    p <- p + geom_line(aes(x = TIME, y = FF),
    data = plotdf[cond1,], colour = "black", size = 0.8, linetype = "dashed")
    p <- p + xlab("\nTime (min)")
    if (is.null(ylim)) {
      p <- p + ylab(NULL)
    } else {
      p <- p + ylab("Concentration (nmol/mL)\n")
    }
    p <- p + coord_cartesian(xlim = c(0, 0.6), ylim = ylim)
    p <- p + facet_wrap(~TISSUE, ncol = 2, scale = "free_y")
    p
  }
  p3a <- plotr3fn(plotdf$TISSUE == "Brain", plotdf2$TISSUE == "Brain", c(4.9, 7.4))
  p3b <- plotr3fn(plotdf$TISSUE == "Lung", plotdf2$TISSUE == "Lung", NULL) 
  p3 <- plot_grid(p3a, p3b, align = "h", rel_widths = c(1.1, 1))
  p3
  
# Plot together using cowplot
  plot_grid(p1, p2, p3, align = "hv", axis = "l", ncol = 1)
  
# Produce Figure5 for the MoBi methods paper
  ggsave("produced_data/mobi_paper_exp4.png", width = 17.4, height = 19.4,
    units = c("cm"))
  # ggsave("produced_data/Figure5.eps", width = 17.4, height = 19.4,
  #   units = c("cm"), dpi = 1200, device = cairo_ps, fallback_resolution = 1200)

  ### THIS PLOT WAS NOT SUCCESSFULLY IMPLEMENTED ###