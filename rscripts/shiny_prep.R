# R script for running models
# -----------------------------------------------------------------------------
# Ready workspace
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list=ls(all=TRUE))
    graphics.off()
    #git.dir <- "E:/Hughes/Git"
    git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    #git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "len_pbpk"
    scriptname <- "run_model.R"
  }

# Set working directory
  master.dir <- paste(git.dir, reponame, sep = "/")
  setwd(master.dir)

# Source observed data
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"), verbose = F)
  source(paste(git.dir, reponame, "rscripts",
    "data_po.R", sep = "/"), verbose = F)
  source(paste(git.dir, reponame, "rscripts",
    "data_ip.R", sep = "/"), verbose = F)

# Load package libraries
  library(mrgsolve)
  library(reshape2)
  library(ggplot2)

# Define a custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 16))
  theme_bw2 <- theme_update(axis.text.x = element_text(angle = 35, hjust = 0.8))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set number of individuals that make up the 95% prediction intervals
  n <- 1
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
  CI95lo <- function(x) quantile(x, probs = 0.025)
  CI95hi <- function(x) quantile(x, probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)
# Set seed for reproducible numbers
  set.seed(123456)

  TIME <- seq(from = 0, to = 500, by = 0.2)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Simulate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source the model
  source(paste(master.dir, "models", "mouse_brown.R", sep = "/"))

# Simulate concentration-time profiles for the population
  ID <- 1:11
  ID2 <- sort(c(rep(ID, times = length(TIME))))
  times <- rep(TIME, times = length(ID))
  input.appdata <- data.frame(
    ID = ID2,
    time = times,
    amt = 0,
    evid = 0,
    rate = 0,
    cmt = 1,
    WT = 28
  )

  dose.times <- 0
  dosedata <- input.appdata[input.appdata$time %in% dose.times, ]
  dosedata$amt <- c(
    c(0.5, 1.5, 5, 10)*unique(input.appdata$WT)*10^3,
    c(0.5, 10)*unique(input.appdata$WT)*10^3,
    c(0.5, 10)*unique(input.appdata$WT)*10^3,
    1, 1, 1
  )
  dosedata$evid <- 1
  dosedata$rate <- max(dosedata$amt)*60
  dosedata$cmt <- c(rep(1, 4), rep(2, 2), rep(4, 2), c(1, 2, 4))

  input.appdata <- rbind(input.appdata, dosedata)
  input.appdata <- input.appdata[with(input.appdata, order(ID, time)), ]

  appdata <- as.data.frame(mrgsim(
    data_set(brown.mod, input.appdata)
  ))  # mrgsim

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Observed Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IV concentrations
  obsdata.iv <- dataiv[c("ID", "TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT",
    "PLA", "BRA", "LVR", "MSC", "HRT", "SPL", "LUN", "KID")]
  names(obsdata.iv)[names(obsdata.iv) == "PLA"] <- "PA"
  names(obsdata.iv)[names(obsdata.iv) == "SPL"] <- "SPR"
  names(obsdata.iv)[names(obsdata.iv) == "LUN"] <- "ART"
  obsdata.iv$SPS <- obsdata.iv$SPR
  obsdata.iv$DATA <- 1
  obsdata.iv.plot <- melt(obsdata.iv,
    c("ID", "TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT", "DATA"),
    variable.name = "COMP", value.name = "C"
  )
  obsdata.iv.plot$TYPE <- 0

# IV average concentrations
  avedata.iv <- dataiv.av[c("TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT", "PLA",
    "BRA","LVR", "MSC", "HRT", "SPL", "LUN", "KID")]
  names(avedata.iv)[names(avedata.iv) == "PLA"] <- "PA"
  names(avedata.iv)[names(avedata.iv) == "SPL"] <- "SPR"
  names(avedata.iv)[names(avedata.iv) == "LUN"] <- "ART"
  avedata.iv$SPS <- avedata.iv$SPR
  avedata.iv$DATA <- 1
  avedata.iv.plot <- melt(avedata.iv,
    c("TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT", "DATA"),
    variable.name = "COMP",
    value.name = "C"
  )
  avedata.iv.plot$TYPE <- 1

# PO concentrations
  obsdata.po.plot <- datapo[c("ID", "TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT")]
  obsdata.po.plot$COMP <- factor("PA")
  obsdata.po.plot$DATA <- 0
  obsdata.po.plot$C <- datapo$DV
  obsdata.po.plot$TYPE <- 0

# PO average concentrations
  avedata.po.plot <- datapo.av[c("TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT")]
  avedata.po.plot$COMP <- factor("PA")
  avedata.po.plot$DATA <- 0
  avedata.po.plot$C <- datapo.av$DV
  avedata.po.plot$TYPE <- 1

# IP Concentrations
  obsdata.ip.plot <- dataip[c("ID", "TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT")]
  obsdata.ip.plot$COMP <- factor("PA")
  obsdata.ip.plot$DATA <- 2
  obsdata.ip.plot$C <- dataip$DV
  obsdata.ip.plot$TYPE <- 0

# IP average concentrations
  avedata.ip.plot <- dataip.av[c("TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT")]
  avedata.ip.plot$COMP <- factor("PA")
  avedata.ip.plot$DATA <- 2
  avedata.ip.plot$C <- dataip.av$DV
  avedata.ip.plot$TYPE <- 1

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Plots for Shiny App
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  appdata$dosemgkg <- factor(appdata$ID)
  levels(appdata$dosemgkg) <- c(0.5, 1.5, 5, 10, 0, 0.5, 10, 0)

  appdata$iv <- factor(appdata$ID)
  levels(appdata$iv) <- c(rep(1, 5), rep(0, 3))

  init.str <- names(brown.mod@init)
  init.n <- length(init.str)
  n.cols <- length(names(appdata))
  appdata.plot <- cbind(
    melt(appdata[c(1:(2+init.n), n.cols -1, n.cols)], c("ID", "time", "dosemgkg", "iv"), variable.name = "tissue"),
    melt(appdata[c(1:2, (init.n+3):(2+init.n*2), n.cols -1, n.cols)], c("ID", "time", "dosemgkg", "iv"))[-(1:5)]
  )
  names(appdata.plot) <- c("ID", "TIME", "DOSEMGKGf", "IVf", "COMP", "A", "C")
  appdata.plot$COMP <- as.factor(appdata.plot$COMP)
  levels(appdata.plot$COMP) <- c(
    toupper(substr(init.str, 2, nchar(init.str)))
  )
  appdata.plot <- appdata.plot[appdata.plot$TIME != 0, ]

  obsdata.iv.plot$DOSEMGKGf <- factor(obsdata.iv.plot$DOSEMGKG)

  p4 <- NULL
  p4 <- ggplot(data = appdata.plot[appdata.plot$COMP == "PA",])
  p4 <- p4 + geom_line(aes(x = TIME, y = C), colour = "blue")
  p4 <- p4 + geom_point(aes(x = TIME, y = C), data = obsdata.iv.plot[obsdata.iv.plot$COMP == "PA",], colour = "red", alpha = 0.2)
  p4 <- p4 + facet_wrap(~ DOSEMGKGf, ncol = 2)
  p4 <- p4 + scale_y_log10("Concentrations (ng/mL)\n", labels = scales::comma)  #scale_y_log10("Concentrations (mg/mL)\n")
  p4 <- p4 + scale_x_continuous("\nTime (mins)", lim = c(0, 100))
  p4

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Data for Shiny App
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  shinydata1 <- obsdata.iv.plot[-1]
  shinydata2 <- avedata.iv.plot
  shinydata3 <- obsdata.po.plot[-1]
  shinydata4 <- avedata.po.plot
  shinydata5 <- obsdata.ip.plot[-1]
  shinydata6 <- avedata.ip.plot
  shinydata <- rbind(shinydata1, shinydata2, shinydata3, shinydata4)
  write.csv(shinydata, "rscripts/modapp-lenapbpk/www/data.csv", row.names = F)
