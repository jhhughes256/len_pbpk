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
#	Simulation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source the model
  source(paste(master.dir, "models", "mouse_brown.R", sep = "/"))

# Simulate concentration-time profiles for the population
  ID <- 1:n
  ID2 <- sort(c(rep(ID, times = length(TIME))))
  times <- rep(TIME, times = length(ID))
  input.simdata <- data.frame(
    ID = ID2,
    time = times,
    amt = 0,
    evid = 0,
    rate = 0,
    cmt = 1
  )

  dose.times <- 0
  dosedata <- input.simdata[input.simdata$time %in% dose.times, ]
  dosedata$amt <- 1
  dosedata$evid <- 1
  dosedata$rate <- 60
  dosedata$cmt <- 1

  input.simdata <- rbind(input.simdata, dosedata)
  input.simdata <- input.simdata[with(input.simdata, order(ID, time)), ]

  simdata <- as.data.frame(mrgsim(
    data_set(brown.mod, input.simdata)
  ))  # mrgsim

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Plot Observed Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  obsdata <- dataiv[c("ID", "TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT", "PLA", "BRA",
    "LVR", "MSC", "HRT", "SPL", "LUN", "KID")]
  names(obsdata)[names(obsdata) == "PLA"] <- "PA"
  names(obsdata)[names(obsdata) == "SPL"] <- "SPR"
  names(obsdata)[names(obsdata) == "LUN"] <- "ART"
  obsdata.plot <- melt(obsdata,
    c("ID", "TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT"),
    variable.name = "COMP", value.name = "C"
  )

  # dose normalised concentrations

  p0 <- NULL
  p0 <- ggplot(data = obsdata.plot)
  p0 <- p0 + geom_point(aes(x = TIME, y = C), colour = "blue")
  p0 <- p0 + facet_wrap(~ COMP, ncol = 4)
  p0 <- p0 + scale_y_log10("Concentrations (ng/mL)\n", labels = scales::comma)  #scale_y_log10("Concentrations (mg/mL)\n")
  p0 <- p0 + scale_x_continuous("\nTime (mins)")
  p0

  avedata <- dataiv.av[c("TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT", "PLA", "BRA",
    "LVR", "MSC", "HRT", "SPL", "LUN", "KID")]
  names(avedata)[names(avedata) == "PLA"] <- "PA"
  names(avedata)[names(avedata) == "SPL"] <- "SPR"
  names(avedata)[names(avedata) == "LUN"] <- "ART"
  avedata$SPS <- avedata$SPR
  avedata.plot <- melt(avedata,
    c("TIME", "TADNOM", "DOSEMGKG", "DOSEMG", "WT"),
    variable.name = "COMP"
  )
  avedata.plot$C <- avedata.plot$value/(avedata.plot$DOSEMG*10^6)

  p1 <- NULL
  p1 <- ggplot(data = avedata.plot)
  p1 <- p1 + geom_point(aes(x = TIME, y = C), colour = "blue")
  p1 <- p1 + facet_wrap(~ COMP, ncol = 4)
  p1 <- p1 + scale_y_log10("Concentrations (ng/mL)\n", labels = scales::comma)  #scale_y_log10("Concentrations (mg/mL)\n")
  p1 <- p1 + scale_x_continuous("\nTime (mins)")
  p1

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Plot Simulated Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Determine number of model compartments
  init.str <- names(brown.mod@init)
  init.n <- length(init.str)
  simdata.plot <- cbind(
    melt(simdata[1:(2+init.n)], c("ID", "time"), variable.name = "tissue"),
    melt(simdata[c(1:2, (init.n+3):(2+init.n*2))], c("ID", "time"))[-(1:3)]
  )
  names(simdata.plot) <- c("ID", "TIME", "COMP", "A", "C")
  simdata.plot$COMP <- as.factor(simdata.plot$COMP)
  levels(simdata.plot$COMP) <- c(
    toupper(substr(init.str, 2, nchar(init.str)))
  )
  simdata.plot <- simdata.plot[simdata.plot$TIME != 0, ]

  p2 <- NULL
  p2 <- ggplot(data = simdata.plot[!simdata.plot$COMP %in% c("GUT", "TUBF", "TUBC", "BOD"),])
  p2 <- p2 + geom_line(aes(x = TIME, y = C), colour = "blue")
  p2 <- p2 + facet_wrap(~ COMP, ncol = 4)
  p2 <- p2 + scale_y_log10("Concentrations (ng/mL)\n", labels = scales::comma)  #scale_y_log10("Concentrations (mg/mL)\n")
  p2 <- p2 + scale_x_continuous("\nTime (mins)")
  p2

  p3 <- NULL
  p3 <- ggplot(data = simdata.plot[!simdata.plot$COMP %in% c("GUT", "TUBF", "TUBC", "BOD"),])
  p3 <- p3 + geom_line(aes(x = TIME, y = C), colour = "blue")
  p3 <- p3 + geom_point(aes(x = TIME, y = C), data = avedata.plot, colour = "red", alpha = 0.2)
  p3 <- p3 + facet_wrap(~ COMP, ncol = 4)
  p3 <- p3 + scale_y_log10("Concentrations (ng/mL)\n", labels = scales::comma)  #scale_y_log10("Concentrations (mg/mL)\n")
  p3 <- p3 + scale_x_continuous("\nTime (mins)", lim = c(0, 100))
  p3

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Plots for Shiny App
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate concentration-time profiles for the population
  ID <- 1:4
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
  dosedata$amt <- c(0.5, 1.5, 5, 10)*unique(input.appdata$WT)*10^3
  dosedata$evid <- 1
  dosedata$rate <- max(dosedata$amt)*60
  dosedata$cmt <- 1

  input.appdata <- rbind(input.appdata, dosedata)
  input.appdata <- input.appdata[with(input.appdata, order(ID, time)), ]

  appdata <- as.data.frame(mrgsim(
    data_set(brown.mod, input.appdata)
  ))  # mrgsim

  appdata$dosemgkg <- factor(appdata$ID)
  levels(appdata$dosemgkg) <- c(0.5, 1.5, 5, 10)

  init.str <- names(brown.mod@init)
  init.n <- length(init.str)
  n.cols <- length(names(appdata))
  appdata.plot <- cbind(
    melt(appdata[c(1:(2+init.n), n.cols)], c("ID", "time", "dosemgkg"), variable.name = "tissue"),
    melt(appdata[c(1:2, (init.n+3):(2+init.n*2), n.cols)], c("ID", "time", "dosemgkg"))[-(1:4)]
  )
  names(appdata.plot) <- c("ID", "TIME", "DOSEMGKGf", "COMP", "A", "C")
  appdata.plot$COMP <- as.factor(appdata.plot$COMP)
  levels(appdata.plot$COMP) <- c(
    toupper(substr(init.str, 2, nchar(init.str)))
  )
  appdata.plot <- appdata.plot[appdata.plot$TIME != 0, ]

  obsdata.plot$DOSEMGKGf <- factor(obsdata.plot$DOSEMGKG)

  p4 <- NULL
  p4 <- ggplot(data = appdata.plot[appdata.plot$COMP == "PA",])
  p4 <- p4 + geom_line(aes(x = TIME, y = C), colour = "blue")
  p4 <- p4 + geom_point(aes(x = TIME, y = C), data = obsdata.plot[obsdata.plot$COMP == "PA",], colour = "red", alpha = 0.2)
  p4 <- p4 + facet_wrap(~ DOSEMGKGf, ncol = 2)
  p4 <- p4 + scale_y_log10("Concentrations (ng/mL)\n", labels = scales::comma)  #scale_y_log10("Concentrations (mg/mL)\n")
  p4 <- p4 + scale_x_continuous("\nTime (mins)", lim = c(0, 100))
  p4

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Data for Shiny App
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  shinydata1 <- obsdata.plot[-c(1,9)]
  shinydata2 <- avedata.plot[-8]
  names(shinydata2)[names(shinydata2) == "value"] <- "C"
  shinydata1$TYPE <- 0  # all data
  shinydata2$TYPE <- 1  # mean data
  shinydata <- rbind(shinydata1, shinydata2)
  write.csv(shinydata, "rscripts/model_app/data.csv", row.names = F)
