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

# Load package libraries
	library(mrgsolve)
	library(reshape2)
	library(ggplot2)

# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

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

	TIME <- seq(from = 0, to = 60, by = 0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Simulation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source the model
	source(paste(master.dir, "models", "mouse_brown.R", sep = "/"))

# Simulate concentration-time profiles for the population
	ID <- 1:n
	ID2 <- sort(c(rep(ID, times = length(TIME))))
	times <- rep(TIME, times = length(ID))
	input.concdata <- data.frame(
		ID = ID2,
		time = times,
		amt = 0,
		evid = 0,
		rate = 0,
		cmt = 1
	)

	dose.times <- 0
	dosedata <- input.concdata[input.concdata$time %in% dose.times, ]
	dosedata$amt <- 10
	dosedata$evid <- 1
	dosedata$rate <- 1
	dosedata$cmt <- 1

	input.concdata <- rbind(input.concdata, dosedata)
	input.concdata <- input.concdata[with(input.concdata, order(ID, time)), ]

	concdata <- as.data.frame(mrgsim(
		data_set(brown.mod, input.concdata)
	))  # mrgsim

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#	Plot Simulated Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Determine number of model compartments
  init.str <- names(brown.mod@init)
  init.n <- length(init.str)
	plotdata <- cbind(
		melt(concdata[1:(2+init.n)], c("ID", "time"), variable.name = "tissue"),
		melt(concdata[c(1:2, (init.n+3):(2+init.n*2))], c("ID", "time"))[-(1:3)]
	)
	names(plotdata) <- c("ID", "TIME", "COMP", "A", "C")
	plotdata$COMP <- as.factor(plotdata$COMP)
	levels(plotdata$COMP) <- c(
    toupper(substr(init.str, 2, nchar(init.str)))
  )

	plotobj <- NULL
	plotobj <- ggplot(data = plotdata)
	plotobj <- plotobj + geom_line(aes(x = TIME, y = C), colour = "blue")
	plotobj <- plotobj + facet_wrap(~ COMP, ncol = 4)
	plotobj <- plotobj + scale_y_continuous("Concentrations (mg/mL)\n")  #scale_y_log10("Concentrations (mg/mL)\n")
	plotobj <- plotobj + scale_x_continuous("\nTime (mins)")
	plotobj
