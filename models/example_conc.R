# R script for simulating a population from a 3 compartment recirculation model
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(dplyr)	# New plyr - required for mrgsolve
	library(mrgsolve)	# Metrum differential equation solver for pharmacometrics
# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# ------------------------------------------------------------------------------
# Set number of individuals that make up the 95% prediction intervals
	n <- 1
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x,probs = 0.05)
	CI90hi <- function(x) quantile(x,probs = 0.95)
# Set seed for reproducible numbers
	set.seed(123456)

	TIME <- seq(from = 0, to = 240, by = 0.1)

# ------------------------------------------------------------------------------
# Define the model parameters and equations
	# Using mrgsolve
	code <- '
$INIT
	Cart = 0,  // Arterial blood
	Chep = 0,  // Liver
	Cbod = 0  // Rest of body
$PARAM
	// Regional blood flow
  Qco = 9,
	Qhep = 2,
	CLhep = 1,
  // Tissue mass balance
  Vlng = 2,
	Vhep = 1,
	Vbod = 25
$MAIN
	double Qbod = Qco - Qhep;
$ODE
	dxdt_Cart = (-Qco*Cart +Qhep*Chep +Qbod*Cbod)/Vlng;
	dxdt_Chep = (Qhep*Cart -Qhep*Chep -CLhep*Cart)/Vhep;
	dxdt_Cbod = (Qbod*Cart - Qbod*Cbod)/Vbod;
$CAPTURE
  Qco Qhep Qbod CLhep Vlng Vhep Vbod
'
	# Compile the model code
	mod <- mcode("3COMPrecirc", code)

# ------------------------------------------------------------------------------
# Simulate concentration-time profiles for the population
	ID <- 1:1
	ID2 <- sort(c(rep(ID, times = length(TIME))))
	time <- rep(TIME, times = length(ID))
	input.conc.data <- data.frame(
		ID = ID2,
		time,
		amt = 0,
		evid = 0,
		rate = 0,
		cmt = 1
	)

	inf.dose.times <- 0
	inf.dose.data <- input.conc.data[input.conc.data$time %in% inf.dose.times,]
	inf.dose.data$amt <- 20000
	inf.dose.data$evid <- 1
	inf.dose.data$rate <- 100
	inf.dose.data$cmt <- 1

	input.conc.data <- rbind(input.conc.data, inf.dose.data)
	input.conc.data <- input.conc.data[with(input.conc.data, order(input.conc.data$ID, input.conc.data$time)), ]

	conc.data <- mod %>% data_set(input.conc.data) %>% mrgsim()
	conc.data <- as.data.frame(conc.data)
