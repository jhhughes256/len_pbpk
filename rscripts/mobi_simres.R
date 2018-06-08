# Testing the relationship between sensitivity index and simulation resolution
# -----------------------------------------------------------------------------
# Ready workspace
  rm(list=ls(all=TRUE))

# Read in data from simulation resolution test
  res <- c(0.33, 0.65, 1.3, 2.6, 5.2)
  si <- c(1.72, 0.88, 0.44, 0.22, 0.11)
  
# Looks asymptotic
  matres <- matrix(1/res, nrow = length(res))
  ones <- matrix(rep(1, times = length(res)), nrow = length(res))
  X <- cbind(ones, matres)
  mod <- lm(si ~ 0 + X)
  summary(mod)
  
# r2 of 1, wow!
# lets plot that
  modpar <- mod$coefficients
  simres <- seq(0, 6, by = 0.05)
  simsi <- modpar[1] + modpar[2]/simres
  plot(si ~ res)
  lines(simsi ~ simres)
  
# alright... so the asymptote is shifted upwards to 0.0028, which is effectively
# zero given our lowest sensitivity index is 0.11
# additionally its predicted with a really poor standard error
  
# doesnt get any better though, when res is equal to 5.2, 6 and 9.9, si is equal
# to 0.11
# -----------------------------------------------------------------------------
# Is different when using linear interpolation as a forcing function
# This forcing function is different in that it uses the "Table" formula type
# Ready workspace
  rm(list=ls(all=TRUE))

# Read in data from simulation resolution test
  res <- c(0.33, 0.65, 1.3, 2.6, 5.2)
  si <- c(1.51, 0.94, 0.48, 0.24, 0.12)
  cmaxsi <- c(0.59, 0.62, 0.19, 0.19, 0.07)
  
  matres <- matrix(1/res, nrow = length(res))
  ones <- matrix(rep(1, times = length(res)), nrow = length(res))
  X <- cbind(ones, matres)
  mod <- lm(si ~ 0 + X)
  summary(mod)
  
# Worse r2 and a larger shift of the asymptote
# lets plot that
  modpar <- mod$coefficients
  simres <- seq(0, 6, by = 0.05)
  simsi <- modpar[1] + modpar[2]/simres
  plot(si ~ res)
  lines(simsi ~ simres)
  
# Cmax is weird because cmax does not occur at the first possible simulated time
# it actually occurs at the closest time to 0.03 hours (the first time point)
# This may also be what's affecting AUC!
# -----------------------------------------------------------------------------
# Lets look at the heart now
# Ready workspace
  rm(list=ls(all=TRUE))

# Read in data from simulation resolution test
  res <- c(0.33, 0.65, 1.3, 2.6, 5.2)
  si <- c(3.72, 1.88, 0.94, 0.47, 0.24)
  
  matres <- matrix(1/res, nrow = length(res))
  ones <- matrix(rep(1, times = length(res)), nrow = length(res))
  X <- cbind(ones, matres)
  mod <- lm(si ~ 0 + X)
  summary(mod)
  
  modpar <- mod$coefficients
  simres <- seq(0, 6, by = 0.05)
  simsi <- modpar[1] + modpar[2]/simres
  plot(si ~ res)
  lines(simsi ~ simres)
  
# -----------------------------------------------------------------------------
# Now kidneys
# Ready workspace
  rm(list=ls(all=TRUE))

# Read in data from simulation resolution test
# 0.33 removed as concentrations are too low!
  res <- c(0.65, 1.3, 2.6, 5.2)
  si <- c(9.08, 4.39, 2.19, 1.09)
  
  matres <- matrix(1/res, nrow = length(res))
  ones <- matrix(rep(1, times = length(res)), nrow = length(res))
  X <- cbind(ones, matres)
  mod <- lm(si ~ 0 + X)
  summary(mod)
  
  modpar <- mod$coefficients
  simres <- seq(0, 6, by = 0.05)
  simsi <- modpar[1] + modpar[2]/simres
  plot(si ~ res)
  lines(simsi ~ simres)
  
# -----------------------------------------------------------------------------
# And lungs
# Ready workspace
  rm(list=ls(all=TRUE))

# Read in data from simulation resolution test
# 0.33 and 0.65 removed as concentrations are too low!
# new resolution of 3.9 added!
  res <- c(1.3, 2.6, 3.9, 5.2)
  si <- c(8.26, 4.01, 2.67, 2.00)
  
  matres <- matrix(1/res, nrow = length(res))
  ones <- matrix(rep(1, times = length(res)), nrow = length(res))
  X <- cbind(ones, matres)
  mod <- lm(si ~ 0 + X)
  summary(mod)
  
  modpar <- mod$coefficients
  simres <- seq(0, 6, by = 0.05)
  simsi <- modpar[1] + modpar[2]/simres
  plot(si ~ res)
  lines(simsi ~ simres)
  
# -----------------------------------------------------------------------------
# Lets see if the same applies in humans!
# Brain first
# Ready workspace
  rm(list=ls(all=TRUE))

# Read in data from simulation resolution test
  res <- c(0.33, 0.65, 1.3, 2.6, 5.2)
  si <- c(5.23, 2.63, 1.32, 0.66, 0.33)
  
  matres <- matrix(1/res, nrow = length(res))
  ones <- matrix(rep(1, times = length(res)), nrow = length(res))
  X <- cbind(ones, matres)
  mod <- lm(si ~ 0 + X)
  summary(mod)
  
  modpar <- mod$coefficients
  simres <- seq(0, 6, by = 0.05)
  simsi <- modpar[1] + modpar[2]/simres
  plot(si ~ res)
  lines(simsi ~ simres)
  