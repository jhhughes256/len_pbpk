# Fitting to single tissues
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Source functions, data and models
  source(paste(git.dir, reponame, "functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))
  # loads dataiv & dataiv.av
  source(paste(git.dir, reponame, "models", "single_tissues",
    "flow_limited.R", sep = "/"))
  source(paste(git.dir, reponame, "models", "single_tissues",
    "memb_limited.R", sep = "/"))

# Load additional libraries
  library(mrgsolve)
  library(splines)

# -----------------------------------------------------------------------------
# Add ID numbers to dataiv.av
  dataiv.av$ID <- factor(dataiv.av$DOSEMGKG)
  levels(dataiv.av$ID) <- 1:4
  dataiv.av$ID <- as.numeric(as.character(dataiv.av$ID))

# Bring together artery concentration data
# Units for concentrations are in ng/mL
  Cart_colnames <- c("ID", "DOSEMGKG", "TADNOM", "DOSEMG", "AMT", "WT", "TIME", "LUN")
  Cart_data <- dataiv.av[Cart_colnames]
  names(Cart_data)[length(Cart_colnames)] <- "DV"

  # Cart_data$DVNORM <- with(Cart_data, DV/AMT)

# Determine coefficients for cubic spline forcing function
  Cart_cub <- ddply(Cart_data, .(ID), function(x) {
    time <- c(0, x$TIME)
    dv <- c(0, x$DV)
    spline <- interpSpline(dv ~ time)
    cbind(signif(spline$knots, 2), signif(spline$coefficients, 2))
  })
  names(Cart_cub)[-1] <- c("CUBT", "COF1", "COF2", "COF3", "COF4")

# Ready data.frame for simulating data
  ID <- 1:length(unique(dataiv.av$ID))
  time_samp <- c(seq(0, 30, by = 2), seq(35, 90, by = 5), seq(120, 480, by = 30))
  ID2 <- sort(c(rep(ID, times = length(time_samp))))
  times <- rep(time_samp, times = length(ID))
  input_simdata <- data.frame(
    ID = ID2,
    time = times
  )

# Ensure that doses 0.5 and 1.5 don't go over 300 mins
  input_simdata <- input_simdata[
    input_simdata$ID < 3 & input_simdata$time <= 300 |
    input_simdata$ID > 2,
  ]

  input_cubdata <- Cart_cub
  input_cubdata$time <- Cart_cub$CUBT

  input_simdata <- merge(input_simdata, input_cubdata, all.x = T, all.y = T)
  input_simdata$cmt <- 1
  input_simdata$CUBT <- locf(input_simdata$CUBT)
  input_simdata$COF1 <- locf(input_simdata$COF1)
  input_simdata$COF2 <- locf(input_simdata$COF2)
  input_simdata$COF3 <- locf(input_simdata$COF3)
  input_simdata$COF4 <- locf(input_simdata$COF4)

  simdata <- as.data.frame(mrgsim(
    data_set(flowmod, input_simdata)
  ))  # mrgsim
