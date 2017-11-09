# Fitting to single tissues
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/WINDOWS/system32",
            "C:/Users/hugjh001/Desktop", "C:/Users/hugjh001/Documents/len_pbpk")
    
    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2] | getwd() == wd[4]) {
      git.dir <- "C:/Users/hugjh001/Documents"
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Load additional libraries
  library(mrgsolve)
  library(splines)
  library(ggplot2)
  library(GA)

# Source functions, data and models
  source(paste(git.dir, reponame, "functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))
  # loads dataiv & dataiv.av
  source(paste(git.dir, reponame, "models", "single_tissues",
    "flow_limited_sumexp.R", sep = "/"))
  source(paste(git.dir, reponame, "models", "single_tissues",
    "memb_limited_sumexp.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "firstorderfit_functions.R", sep = "/"))

# -----------------------------------------------------------------------------
# Add ID numbers to dataiv.av
  dataiv.av$ID <- factor(dataiv.av$DOSEMGKG)
  levels(dataiv.av$ID) <- 1:4
  dataiv.av$ID <- as.numeric(as.character(dataiv.av$ID))

# Create input data.frame but also...
# Determine coefficients for linear forcing function
  input_sumexp <- ddply(dataiv.av, .(ID), function(x) {
    time <- c(0, x$TIME)
    pla <- c(0, x$PLA)
    tis <- c(NaN, x$HRT)
    res <- optim.sumexp(data.frame(time, pla), TRUE, 1)$par[[1]]
    cbind(round(time, 2), signif(res[1], 2), 
      signif(res[2], 2), signif(res[3], 2), round(tis, 4)
    )
  })
  names(input_sumexp)[-1] <- c("time", "M1", "M2", "B", "dv")

# Ready data.frame for simulating data
  ID <- 1:length(unique(dataiv.av$ID))
  time_samp <- c(seq(1, 30, by = 0.5), seq(32, 90, by = 1), seq(95, 480, by = 2.5))
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

  input_simdata <- merge(input_simdata, input_sumexp, all.x = T, all.y = T)
  input_simdata$cmt <- 1
  input_simdata$M1 <- locf(input_simdata$M1)
  input_simdata$M2 <- locf(input_simdata$M2)
  input_simdata$B <- locf(input_simdata$B)
  # input_simdata$V <- 0.01
  # input_simdata$V1 <- 0.015
  # input_simdata$V2 <- 0.008
  # input_simdata$PS <- 0.05

  simdata <- as.data.frame(mrgsim(
    data_set(membmod, input_simdata)
  ))  # mrgsim

  p <- NULL
  p <- ggplot(data = simdata)
  p <- p + geom_line(aes(x = time, y = Cart),
    size = 1, alpha = 0.5, colour = "red")
  p <- p + geom_line(aes(x = time, y = Cven),
    size = 1, alpha = 0.5, colour = "blue")
  p <- p + geom_line(aes(x = time, y = Ctis),
    size = 1, alpha = 0.5, colour = "green4")
  p <- p + geom_point(aes(x = time, y = dv), data = input_sumexp)
  # p <- p + scale_y_log10()
  p <- p + scale_x_continuous(lim = c(0, 50))
  p <- p + facet_wrap(~ID, scales = "free")
  p
