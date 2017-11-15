# Fitting to single tissues
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/WINDOWS/system32",
            "C:/Users/hugjh001/Desktop", "C:/Users/hugjh001/Documents/len_pbpk")
    
    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
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

# Source model
 source(paste(git.dir, reponame, "rscripts", "degradation-app",
    "model.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set up simulation dataset
  ID <- 1
  time_seq <- seq(2, 480, by = 2)
  ID_samp <- rep(ID, each = length(time_seq))
  time_samp <- rep(time_seq, times = length(ID))
  input_simdata <- data.frame(
    ID = ID_samp,
    time = time_samp,
    amt = 0,
    cmt = 1,
    evid = 0
  )
  
  input_amtdata <- data.frame(
    ID = unique(ID),
    time = 0,
    amt = 10,
    cmt = 1,
    evid = 1
  )
  
  input_simdata <- merge(input_simdata, input_amtdata, all = T)
  
  simdata <- as.data.frame(mrgsim(
    data_set(degmod, input_simdata)
  ))  # mrgsim
  
  simdata$Ctot <- simdata$Cpu + simdata$Cpb
  simdata$fu <- with(simdata, 1-(Ctot - Cdu)/Ctot)
  
# 
  p <- NULL
  p <- ggplot(data = simdata)
  p <- p + geom_line(aes(x = simdata$time, y = fu.val))
  
  p <- NULL
  p <- ggplot(data = simdata)
  p <- p + geom_line(aes(x = time, y = Cpu),
    size = 1, alpha = 1, colour = "red")
  # p <- p + geom_line(aes(x = time, y = Cpd),
  #   size = 1, alpha = 0.3, colour = "red")
  p <- p + geom_line(aes(x = time, y = Cdu),
    size = 1, alpha = 1, colour = "blue")
  # p <- p + geom_line(aes(x = time, y = Cdd),
  #   size = 1, alpha = 0.3, colour = "blue")
  p <- p + geom_line(aes(x = time, y = Cpb),
    size = 1, alpha = 0.3, colour = "green4")
  p <- p + geom_line(aes(x = time, y = Ctot),
    size = 1, alpha = 1, colour = "orange")
  # p <- p + scale_y_log10()
  # p <- p + facet_wrap(~ID, scales = "free")
  p