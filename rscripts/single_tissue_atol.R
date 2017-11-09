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
  
# Source functions, data and models
  source(paste(git.dir, reponame, "functions",
               "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
               "data_iv.R", sep = "/"))
# loads dataiv & dataiv.av
  source(paste(git.dir, reponame, "models", "single_tissues",
               "flow_limited_linear.R", sep = "/"))
  source(paste(git.dir, reponame, "models", "single_tissues",
               "memb_limited_linear.R", sep = "/"))
  source(paste(git.dir, reponame, "models", "single_tissues",
               "linear.R", sep = "/"))

# -----------------------------------------------------------------------------
# Add ID numbers to dataiv.av
  dataiv.av$ID <- factor(dataiv.av$DOSEMGKG)
  levels(dataiv.av$ID) <- 1:4
  dataiv.av$ID <- as.numeric(as.character(dataiv.av$ID))

# Bring together artery concentration data
# Units for concentrations are in ng/mL
  Cart_colnames <- c("ID", "DOSEMGKG", "TADNOM", "DOSEMG", "AMT", "WT", "TIME", "PLA")
  Cart_data <- dataiv.av[Cart_colnames]
  names(Cart_data)[length(Cart_colnames)] <- "DV"

# Cart_data$DVNORM <- with(Cart_data, DV/AMT)

# Determine coefficients for linear forcing function
  input_lindata <- ddply(Cart_data, .(ID), function(x) {
    time <- c(0, x$TIME)
    dv <- c(0, x$DV)
    last <- length(dv)
    slope <- diff(dv)/diff(time)
    int <- dv[-last]-slope*time[-last]
    cbind(signif(time[-last], 2), signif(slope, 2), signif(int, 2))
  })
  names(input_lindata)[-1] <- c("time", "M", "B")

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
  
  input_simdata <- merge(input_simdata, input_lindata, all.x = T, all.y = T)
  input_simdata$cmt <- 1
  input_simdata$M <- locf(input_simdata$M)
  input_simdata$B <- locf(input_simdata$B)
  
# Determine correct value for forcing function
  forcedata <- ddply(input_simdata, .(ID), function(df) {
    data.frame(time = df$time, Cart = df$M*df$time + df$B)
  })
  
# Determine required atol or rtol for linmod
  lindata <- as.data.frame(mrgsim(
    data_set(linmod, input_simdata), rtol = 1e-12
  ))  # mrgsim
  
  flowdata <- as.data.frame(mrgsim(
    data_set(flowmod, input_simdata), rtol = 1e-12
  ))  # mrgsim
  
  membdata <- as.data.frame(mrgsim(
    data_set(membmod, input_simdata), rtol = 1e-12
  ))  # mrgsim
  
  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = time, y = Cart), data = lindata,
                     size = 1, alpha = 0.5, colour = "red")
  p <- p + geom_line(aes(x = time, y = Cart), data = forcedata,
                     size = 1, alpha = 0.5, colour = "black")
  p <- p + geom_line(aes(x = time, y = Cart), data = flowdata,
                     size = 1, alpha = 0.5, colour = "blue")
  p <- p + geom_line(aes(x = time, y = Cart), data = membdata,
                     size = 1, alpha = 0.5, colour = "green4")
  p <- p + scale_y_log10()
  p <- p + scale_x_continuous(lim = c(0, 50))
  p <- p + facet_wrap(~ID, scales = "free")
  p
