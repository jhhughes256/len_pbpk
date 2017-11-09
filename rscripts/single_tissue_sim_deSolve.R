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
  library(deSolve)
  library(ggplot2)

# Source functions, data and models
  source(paste(git.dir, reponame, "functions",
    "utility.R", sep = "/"))
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set up model functions
  DESflow <- function(T, C, THETAin) {
    dC <- vector(len = 1)
    Cart <- Cart_int(T)
    dC[1] <- Q*(Cart-C[1])/V1
    
    list(dC, "Cart" = Cart)
  }
  
  DESmemb <- function(T, C, THETAin) {
    dC <- vector(len = 2)
    
    Cart <- Cart_int(T)
    dC[1] <- (Q*(Cart  -C[1]) +PS*(C[2] -C[1]))/V1  
    dC[2] <- PS*(C[1] -C[2])/V2  
    
    list(dC, "Cart" = Cart)
  }

# Add ID numbers to dataiv.av
  dataiv.av$ID <- factor(dataiv.av$DOSEMGKG)
  levels(dataiv.av$ID) <- 1:4
  dataiv.av$ID <- as.numeric(as.character(dataiv.av$ID))
  
# Create parameter list
  init_par <- c("V1" = 0.125)
  fixed_par <- c("Q" = 0.9227, "V" = 0.125)
  paramlist <- list(Q = 0.9227, V1 = 0.125)
  C_0 <- c(Cven = 0)
  simdata <- ddply(dataiv.av, .(ID), function(df) {
    Cart_int <<- approxfun(df$TIME, df$PLA, method = "linear", rule = 2)
    simdata <- lsoda(C_0, c(0, df$TIME), DESflow, paramlist)
  })
  Vrat <- paramlist["V"]
  simdata$Ctis <- 
  
# Plot data
  p <- NULL
  p <- ggplot(data = simdata)
  p <- p + geom_line(aes(x = time, y = Cart),
    size = 1, alpha = 0.5, colour = "red")
  p <- p + geom_line(aes(x = time, y = Cven),
    size = 1, alpha = 0.5, colour = "blue")
  p <- p + geom_line(aes(x = time, y = Ctis),
    size = 1, alpha = 0.5, colour = "green4")
  # p <- p + scale_y_log10()
  p <- p + scale_x_continuous(lim = c(0, 50))
  p <- p + facet_wrap(~ID, scales = "free")
  p
  
  