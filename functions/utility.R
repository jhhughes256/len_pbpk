# Process NONMEM SIM data
  processSIMdata <- function(model.path.file.ext) {
  # Work out the nonmem model dir and model file name
    model.path.file <- gsub(".ctl", "", model.path.file.ext)
    model.path <- dirname(model.path.file)
    model.file.ext <- basename(model.path.file.ext)
    model.file <- gsub(".ctl", "", model.file.ext)

    model.dir <- paste(master.dir, "/", model.path, sep = "")

  # Work out the nonmem output dir
    output.dir <- paste(master.dir, "/", model.path.file, dir.extension, sep = "")
    setwd(output.dir)

  # Work out the SIM fit filename
    SIM.file.ext <- paste(model.file, ".fit", sep = "")
  # Work out the SIM fit filename
    SIM.file.ext.out <- paste(model.file, ".fit.csv", sep = "")

  # Need to remove lines with "TABLE NO. 1 in them" and the header lines for each new subject
  # Strip the unwanted lines
    indata <- readLines(SIM.file.ext)
    tablelines <- grep("TABLE NO.  1", indata)  #May be installation specific
    headerlines <- grep(" ID", indata)          #May be installation specific

  # Extract the information in the header line for column names
    header <- indata[headerlines[1]]
    header <- scan(textConnection(header), what = "char")
    colnum <- length(header)

  # Strip out the unwanted lines
    badlines <- c(tablelines,headerlines)
    indata <- indata[-badlines]

  # replace white space with commas
    for (i in 1:length(indata)) {
      indata[i] <- gsub('[[:space:]]+', ',', indata[i])
    }
    
  # write to a file
    writeLines(indata, "SIMtemp.txt")

  # read again as csv file
    SIMdata <- read.csv("SIMtemp.txt", header = F)
    SIMdata <- SIMdata[, -1]     #delete first blank column
    names(SIMdata) <- header

  # The NONMEM output does not make a new ID number for each simulation
  # Therefore make a list of SIM numbers
    nsims <- length(tablelines)
    numtimepoints <- length(SIMdata$ID)
    numtimespersubject <- numtimepoints/nsims
    SIMdata$SIM <- rep(1:nsims, each = numtimespersubject)

  # Write the SIMdata to a file
    write.csv(SIMdata, SIM.file.ext.out, row.names = F)

  # tidy up
  # consider deleting the original fit file?
    file.remove("SIMtemp.txt")
    setwd(master.dir)
}

#--------------------------------------------------------------------------------------------
# DATA SUMMARY FUNCTIONS

  oneperID <- function(x) {
    # returns a single value for each ID of dataframe
    ID <- head(x, n = 1)
    ID
  }

# 90% prediction intervals
  CI80lo <- function(x) quantile(x, probs = 0.1)
  CI80hi <- function(x) quantile(x, probs = 0.9)
# 95% prediction intervals
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)

# Define a function for geometric mean
  geomean <- function(x, na.rm = F) {
    if (na.rm == T) x <- x[is.na(x) == F]
    exp(mean(log(x)))
  }
# Note x cannot be negative, zero

#Median and 90% tolerance intervals
  sumfunc90 <- function(x)
  {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  quantile(x, probs = 0.05, na.rm = T, names = F)  #90%CI
    stat3 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    stat4 <-  lengthNA(x)
    result <- c("median" = stat1, "low90" = stat2, "hi90" = stat3, "n" = stat4)
    result
  }

 #Median and 95% tolerance intervals
  sumfunc95 <- function(x)
  {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  quantile(x, probs = 0.025, na.rm = T, names = F)  #95%CI
    stat3 <-  quantile(x, probs = 0.975, na.rm = T, names = F)
    stat4 <-  lengthNA(x)
    result <- c("median" = stat1, "low95" = stat2, "hi95" = stat3, "n" = stat4)
    result
  }

#Define a function for geometric mean and 90% CI of the sem
  geomeansemCI <- function(x, na.rm = F) {
  # Note x cannot be negative, zero
    logx <- log(x)
    logmean <- mean(logx)
    n <- length(x)
    logsem <- sd(logx)/sqrt(n)
  # Critical value of the t-distribution for two one-sided p=0.05
    critt <- qt(.95, df = (n - 1))
    loglo95 <- logmean - critt*logsem
    loghi95 <- logmean + critt*logsem
    gmean <- exp(logmean)
    glo95 <- exp(loglo95)
    ghi95 <- exp(loghi95)
    result <- c(
      "gmean" = gmean, "glo95" = glo95, "ghi95" = ghi95, "crit.t" = critt
    )  # result
    result
  }

### General R Functions
# Allows combination of two or more vectors such that the first occurrence of a
# non NA is taken
# > a <- c(1,  2,  NA, NA, NA)
# > b <- c(NA, NA, NA, NA, 6)
# > c <- c(7,  8,  NA, 9, 10)
# > coalesce2(a, b, c)
# # [1]  1  2 NA  9  6
  coalesce2 <- function(...) {
    Reduce(
      function(x, y) {
        i <- which(is.na(x))
        x[i] <- y[i]
        x
      },
      list(...)
    )  # Reduce
  }

  gsub.all <- function(pat,repl,x) {
    vec <- x
    if (length(pat) == length(repl)) {
      for(i in 1:length(pat)) {
	      vec <- gsub(pat[i], repl[i], vec, fixed = TRUE)
	    }
	    vec
	  } else {
	    warning("length(pattern)!=length(replacement): Amount of values to be
        replaced is not equal to amount of values given")
	  }
  }

# -----------------------------------------------------------------------------
# Function for filling in gaps in covariate data for one ID - use ddply for all
# Only applies for covariates that don't change with time
# test0 <- c(0.05,0.05,0.05,0.05,0.05,0.05,0.05)
# test1 <- c(NA,0.05,0.05,NA,0.05,NA,0.05)
# test2 <- c(NA,"2_0.0833_0.05_1","2_0.0833_0.05_1",NA,"2_0.0833_0.05_1",NA,"2_0.0833_0.05_1")
# test3 <- c(NA,NA,NA,NA)
# test4 <- c(0.05,0.01,0.05,0.05,NA,0.05,0.05)

  fill.missing <- function(x) {
    all.fill <- unique(x[is.na(x)==F])  #All non-missing values
    fill <- all.fill[1]  #value to replace missing values
    if (length(all.fill) > 1) {
      stop("Can't fill with non-unique values")
    } else {
      x[is.na(x) == T] <- fill
    }
    x
  }

# fill.missing(test0)
# fill.missing(test1)
# fill.missing(test2)
# fill.missing(test3)
# fill.missing(test4)

#----------------------------------------------------------------------------------------------

###Imputation functions###

locf <- function (x)
#Last observation carried forward
#Finds an NA and carries forward the previous value
  {
    good <- !is.na(x)
    positions <- seq(length(x))
    good.positions <- good * positions
    last.good.position <- cummax(good.positions)
    last.good.position[last.good.position == 0] <- NA
    x[last.good.position]
}


nocb <- function (x)
#Next observation carried backward
#Reverses the locf function
  {
   rev(locf(rev(x)))
  }


impute <- function(x)
#Function that runs locf first, then nocb
{
    x <- locf(x)
    nocb(x)
    x
}

#-------------------------------------------------------------------------------------------
###Functions for counting missing data###

#see also lengthNA

#Function for calculating percent missing - use apply to do column by column
  calculate.percent.missing <- function(x) {
    length(x[is.na(x)==T])/length(x)*100
  }

# #Testing
# Concentration <- c(1,2,3,4)
# percent.missing(Concentration)

# Concentration <- c(1,2,3,NA,4)
# percent.missing(Concentration)

# Concentration <- c(NA,NA,NA,NA,NA)
# percent.missing(Concentration)
