  fig4c <- function(k) {
  # Define simulation resolution and forcing functions
    tal <- 5.2
    forcingFunction <- function(t) {
      0.015*(exp(-3.81E-7*t - 1.34) + exp(-0.0145*t + 3.02) + exp(-0.101*t + 5.99))
    }
  # Create base dataset times
    df <- data.frame(TIME = seq(0, 4/tal, by = 1/tal))
      
  # Define starting concentrations
  # For normal forcing function
    df$FFlo <- forcingFunction(df$TIME)
    df$FFlolag <- forcingFunction(df$TIME + 1/tal)
    
  # For decay-corrected forcing function
    df$FFhi <- df$FFlo/exp(-k/tal)  
    df$FFhilag <- df$FFlolag/exp(-k/tal) 
    
  # For mean method
    df$FFav <- (df$FFlo + df$FFhi)/2
    df$FFavlag <- (df$FFlolag + df$FFhilag)/2 
    
  # For geometric mean method
    df$FFgm <- exp((log(df$FFlo) + log(df$FFhi))/2) 
    df$FFgmlag <- exp((log(df$FFlolag) + log(df$FFhilag))/2)  
    
  # For mixed method
  # Mixed method finds the mean of normal w/ no delay and corrected w/ delay
    df$FFnd <- forcingFunction(df$TIME - 1/tal) 
    df$FFndlag <- df$FFlo
    df$FFmx <- (df$FFnd + df$FFhi)/2       
    df$FFmxlag <- (df$FFndlag + df$FFhilag)/2
    
  # Now we define our reported concentrations
  # For normal and decay-corrected forcing function
    df$PREDhi <- df$FFlo
    df$PREDhi[df$TIME == 0] <- 0
    df$PREDlo <- df$PREDhi*exp(-k/tal)
    
  # For mean and geomean methods
    df$PREDav <- (df$PREDlo + df$PREDhi)/2 # reported conc. for average
    df$PREDgm <- exp((log(df$PREDlo) + log(df$PREDhi))/2)
    
  # For mixed method
    df$PREDnd <- df$FFnd*exp(-k/tal)     # reported conc. for uncorrected
    df$PREDnd[df$TIME == 0] <- 0
    df$PREDmx <- (df$PREDnd + df$PREDhi)/2 # reported conc. for average
    
  # Using ddply create exponential decay
  # Do not use the final row, as we only want 10 decays
    plotdf <- ddply(df[-dim(df)[1],], .(TIME), function(dv) {
      t <- dv$TIME
    # use a really small time to represent assignment occurring after reporting
      tdecay <- c(0, 1e-9, seq(0.1/tal, 0.9/tal, by = 0.1/tal))
      out <- data.frame(
        TIME = t + tdecay 
      )
      out$DV <- forcingFunction(out$TIME)
    # do not use first row, as this is defined by previous decay
      out$PREDlo <- c(dv$PREDlo, dv$FFlolag*exp(-k*tdecay[-1]))
      out$PREDhi <- c(dv$PREDhi, dv$FFhilag*exp(-k*tdecay[-1]))
      out$PREDav <- c(dv$PREDav, dv$FFavlag*exp(-k*tdecay[-1]))
      out$PREDgm <- c(dv$PREDgm, dv$FFgmlag*exp(-k*tdecay[-1]))
      out$PREDnd <- c(dv$PREDnd, dv$FFndlag*exp(-k*tdecay[-1]))
      out$PREDmx <- c(dv$PREDmx, dv$FFmxlag*exp(-k*tdecay[-1]))
      out
    })
    
  # Need to bind the final values to the end of the dataset
  # PREDhi should equal DV at the final time
    lastrow <- tail(df, 1)[, !str_detect(names(df), "FF")]
    lastrow$DV <- lastrow$PREDhi
    plotdf <- rbind(plotdf, lastrow)
  
  # Create mean concentration dataset
  # create for loop to look at each group of 10 numbers from tdecay = 1e-9 to tal
  # +2 is + 1 + 1, one is for indexing starting at 1, other is to skip TIME == 0\
    meandf <- data.frame(NULL)
    for (i in 0:3*11+2) {
      loopdf <- plotdf[i:(i+10),]
      looptmp <- data.frame(
        TIME = c(head(loopdf$TIME, 1), tail(loopdf$TIME, 1)),
        mDV = mean(loopdf$DV),
        mlo = mean(loopdf$PREDlo),
        mhi = mean(loopdf$PREDhi),
        mav = mean(loopdf$PREDav),
        mgm = mean(loopdf$PREDgm),
        mnd = mean(loopdf$PREDnd),
        mmx = mean(loopdf$PREDmx),
        gmDV = exp(mean(log(loopdf$DV))),
        gmlo = exp(mean(log(loopdf$PREDlo))),
        gmhi = exp(mean(log(loopdf$PREDhi))),
        gmav = exp(mean(log(loopdf$PREDav))),
        gmgm = exp(mean(log(loopdf$PREDgm))),
        gmnd = exp(mean(log(loopdf$PREDnd))),
        gmmx = exp(mean(log(loopdf$PREDmx)))
      )
      meandf <- rbind(meandf, looptmp)
    }
    return(plotdf)
  }
  
# Run function on different elimination rate constants
  plotbrain <- fig4c(0.00017816/0.00031271)
  plotbrain$TISSUE <- "Brain"
  # plotheart <- fig4c(0.00038372/0.00031271)
  # plotheart$TISSUE <- "Heart"
  # plotkidney <- fig4c(0.0017816/0.00031271)
  # plotkidney$TISSUE <- "Kidney"
  plotlung <- fig4c(0.0075004/0.00071923)
  plotlung$TISSUE <- "Lung"

# Rbind the plots together
  plotdf2 <- rbind(plotbrain, plotlung)
  