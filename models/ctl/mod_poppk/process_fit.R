# R script for processing NONMEM .fit files
# Specifically for mouse PBPK data
# -----------------------------------------------------------------------------
# Ready workspace
# -----------------------------------------------------------------------------
# Identify if git dir. exists and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list=ls(all=T))
	graphics.off()
	git.dir <- "J:/Git"
	#git.dir <- "E:/Hughes/Git"
	#git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
	#git.dir <- "C:/Users/hugjh001/Documents"
	reponame <- "len_pbpk"
	scriptname <- "process_fit.R"
  }
  
# Set working directory
  master.dir <- paste(git.dir, reponame, sep = "/")
  setwd(master.dir)
  
# Load library packages
  library(ggplot2)
  
# Source observed data
  nmprep <- read.csv("J:/NM_POPPK/nmprep.csv")
  
 #Load libraries
  library(ggplot2)
  library(doBy)
  library(plyr)
  library(stringr)
  library(R2HTML)
  library(Hmisc)
  library(grid)
  library(reshape2)
  library(GGally)

# -----------------------------------------------------------------------------
# Read in fit file
# -----------------------------------------------------------------------------
# Set the name of the required file and set the working directory
  cat("Select one of files in directory to process:\n")
  path <- gsub("\\\\", "/", file.choose())
  base.path <- dirname(path)
  setwd(base.path)
  file.name.in <- basename(path)
  file.name.out <- paste(file.name.in, ".csv", sep = "")
  runfolder <- base.path  # picks up the folder of the run being analysed

# Read *.fit file and attach, so column names are available
  fitdata <- read.table(file=file.name.in, sep = "", 
	skip = 1, header = T, na.strings = c("NA","***********","1.#INFE+00")
  )  # read.table
  
# Write to file
  # write.csv(fitdata, file = file.name.out)

# Remove dose events & missing values
  fitdata <- fitdata[fitdata$MDV != 1, ]
  
# Set factors and levels
# Note: not all categorical covariates have been included
  fitdata$IDf <- as.factor(fitdata$ID)

  fitdata$TADNOMf <- as.factor(fitdata$TADNOM)
  
  fitdata$DOSEMGKGf <- as.factor(fitdata$DOSEMGKG)
  levels(fitdata$DOSEMGKGf) <- paste(levels(fitdata$DOSEMGKGf), "mg/kg")
  
  fitdata$ROUTEf <- as.factor(fitdata$ROUTE)
  levels(fitdata$ROUTEf) <- c("IV", "PO")

# -----------------------------------------------------------------------------
# Diagnostic plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Goodness of Fit Plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Combine four ggplot2 graphs in a grid layout
  filename <- "diagnostic_dashboard.png"
  
### DEVICE ON ###
  png(filename, width = 720, height = 920)  # pixels = 2.4 in x 3.1 in @ 300dpi
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(4, 4)))

#  Plot 1  PRED vs DV
  max.OBS1 <- max(c(fitdata$PRED, fitdata$DV), na.rm = T)
  
  plotobj1 <- NULL
  plotobj1 <- ggplot(fitdata)
  plotobj1 <- plotobj1 + geom_point(aes(x = PRED, y = DV), shape = 1)
  plotobj1 <- plotobj1 + geom_abline(aes(x = PRED, y = DV), 
    intercept = 0, slope = 1, colour = "black")  #Add line of identity
  plotobj1 <- plotobj1 + geom_smooth(aes(x = PRED, y = DV), 
    method = loess, se = T, colour = "red")  #Add loess smoothing line
  plotobj1 <- plotobj1 + scale_x_continuous("Population Predicted conc (ug/mL)",
    lim = c(-0.1, max.OBS1))
  plotobj1 <- plotobj1 + scale_y_continuous("Observed conc (ug/mL)",
    lim = c(-0.1, max.OBS1))
  plotobj1 <- plotobj1 + scale_colour_brewer("Dose Level", palette = "Set1")
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  print(plotobj1, vp=vplayout(1:2, 1:2))

#  Plot 2  IPRE vs DV
  max.OBS2 <- max(c(fitdata$IPRED, fitdata$DV), na.rm = T)

  plotobj2 <- NULL
  plotobj2 <- ggplot(fitdata)
  plotobj2 <- plotobj2 + geom_point(aes(x = IPRE, y = DV), shape = 1)
  plotobj2 <- plotobj2 + geom_abline(aes(x = IPRE, y = DV), 
    intercept = 0, slope = 1, colour = "black")  #Add line of identity
  plotobj2 <- plotobj2 + geom_smooth(aes(x = IPRE, y = DV), 
    method = loess, se = T, colour = "red")  #Add loess smoothing line
  plotobj2 <- plotobj2 + geom_abline(aes(x = IPRE, y = DV), 
    intercept = 0, slope = 2, colour = "darkgreen", linetype = "dashed")
  plotobj2 <- plotobj2 + geom_abline(aes(x = IPRE, y = DV), 
    intercept = 0, slope = 0.5, colour = "darkgreen", linetype = "dashed")
  plotobj2 <- plotobj2 + scale_x_continuous("IPRED conc (ug/mL)", 
    lim = c(-0.1,max.OBS2))
  plotobj2 <- plotobj2 + scale_y_continuous("Observed conc (ug/mL)", 
    lim = c(-0.1,max.OBS2))
  plotobj2 <- plotobj2 + scale_colour_brewer("Dose Level", palette = "Set1")
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  print(plotobj2, vp = vplayout(1:2, 3:4))

#  Plot 3  CWRES vs TAFDE
  max.CWRES <- max(abs(fitdata$CWRES), na.rm = T)
  if (max.CWRES < 0.1) max.CWRES <- 0.1 

  plotobj3 <- NULL
  plotobj3 <-  ggplot(fitdata)
  plotobj3 <- plotobj3 + geom_point(aes(x = TIME, y = CWRES), shape = 1)
  plotobj3 <- plotobj3 + geom_abline(aes(x = TIME, y = CWRES),
    intercept = 0, slope = 0, colour = "black")  #Add zero line
  plotobj3 <- plotobj3 + geom_smooth(aes(x = TIME, y = CWRES),
    method = loess, se = T, colour = "red")  #Add loess smoothing line
  plotobj3 <- plotobj3 + scale_x_continuous("Time after dose (h)")
  plotobj3 <- plotobj3 + scale_y_continuous("CWRES", lim = c(-max.CWRES, max.CWRES))
  plotobj3 <- plotobj3 + scale_colour_brewer("Dose Level", palette = "Set1")
  print(plotobj3, vp = vplayout(3, 1:4))

# Plot 4   CWRES vs PRED
  max.CWRES <- max(abs(fitdata$CWRES), na.rm = T)
  if (max.CWRES < 0.1) max.CWRES <- 0.1

  plotobj4 <- NULL
  plotobj4 <-  ggplot(fitdata)
  plotobj4 <- plotobj4 + geom_point(aes(x = PRED, y = CWRES), shape = 1)
  plotobj4 <- plotobj4 + geom_abline(aes(x = PRED, y = CWRES),
    intercept = 0, slope = 0, colour = "black")  #Add zero line
  plotobj4 <- plotobj4 + geom_smooth(aes(x = PRED, y = CWRES),
    method = loess, se = T, colour = "red")  #Add loess smoothing line
  plotobj4 <- plotobj4 + scale_x_continuous("Population Predicted conc (ug/mL)")
  plotobj4 <- plotobj4 + scale_y_continuous("CWRES", 
    lim = c(-max.CWRES, max.CWRES))
  plotobj4 <- plotobj4 + scale_colour_brewer("Dose Level", palette = "Set1")
  print(plotobj4, vp = vplayout(4, 1:4))

  dev.off()
### DEVICE OFF ###

# Individual diagnostic plots
  ggsave("GoF_DV_PRED.png", plotobj1, width = 2.2, height = 1.7)
  ggsave("GoF_DV_IPRED.png", plotobj2, width = 2.2, height = 1.7)
  ggsave("GoF_CWRES_TAFDE.png", plotobj3, width = 2.2, height = 1.7)
  ggsave("GoF_CWRES_PRED.png", plotobj4, width = 2.2, height = 1.7)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Goodness of Fit Plots - grouped by factor
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  diag.plot <- function(df1, icol, dlow, dup) {
    filename <- paste0("diagnostic_dashboard_", 
	  icol, dlow, ifelse(dlow != dup, dup, ""), ".png")
    col.interest <- which(names(fitdata) %in% icol)
    fitdata <- df1[df1[col.interest] <= dup & df1[col.interest] >= dlow, ]
    png(filename, width = 720, height = 920)  # pixels = 2.4 in x 3.1 in @ 300dpi
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 4)))

  # Plot 1  PRED vs DV
    max.OBS1 <- max(c(fitdata$PRED, fitdata$DV), na.rm = T)
  
    plotobj1 <- NULL
    plotobj1 <- ggplot(fitdata)
    plotobj1 <- plotobj1 + geom_point(aes(x = PRED, y = DV), shape = 1)
    plotobj1 <- plotobj1 + geom_abline(aes(x = PRED, y = DV), 
      intercept = 0, slope = 1, colour = "black")  #Add line of identity
    plotobj1 <- plotobj1 + geom_smooth(aes(x = PRED, y = DV), 
      method = loess, se = T, colour = "red")  #Add loess smoothing line
    plotobj1 <- plotobj1 + scale_x_continuous("Population Predicted conc (ug/mL)",
      lim = c(-0.1, max.OBS1))
    plotobj1 <- plotobj1 + scale_y_continuous("Observed conc (ug/mL)",
      lim = c(-0.1, max.OBS1))
    plotobj1 <- plotobj1 + scale_colour_brewer("Dose Level", palette = "Set1")
    plotobj1 <- plotobj1 + theme(legend.position = "none")
    print(plotobj1, vp=vplayout(1:2, 1:2))

    # Plot 2  IPRE vs DV
    max.OBS2 <- max(c(fitdata$IPRED, fitdata$DV), na.rm = T)

    plotobj2 <- NULL
    plotobj2 <- ggplot(fitdata)
    plotobj2 <- plotobj2 + geom_point(aes(x = IPRE, y = DV), shape = 1)
    plotobj2 <- plotobj2 + geom_abline(aes(x = IPRE, y = DV), 
      intercept = 0, slope = 1, colour = "black")  #Add line of identity
    plotobj2 <- plotobj2 + geom_smooth(aes(x = IPRE, y = DV), 
      method = loess, se = T, colour = "red")  #Add loess smoothing line
    plotobj2 <- plotobj2 + geom_abline(aes(x = IPRE, y = DV), 
      intercept = 0, slope = 2, colour = "darkgreen", linetype = "dashed")
    plotobj2 <- plotobj2 + geom_abline(aes(x = IPRE, y = DV), 
      intercept = 0, slope = 0.5, colour = "darkgreen", linetype = "dashed")
    plotobj2 <- plotobj2 + scale_x_continuous("IPRED conc (ug/mL)", 
      lim = c(-0.1,max.OBS2))
    plotobj2 <- plotobj2 + scale_y_continuous("Observed conc (ug/mL)", 
      lim = c(-0.1,max.OBS2))
    plotobj2 <- plotobj2 + scale_colour_brewer("Dose Level", palette = "Set1")
    plotobj2 <- plotobj2 + theme(legend.position = "none")
    print(plotobj2, vp = vplayout(1:2, 3:4))

  # Plot 3  CWRES vs TAFDE
    max.CWRES <- max(abs(fitdata$CWRES), na.rm = T)
    if (max.CWRES < 0.1) max.CWRES <- 0.1 

    plotobj3 <- NULL
    plotobj3 <-  ggplot(fitdata)
    plotobj3 <- plotobj3 + geom_point(aes(x = TIME, y = CWRES), shape = 1)
    plotobj3 <- plotobj3 + geom_abline(aes(x = TIME, y = CWRES),
      intercept = 0, slope = 0, colour = "black")  #Add zero line
    plotobj3 <- plotobj3 + geom_smooth(aes(x = TIME, y = CWRES),
      method = loess, se = T, colour = "red")  #Add loess smoothing line
    plotobj3 <- plotobj3 + scale_x_continuous("Time after dose (h)")
    plotobj3 <- plotobj3 + scale_y_continuous("CWRES", lim = c(-max.CWRES, max.CWRES))
    plotobj3 <- plotobj3 + scale_colour_brewer("Dose Level", palette = "Set1")
    print(plotobj3, vp = vplayout(3, 1:4))

  # Plot 4   CWRES vs PRED
    max.CWRES <- max(abs(fitdata$CWRES), na.rm = T)
    if (max.CWRES < 0.1) max.CWRES <- 0.1

    plotobj4 <- NULL
    plotobj4 <-  ggplot(fitdata)
    plotobj4 <- plotobj4 + geom_point(aes(x = PRED, y = CWRES), shape = 1)
    plotobj4 <- plotobj4 + geom_abline(aes(x = PRED, y = CWRES),
      intercept = 0, slope = 0, colour = "black")  #Add zero line
    plotobj4 <- plotobj4 + geom_smooth(aes(x = PRED, y = CWRES),
      method = loess, se = T, colour = "red")  #Add loess smoothing line
    plotobj4 <- plotobj4 + scale_x_continuous("Population Predicted conc (ug/mL)")
    plotobj4 <- plotobj4 + scale_y_continuous("CWRES", 
      lim = c(-max.CWRES, max.CWRES))
    plotobj4 <- plotobj4 + scale_colour_brewer("Dose Level", palette = "Set1")
    print(plotobj4, vp = vplayout(4, 1:4))

    dev.off()
  }  # diag.plot
  
  diag.plot(fitdata, "DOSEMGKG", 0.5, 0.5)
  diag.plot(fitdata, "DOSEMGKG", 1.5, 1.5)
  diag.plot(fitdata, "DOSEMGKG", 5, 5)
  diag.plot(fitdata, "DOSEMGKG", 10, 10)
  diag.plot(fitdata, "TADNOM", 2, 20)
  diag.plot(fitdata, "TADNOM", 40, 90)
  diag.plot(fitdata, "TADNOM", 180, 480)
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Residual Error Distribution Plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QQ plot
  plotobj <- NULL
  plotobj <- ggplot(fitdata)
  plotobj <- plotobj + ggtitle("QQ normal plot of residuals\n")
  plotobj <- plotobj + scale_x_continuous("Theoretical", lim = c(-3, 3)) 
  plotobj <- plotobj + scale_y_continuous("Observed", lim = c(-3, 3))
  plotobj <- plotobj + geom_abline(aes(sample = CWRES),
    intercept = 0, slope = 1, colour = "black")  #Add line of identity
	
# Condition qq plot for all data
  plotobj1 <- NULL
  plotobj1 <- plotobj + stat_qq(aes(sample = CWRES), geom = "point")
  #plotobj1
  ggsave("residual_qq.png", plotobj1, width = 2.2, height = 1.7)

# Condition qq plot on dose group
  plotobj2 <- NULL
  plotobj2 <- plotobj + stat_qq(aes(sample = CWRES, colour = DOSEMGKGf),
    geom = "point", alpha = 0.5)
  plotobj2 <- plotobj2 + scale_colour_brewer(name = "Dose Group", palette = "Set1")
  #plotobj2
  ggsave("residual_qq_group.png", plotobj2, width = 2.2, height = 1.7)

# Condition qq plot on nominal time after dose group
  plotobj3 <- NULL
  plotobj3 <- plotobj + stat_qq(aes(sample = CWRES, colour = TADNOMf),
    geom = "point", alpha = 0.5)
  plotobj3 <- plotobj3 + scale_colour_brewer(name = "Nominal TAD", palette = "Set1")
  #plotobj3
  ggsave("residual_qq_study.png", plotobj3, width = 2.2, height = 1.7)


# Residual density plot (instead of histogram)
  plotobj <- NULL
  plotobj <- ggplot(fitdata)
  plotobj <- plotobj + ggtitle("Distribution density of residuals\n")
  plotobj <- plotobj + geom_density(aes(x=CWRES, y=..density..))
  plotobj <- plotobj + scale_x_continuous("Residual", lim = c(-4, 4))
  plotobj <- plotobj + scale_y_continuous("distribution density")
  #plotobj

  ggsave("residual_density.png", plotobj, width = 2.2, height = 1.7)


# Condition density plot on group
  plotobj1 <- plotobj + geom_density(aes(x = CWRES, y = ..density..,
    colour = DOSEMGKGf))
  plotobj1 <- plotobj1 + scale_colour_brewer("Group", palette = "Set1")
  #plotobj1

  ggsave("residual_density_group", plotobj, width = 2.2, height = 1.7)

#--------------------------------------------------------------------------------------------------
#Correlation matrix between eta's
#Subset the eta's from file - make sure they have been tabled in the *.fit file

#ETA correlation--------------------------------------------------------------------------------------------------------------------

#Make a dataframe with only 1 line per subject (avoids duplications of ETA as it is the same for all observations in a subject)
  fitdataone <- ddply(fitdata, .(ID), oneperID)
#Can't use this as there are several dose levels per patient

  eta.cols <-  grep(glob2rx("ETA*"), names(fitdataone))
  etabov.cols <- c("ETA4","ETA5","ETA6","ETA7")  #this is hard-coded for BOV
  #eta.cols <- c("ETA4","ETA5")  #hard-coded#
  eta.cols <- eta.cols[eta.cols %in% etabov.cols==F]
  eta.cols <- names(fitdataone)[eta.cols]
  etadata <- subset(fitdataone, select=eta.cols)

if (ncol(etadata)>1)  #more than 1 ETA is scatterplot matrix
{
  plotobj <- NULL
  plotobj <-  ggpairs(etadata, title ="Correlation between random effects")
  #plotobj
  to.png(plotobj,"etascatter")
} else               #only 1 ETA  is density plot
{
  plotobj <- NULL
  plotobj  <- ggplot(data=etadata)
  plotobj <- plotobj + geom_density(aes(x=ETA1, y=..density..), colour="black")
  plotobj <- plotobj +  scale_x_continuous("ETA1")
  plotobj <- plotobj +  scale_y_continuous("distribution density")
  to.png(plotobj,"etascatter")
}

#--------------------------------------------------------------------------------------------------
#Look at ETA versus covariate relationships
#Make a dataframe with only 1 line per subject (avoids duplications of ETA as it is the same for all observations in a subject)
#uses above fitdataone <- ddply(fitdata, .(ID), oneperID)

#eta.cols <-  grep(glob2rx("ETA*"), names(fitdataone))
#etabov.cols <- c("ETA1","ETA2","ETA3")  #this is hard-coded for BOV
#eta.cols <- c("ETA4","ETA5")  #hard-coded#
#eta.cols <- eta.cols[eta.cols %in% etabov.cols==F]

#Get the columns with categorical covariates
covcat.cols <- c("STUDYf","DOSELVLf","GENDf","RACEf","DXCATf","OCCf")

#Get the columns with continuous covariates
covcont.cols <- c("AGE","WT","HT","SECR","CRCL","IBW")

#-------------------------------------------------------------------------------------------------------
#Categorical covariates
covcatdf <- expand.grid(eta.cols,covcat.cols,stringsAsFactors = F)
names(covcatdf) <- c("ETAname","covname")
#covcatdf


#Function to count numbers in a boxplot and return as a coardinate,label pair
boxplot.give.n <- function(x)
{
  return(c(y = median(x), label = length(x)))
}

#Function to plot effect of categorical covariate on ETA
ETACovariatePlotCAT <- function(ETAname,covname)
{

  plotobj <- NULL
  ETAtext <- ETAname  #This could be made an argument to the function
  covtext <- covname  #This could be made an argument to the function
  plotobj <-  ggplot(data=fitdataone)
  plotobj <- plotobj + geom_boxplot(aes_string(x = covname, y = ETAname), position=position_dodge(width=0.9))
  plotobj <- plotobj + stat_summary(aes_string(x = covname, y = ETAname), fun.data = boxplot.give.n, geom = "text", size = 6, colour="red")
  plotobj <- plotobj + scale_x_discrete(paste(covtext))
  plotobj <- plotobj + scale_y_continuous(paste(ETAtext))
  #plotobj  <- plotobj + ggtitle("Final PK model\n")  #legend.position="none",
  #plotobj <- plotobj + ggtitle("Base PK model\n")  #legend.position="none",
  plotobj <- plotobj + ggtitle("ETA PLOT\n")
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  #plotobj

  png.file.name <- paste(ETAname,"_vs_",covname,sep="")

  to.png.sqr(plotobj,png.file.name)

}

#Apply the plotting function - stand back and watch the magic!
#This runs the ETACovariatePlotCAT plotting function, taking each row of covcatdf as the input to the function
mdply(covcatdf, ETACovariatePlotCAT)

#-------------------------------------------------------------------------------------------------------
#Continuous covariates
#Make a little dataframe with each row being the arguments for our covariate plotting functions
covcontdf <- expand.grid(eta.cols,covcont.cols,stringsAsFactors = F)
names(covcontdf) <- c("ETAname","covname")
covcontdf

#Function to plot effect of continuous covariate on ETA
ETACovariatePlotCONT <- function(ETAname,covname)
{

  plotobj <- NULL
  ETAtext <- ETAname  #This could be made an argument to the function
  covtext <- covname  #This could be made an argument to the function
  plotobj <-  ggplot(data=fitdataone)
  plotobj <- plotobj + geom_point(aes_string(x = covname, y = ETAname), colour="blue", size=2)
  plotobj <- plotobj + geom_smooth(aes_string(x = covname, y = ETAname), method=loess, se=T, linetype=2)
  plotobj <- plotobj + geom_smooth(aes_string(x = covname, y = ETAname), method=lm, se=F, colour="red")
  plotobj <- plotobj + geom_abline(intercept=0, slope=0, colour="black")
  plotobj <- plotobj + scale_x_continuous(paste(covtext))
  plotobj <- plotobj + scale_y_continuous(paste(ETAtext))
  #plotobj <- plotobj + theme(legend.position="none") + ggtitle("Final PK model\n")
  #plotobj <- plotobj + ggtitle("Base PK model\n")
  plotobj <- plotobj + ggtitle("ETA PLOT\n")
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  #plotobj

  png.file.name <- paste(ETAname,"_vs_",covname,sep="")
  to.png.sqr(plotobj,png.file.name)

}

#Apply the plotting function - stand back and watch the magic!
#This runs the ETACovariatePlotCAT plotting function, taking each row of covcontdf as the input to the function
	mdply(covcontdf, ETACovariatePlotCONT)

  param.cols <- c("CL", "V2", "KTR")

  covcatpf <- expand.grid(param.cols,covcat.cols,stringsAsFactors = F)
  names(covcatpf) <- c("ETAname","covname")
	mdply(covcatpf, ETACovariatePlotCAT)

  covcontpf <- expand.grid(param.cols,covcont.cols,stringsAsFactors = F)
  names(covcontpf) <- c("ETAname","covname")
	mdply(covcontpf, ETACovariatePlotCONT)

#Covariate correlation----------------------------------------------------------

	cov.cols <- c(covcat.cols,covcont.cols)
	cov.cols <- c(cov.cols,"CRCL")
	covsubdata <- subset(fitdataone, select=cov.cols)
	plotobj <- NULL
  plotobj <- ggpairs(covsubdata, title ="Correlation between covariates")

	#ggsave("covscatter.png", width=5, height=5, units=c("cm"))

#--------------------------------------------------------------------------------------------------
#Complex ETA grids

	#plotobj1 <- NULL
	#plotobj1 <- ggplot(fitdata)
	#titletext <- expression(atop("ETA1 v's ROUTE",
	#												atop("BY ROUTE",
	#														 "coloured by DOSING")))
	#plotobj1 <- plotobj1 + ggtitle(titletext)
	#plotobj1 <- plotobj1 + geom_point(aes(x=ROUTE, y=ETA1, colour = DOSING), shape=1)
	#plotobj1 <- plotobj1+ scale_x_discrete(name="ROUTE")
	#plotobj1 <- plotobj1+ scale_y_continuous(name="ETA1")
	#plotobj1 <- plotobj1 + facet_wrap(~DOSING)
	#plotobj1

	#png.file.name <- paste("ETA1_grid")
	#to.png.sqr(plotobj1,png.file.name)

	#plotobj2 <- NULL
	#plotobj2 <- ggplot(fitdata)
	#titletext <- expression(atop("ETA2 v's ROUTE",
	#												atop("BY ROUTE",
	#														 "coloured by DOSING")))
	#plotobj2 <- plotobj2 + ggtitle(titletext)
	#plotobj2 <- plotobj2 + geom_point(aes(x=ROUTE, y=ETA2, colour = DOSING), shape=1)
	#plotobj2 <- plotobj2+ scale_x_discrete(name="ROUTE")
	#plotobj2 <- plotobj2+ scale_y_continuous(name="ETA2")
	#plotobj2 <- plotobj2 + facet_wrap(~DOSING)
	#plotobj2

	#png.file.name <- paste("ETA2_grid")
	#to.png.sqr(plotobj2,png.file.name)

#--------------------------------------------------------------------------------------------------
# Plot individual fits and population fits by individual
  plotobj <- NULL
  plotobj <- ggplot(fitdata[fitdata$DV>0.001,])
  titletext <- expression(atop("Individual fits",
													atop("Black = observed",
															 "Red = individual prediction, blue = population prediction")))
	plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + geom_point(aes(x=TAD, y=DV), shape=1)
  #plotobj <- plotobj + geom_point(aes(x=TIME, y=DV), shape=1)
  plotobj <- plotobj + geom_line(aes(x=TAD, y=PRED), colour = "blue", data=fitdata[fitdata$DV>0.001&fitdata$STUDY==8056,])	#to account for non day 1 dataset
	plotobj <- plotobj + geom_line(aes(x=TAD, y=IPRED), colour = "red", data=fitdata[fitdata$DV>0.001&fitdata$STUDY==8056,])	#to account for non day 1 dataset
	plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour = "blue")		#not TAD due to graphical errors due to two IPREDs
	plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour = "red")		#not TAD due to graphical errors due to two PREDs
	plotobj <- plotobj + geom_hline(yintercept=0.25, linetype = 2, colour = "darkgreen")
 #plotobj <- plotobj + geom_hline(yintercept=0.3, linetype = 2, colour = "brown")
 #plotobj <- plotobj + geom_hline(yintercept=0.2, linetype = 2, colour = "purple")
	plotobj <- plotobj + facet_wrap(~ ID)
  plotobj <- plotobj+ scale_x_continuous(name="Time after dose (h)", lim = c(0,6))
  plotobj <- plotobj+ scale_y_continuous(name="Lenalidomide (ug/mL)")
  #plotobj

	ggsave("Individual_concs_by_dose_ID_DOSELVL.png", width=90, height=50, units=c("cm"))

	plotobj <- plotobj+ scale_y_log10(name="Lenalidomide Conc (ug/mL)", lim=c(0.01, 5))

	ggsave("Individual_concs_by_dose_ID_DOSELVL_log.png", width=90, height=50, units=c("cm"))

#--------------------------------------------------------------------------------------------------
# Plot 20 best/worst individual fits based on the median prediction errors.
# Not sure I like the use of median PE.  Should really be the sum of the PE's

#First calculate the Prediction Error
  fitdata$PE <- (fitdata$DV-fitdata$IPRED)/fitdata$IPRED
  #head(fitdata)

#Calculate and sort by Median Absolute PE
  #First calculate and order the fitdata by MDAPE
  #Subset for subjects with more than 1 observation to stop graphs crashing
  median.abs <- function(x) median(abs(x), na.rm=T)
  medianPE.df <- summaryBy(PE ~ ID, data=fitdata, FUN=c(median.abs,length))

  medianPE.df.sorted <- orderBy(~ PE.median.abs, data=medianPE.df)
  medianPE.df.sorted <- subset(medianPE.df.sorted, PE.length > 1)

	write.csv(medianPE.df.sorted,file="median.df.sorted.csv", row.names=F)

#Calculate and sort by SSE (sum of PE's squared)
	fitdata$SE <-  fitdata$PE*fitdata$PE
	#fitdata
	sum.abs <- function(x) sum(abs(x), na.rm=T)
  sumSE.df <- summaryBy(SE ~ ID, data=fitdata, FUN=c(sum.abs,length))

  sumSE.df.sorted <- orderBy(~ SE.sum.abs, data=sumSE.df)
  sumSE.df.sorted <- subset(sumSE.df.sorted, SE.length > 1)

	write.csv(sumSE.df.sorted,file="sum.df.sorted.csv", row.names=F)

#Calculate and sort by MSE
	MSE.df.sorted <- sumSE.df.sorted
	MSE.df.sorted <- rename(MSE.df.sorted, c("SE.sum.abs"="PEsum","SE.length"="nData")) #just makes it clearer to me
	MSE.df.sorted$nPar <- 3 # hard coded for the number of fixed effects
	MSE.df.sorted$DegFree <- MSE.df.sorted$nData - MSE.df.sorted$nPar
	MSE.df.sorted$MSE <- MSE.df.sorted$PEsum/MSE.df.sorted$DegFree

  MSE.df.sorted <- orderBy(~ MSE, data=MSE.df.sorted)
  MSE.df.sorted <- subset(MSE.df.sorted, nData > 1)

	write.csv(MSE.df.sorted,file="MSE.df.sorted.csv", row.names=F)

#Plot top 2 fits medianPE---------------------------
  IDtop <- head(subset(medianPE.df.sorted,PE.length>3)$ID, 2)
  fitdatatop <- fitdata[fitdata$ID %in% IDtop,]

	plotobjtop <- NULL
  plotobjtop <- ggplot(fitdatatop)
	titletext <- expression(atop("Best individual fits by medianPE",
													atop("Black = observed",
															 "Red = individual prediction, blue = population prediction")))
	plotobjtop <- plotobjtop + ggtitle(titletext)
	plotobjtop <- plotobjtop + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjtop <- plotobjtop + geom_line(aes(x=TIME, y=PRED), colour = "blue")
	plotobjtop <- plotobjtop + geom_line(aes(x=TIME, y=IPRED), colour = "red")
	plotobjtop <- plotobjtop + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
	plotobjtop <- plotobjtop + scale_x_continuous(name="Time after dose (h)",lim=c(0,24))
  plotobjtop <- plotobjtop + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjtop <- plotobjtop + facet_grid(~ ID)
	#plotobjtop

  to.png(plotobjtop,"CONC_vs_TAD_by_ID_best_medianPE")

	plotobjtop <- plotobjtop + scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjtop,"LOG_CONC_vs_TAD_by_ID_best_medianPE")

 #Plot bottom 2 fits medianPE------------------------
  IDbottom <- tail(medianPE.df.sorted$ID, 2)
  fitdatabottom <- fitdata[fitdata$ID %in% IDbottom,]

	plotobjbottom <- NULL
  plotobjbottom <- ggplot(fitdatabottom)
	titletext <- expression(atop("Worst individual fits by medianPE",
													atop("Black = observed",
															 "Red  = individual prediction, blue = population prediction")))
	plotobjbottom <- plotobjbottom + ggtitle(titletext)

	plotobjbottom <- plotobjbottom + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjbottom <- plotobjbottom + geom_line(aes(x=TIME, y=PRED), colour = "blue")
	plotobjbottom <- plotobjbottom + geom_line(aes(x=TIME, y=IPRED), colour = "red")
  plotobjbottom <- plotobjbottom + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
	plotobjbottom <- plotobjbottom + scale_x_continuous(name="Time after dose (h)",lim=c(0,24))
  plotobjbottom <- plotobjbottom + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjbottom <- plotobjbottom + facet_grid( ~ ID)
	#plotobjbottom

  to.png(plotobjbottom,"CONC_vs_TAD_by_ID_worst_medianPE")

	plotobjbottom <- plotobjbottom+ scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjbottom,"LOG_CONC_vs_TAD_by_ID_worst_medianPE")

#Plot top 2 fits MSE---------------------------
  IDtop <- head(subset(MSE.df.sorted,nData>3)$ID, 2)
  fitdatatop <- fitdata[fitdata$ID %in% IDtop,]

	plotobjtop <- NULL
  plotobjtop <- ggplot(fitdatatop)
	titletext <- expression(atop("Best individual fits by MSE",
													atop("Black = observed",
															 "Red = individual prediction, blue = population prediction")))
	plotobjtop <- plotobjtop + ggtitle(titletext)
	plotobjtop <- plotobjtop + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjtop <- plotobjtop + geom_line(aes(x=TAD, y=PRED), colour = "blue")
	plotobjtop <- plotobjtop + geom_line(aes(x=TAD, y=IPRED), colour = "red")
	plotobjtop <- plotobjtop + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
	plotobjtop <- plotobjtop + scale_x_continuous(name="Time after dose (h)")
  plotobjtop <- plotobjtop + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjtop <- plotobjtop + facet_wrap(~ID)
	#plotobjtop

  to.png(plotobjtop,"CONC_vs_TAD_by_ID_best_MSE")

	plotobjtop <- plotobjtop + scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjtop,"LOG_CONC_vs_TAD_by_ID_best_MSE")


#Plot bottom 2 fits MSE------------------------
  IDbottom <- tail(MSE.df.sorted$ID, 2)
  fitdatabottom <- fitdata[fitdata$ID %in% IDbottom,]

	plotobjbottom <- NULL
  plotobjbottom <- ggplot(fitdatabottom)
	titletext <- expression(atop("Worst individual fits by MSE",
													atop("Black = observed",
															 "Red  = individual prediction, blue = population prediction")))
	plotobjbottom <- plotobjbottom + ggtitle(titletext)

	plotobjbottom <- plotobjbottom + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjbottom <- plotobjbottom + geom_line(aes(x=TAD, y=PRED), colour = "blue")
	plotobjbottom <- plotobjbottom + geom_line(aes(x=TAD, y=IPRED), colour = "red")
  plotobjbottom <- plotobjbottom + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
  plotobjbottom <- plotobjbottom + scale_x_continuous(name="Time after dose (h)")
  plotobjbottom <- plotobjbottom + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjbottom <- plotobjbottom + facet_wrap(~ID)
	#plotobjbottom

  to.png(plotobjbottom,"CONC_vs_TAD_by_ID_worst_MSE")

	plotobjbottom <- plotobjbottom+ scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjbottom,"LOG_CONC_vs_TAD_by_ID_worst_MSE")
