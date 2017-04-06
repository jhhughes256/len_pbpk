###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_nonclin"
  setwd(master.dir)

# Load libraries
  library(ggplot2)
  library(doBy)
  library(Hmisc)
  library(plyr)
  library(grid)
  library(reshape)
  library(stringr)
  library(readxl)

# Source utility functions file
  source("E:/Hughes/functions_utility.r")

# Customize ggplot2 theme - R 2.15.3
  setthemebw2.1()

# Organise working and output directories
  working.dir <- paste(master.dir,"RAW_NonClinical",sep="/")
  workspacefilename <- paste(getwd(),"/",scriptname,".RData", sep="")

### ------------------------------------- Non-Clinical Data ------------------------------------- ###
### Updated from datacheck_front.r 			#reproducible
  file.name.in <- "RAW_NonClinical/All_Tissue_PK_Data_Summary.xls"										# read_excel used due to .xls format
  dataraw <- read_excel(file.name.in, col_types=c(rep("text",2),rep("numeric",19)) ,sheet="Data")		# col_types explicity stated due to "text" appearing in first two columns of data
  colnames(dataraw) <- gsub.all(c(" ","(",")","/"),rep(".",4),colnames(dataraw))						# remove troublesome characters from column names
  datanew <- dataraw[!(is.na(dataraw$Dose..mg.kg.)&is.na(dataraw$Sample.ID)),]							# remove junk rows
  datanew$Dose..mg.kg. <- as.numeric(datanew$Dose..mg.kg.)												# change class from character to numeric
  datanew[is.na(str_detect(datanew$Sample.ID,"_")),1] <- ""											# change NAs into blank strings (allows end.splitter to work)
  rowsplit <- c(which(str_detect(datanew$Sample.ID,"Data")),length(datanew$Sample.ID))					# find rows that split the data up (and also end of the data)

  dataiv <- datanew[c((rowsplit[1]+1):(rowsplit[2]-1)),]												# +1 & -1 used to avoid using the rows that split up the data as they contain no data
  averiv <- dataiv[!is.na(dataiv$Avg.time..min.),]														# rows that contain only average data

  datapo <- datanew[c((rowsplit[2]+1):(rowsplit[3]-1)),][1:7]											# square brackets remove junk columns
  averpo <- datapo[!is.na(datapo$Plasma.Avg..ng.mL.),]

  datapr <- datanew[c((rowsplit[3]+1):(rowsplit[4])),][1:7]											# no -1 here as its the end of the data
  averpr <- datapr[!is.na(datapr$Plasma.Avg..ng.mL.),]

### ---------------------------------------------IV---------------------------------------------- ###
  output.dir <- paste(working.dir,"/",scriptname,"_OutputIV",sep="")
 if(!file.exists(output.dir))
  {
   dir.create(output.dir)
  }

#-------------------------------------------------------------------------------------
 ### Column names
  # As presented
  names(dataiv)
  # Sorted
  sort(names(dataiv))
  # Structure
  str(dataiv)
#-------------------------------------------------------------------------------------
 ### Plot PK data
  with(dataiv, table(Dose..mg.kg., useNA = "always"))				#NA's from Repeat measures
  with(dataiv, table(Mouse.Wt..g., useNA = "always"))				#NA's from Repeat measures

  # Plot Conc vs TAFD
  CvTplot <- function(name,conc,dat)
  {
  plotobj <- NULL
  titletext <- paste("NonClin - Observed",name,"Concentrations\n")
  plotobj <- ggplot(data=dat)  #, colour=AMTMGf
  plotobj <- plotobj + geom_point(aes_string(x="Time..min.", y=conc), size=3, alpha=0.5, colour="blue")  #, colour=DOSEMGM2f
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj + scale_x_continuous("Time after first dose (minutes)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + facet_wrap(~Dose..mg.kg.)
  plotobj
  }
  CvTplot("Plasma","Plasma.DV..ng.mL.",dataiv)
  CvTplot("Brain","Brain..ng.mL.",dataiv)
  CvTplot("Liver","Liver..ng.mL.",dataiv)
  CvTplot("Muscle","Mscl..ng.mL.",dataiv)
  CvTplot("Heart","Hrt..ng.mL.",dataiv)
  CvTplot("Spleen","Spln..ng.mL.",dataiv)
  CvTplot("Lung","Lung..ng.mL.",dataiv)
  CvTplot("Kidney","Kidney..ng.mL.",dataiv)

  # Number of samples (includes repeats)
  nsub <- length(unique(dataiv$Sample.ID))
  nsub

#-------------------------------------------------------------------------------------
 ### Convert dataiv to standard format
  dataiv2 <- data.frame("SAMP" = dataiv$Sample.ID,"DOSEMGKG" = dataiv$Dose..mg.kg.)

  dataiv2$DOSEMG <- dataiv2$DOSEMGKG*dataiv$Mouse.Wt..g./1000

  dataiv2$AMT <- dataiv2$DOSEMG*1000	#dose in ug

  dataiv2$WT <- dataiv$Mouse.Wt..g.

  with(dataiv, table(Avg.time..min., useNA = "always"))
  dataiv2$TIME <- dataiv$Time..min.

  dataiv2$PLA <- dataiv$Plasma.DV..ng.mL.

  dataiv2$BRA <- dataiv$Brain..ng.mL.

  dataiv2$LVR <- dataiv$Liver..ng.mL.

  dataiv2$MSC <- dataiv$Mscl..ng.mL.

  dataiv2$HRT <- dataiv$Hrt..ng.mL.

  dataiv2$SPL <- dataiv$Spln..ng.mL.

  dataiv2$LUN <- dataiv$Lung..ng.mL.

  dataiv2$KID <- dataiv$Kidney..ng.mL.


#-------------------------------------------------------------------------------------
 ### Obtain standardised ID
end.splitter <- function(x) {
	#Debug
	#Still requires code for PO, PR datasets
	#x<-c("5hr_1","45min_2","1_5mg_1h_4","1_5mg_45m_3","1_5mg_45m_3_Repeat","1_5mg_20m_5_090925052635","300min_1","1_5hr_5","0_5mg_Gav_1_5h_1","Gavage_16h_1","")
  xsplit.t <- 0									#placeholder vector, first value to be removed at end of script
	tsplit.t <- 0
	uid.t <- 0
	uid.mark <- 0	 #to be subtracted from i to determine unique ID
  for(i in 1:length(x)) {  #START
    xstr <- str_sub(x[i], -1)  #final character in string
		ystr <- str_sub(x[i], -2, -2)  #penultimate character in string
		zstr <- str_sub(x[i], -3, -3)  #final character of time string (potentially)
    if(length(unique(0:9 %in% xstr)) == 2) {  #separates numbered values from text values
			if(match("_", ystr, nomatch=FALSE) == 1) {  #separates proper ID from junk ID - data here is normal ID
			  xsplit.t <- c(xsplit.t, xstr)  #save ID to vector
			  if(match("r", zstr, nomatch = FALSE) == 0) {  #remove hr from the pool
			    if(match("n", zstr, nomatch = FALSE) == 0) {  #remove min from the pool
				    suppressWarnings(avec <- as.numeric(str_sub(x[i], -5, -4)))
					  tstr <- str_sub(x[i], -4, -4)
				    if(match("h", zstr, nomatch = FALSE) == 0) {	#double digit m
						  if(!is.na(avec) == TRUE) {
								tsplit.t <- c(tsplit.t, as.character(avec/60))
							}else{  #single digit m
								tvec <- as.numeric(tstr)
								tsplit.t <- c(tsplit.t, as.character(tvec/60))
							}
						}else{  #non split and split h
							astr <- str_sub(x[i], -6, -6)
							if(length(unique(0:9 %in% astr)) == 2) {	#split time h
								tsplit.t <- c(tsplit.t, paste(astr, tstr, sep = "."))
							}else{  #single or double digit h
								if(!is.na(avec) == TRUE) {  #double digit h
									tsplit.t <- c(tsplit.t, avec)
								}else{  #single digit h
									tsplit.t <- c(tsplit.t, tstr)
								}
							}
						}
					}else{  #min from here on
						suppressWarnings(avec <- as.numeric(str_sub(x[i],-7,-6)))
						if(!is.na(avec) == TRUE){  #double or triple digit min
							suppressWarnings(tvec <- as.numeric(str_sub(x[i],-8,-6)))
							if(!is.na(tvec) == TRUE){  #triple digit min
								tsplit.t <- c(tsplit.t,as.character(tvec/60))
							}else{							#double digit min
								tsplit.t <- c(tsplit.t,as.character(avec/60))
							}
						}else{								#single digit min
							tvec <- as.numeric(str_sub(x[i],-6,-6))
							tsplit.t <- c(tsplit.t,as.character(tvec/60))
						}
					}
				}else{														#hr from here on
					tstr <- str_sub(x[i],-5,-5)
					astr <- str_sub(x[i],-7,-7)
					if(length(unique(0:9 %in% astr)) == 2){
						tsplit.t <- c(tsplit.t,paste(astr,tstr,sep="."))	#split time hr
					}else{
						tsplit.t <- c(tsplit.t,tstr)						#single digit hr
					}
				}
			}else{											#data here is junk ID attached to end of normal ID
				temp <- as.numeric(tail(xsplit.t,1))+1
				xsplit.t <- c(xsplit.t,temp)					#set new number to replace junk with last value + 1
				tvec <- as.numeric(str_sub(x[i],-18,-17))	#time must be of format "double digit min"
				tsplit.t <- c(tsplit.t,as.character(tvec/60))
        }
		}else{
		   if(str_detect(x[i],"_")){					 #data here is for the repeated concentrations
				xsplit.t <- c(xsplit.t,tail(xsplit.t,1))  #use last number as this is a repeat
				uid.mark <- uid.mark+1
				tsplit.t <- c(tsplit.t,tail(tsplit.t,1))
			}else{										 #data here is missing sample names
				temp <- as.numeric(tail(xsplit.t,1))+1
				xsplit.t <- c(xsplit.t,temp)			    #set new number to replace junk with last value + 1
				temp2 <- tail(tsplit.t,1)
				tsplit.t <- c(tsplit.t,temp2)			 #use last value as all should be from the same group
			}
      }
			uid.t <- c(uid.t,i-uid.mark)
   }
	xsplit <- xsplit.t[-1]	  #strip first value from vector to leave desired output
	tsplit <- tsplit.t[-1]
	uid <- uid.t[-1]
	output <- data.frame(uid,xsplit,as.numeric(tsplit)*60)
	colnames(output) <- c("UID","ID","TADNOM")
	output

#Debug Output
#  UID ID TADNOM
#1  1  1   300
#2  2  2    45
#3  3  4    60
#4  4  3    45
#5  4  3    45
#6  5  4    20
#7  6  1   300
#8  7  5    90
#9  8  1    90
#10  9  1   960
#11 10  2   960
}

  IDiv <- end.splitter(dataiv2$SAMP)
  dataiv3 <- cbind(IDiv,dataiv2[-1])
  ivdfall <- orderBy(~DOSEMGKG+TIME, data=dataiv3)

#-------------------------------------------------------------------------------------
 ### Collate average values in standard format

  averiv2 <- data.frame("DOSEMGKG" = averiv$Dose..mg.kg., "TIME" = averiv$Avg.time..min.)

  averiv2$TADNOM <- IDiv$TADNOM[!is.na(dataiv$Avg.time..min.)]

  averiv2$PLA <- averiv$Plasma.Avg..ng.mL.

  averiv2$BRA <- averiv$Brain.Avg..ng.mL.

  averiv2$LVR <- averiv$Liver.Avg..ng.mL.

  averiv2$MSC <- averiv$Mscl.Avg..ng.mL.

  averiv2$HRT <- averiv$Heart.Avg..ng.mL.

  averiv2$SPL <- averiv$Spln.Avg..ng.mL.

  averiv2$LUN <- averiv$.Lung.Avg..ng.mL.

  averiv2$KID <- averiv$.Kidney.Avg..ng.mL.

  ivdfave <- orderBy(~DOSEMGKG+TIME, data=averiv2)

#-------------------------------------------------------------------------------------
  # Check subject numbers
  with(ivdfall, table(ID))
  with(ivdfall, table(UID))

  # Check the PK sample times
  with(ivdfall, table(DOSEMGKG,TADNOM))
  #!!! 25min TADNOM outlier !!!#

  # Check the dose columns
  with(ivdfall, table(DOSEMG))
  with(ivdfall, table(DOSEMG,DOSEMGKG))	#dose by dose group

#-------------------------------------------------------------------------------------
 ### Calculate dose normalized concentrations
  # Calculate dose normalised DV
  # Units are ng/ml per mg

  ivdfall2 <- melt(ivdfall, id=c("UID","ID","TADNOM","DOSEMGKG","DOSEMG","AMT","WT","TIME"))
  ivdfave2 <- melt(ivdfave, id=c("DOSEMGKG","TIME","TADNOM"))
  colnames(ivdfall2) <- c(head(colnames(ivdfall2),8),"TISSUE","DV")
  colnames(ivdfave2) <- c(head(colnames(ivdfave2),3),"TISSUE","DV")
  # Check distribution of DV
  distplot <- function(input,output) {
    plotobj <- ggplot(input, aes(DV))
    plotobj <- plotobj + geom_histogram()
    plotobj <- plotobj + facet_wrap(~TISSUE,ncol=4)
    filename.out <- paste(output.dir,paste("Histogram_DV",output,sep="_"),sep="/")
    to.png(plotobj,filename.out)
    plotobj

    plotobj <- ggplot(input, aes(log(DV)))
    plotobj <- plotobj + geom_histogram()
    plotobj <- plotobj + facet_wrap(~TISSUE,ncol=4)
    filename.out <- paste(output.dir,paste("Histogram_DVlog",output,sep="_"),sep="/")
    to.png(plotobj,filename.out)
    plotobj
  }
  distplot(ivdfall2,"ALL")
  distplot(ivdfave2,"AVE")

  ivdfall2$DVNORM <- ivdfall2$DV/ivdfall2$DOSEMG

  ivdfall2$MDV <- ifelse(is.na(ivdfall2$DV),1,0)

#-------------------------------------------------------------------------------------
  # Count missing covariate data
  # Missing by Study
  covdata <- subset(ivdfall, select=c("UID","ID","WT","DOSEMGKG"))

  # Reassign missing
  covdata[covdata==-1] <- NA

  # Finish off
  missingbysample <- ddply(covdata, .(ID), colwise(calculate.percent.missing))

  filename.out <- paste(output.dir,"Missing_by_Sample.csv",sep="/")
  write.csv(missingbysample, file=filename.out, row.names=F)

  # Missing by Subject
  # Requires classify.missing function before it works
  missingbysubject <- ddply(covdata, .(ID,UID), colwise(calculate.percent.missing))
  missingsummary <- missingbysubject[missingbysubject$DOSEMGKG != 0,]
  missingsummary <- orderBy(~UID,missingsummary)

  filename.out <- paste(output.dir,"Missing_by_Subject.csv",sep="/")
  write.csv(missingsummary, file=filename.out, row.names=F)

#-------------------------------------------------------------------------------------
  # Subset covariates

  #Keeps missing as -1 - use for categorical summary
  dataallone.iv <- lapplyBy(~UID, data=ivdfall,  oneperID)
  dataallone.iv <- bind.list(dataallone.iv)
  dim(dataallone.iv)

  #Sets missing to NA - use for continuous summary
  covdataone <- lapplyBy(~UID, data=covdata,  oneperID)
  covdataone <- bind.list(covdataone)
  dim(covdataone)

  dataallone.iv$UIDf <- as.factor(dataallone.iv$UID)
  dataallone.iv$IDf <- as.factor(dataallone.iv$ID)
  dataallone.iv$DOSEMGKGf <- as.factor(dataallone.iv$DOSEMGKG)

  dataallone.iv2 <- melt(dataallone.iv, id=c("UID","ID","TADNOM","DOSEMGKG","DOSEMG","AMT","WT","TIME","UIDf","IDf","DOSEMGKGf"))
  colnames(dataallone.iv2) <- c(head(colnames(dataallone.iv),11),"TISSUE","DV")

#-------------------------------------------------------------------------------------
  # Summary of study characteristics
  # Do all subjects have PK data
  testfunc <- function(DV,dat) {
    DVtest <- summaryBy(as.formula(paste(DV,"~UID",sep="")), data=dat, FUN=mean, na.rm=T)
    DVtest <- DVtest[is.na(DVtest[2])==T,]
    DVtestID <- DVtest$UID
    DVtestID
  }
  testfunc("PLA",ivdfall)
  testfunc("BRA",ivdfall)
  testfunc("LVR",ivdfall)
  testfunc("MSC",ivdfall)
  testfunc("HRT",ivdfall)
  testfunc("SPL",ivdfall)
  testfunc("LUN",ivdfall)
  testfunc("KID",ivdfall)
  # Do all subjects have dose data
  testfunc("AMT",ivdfall)	#repeat samples

  # DV count by Dosegroup
  # Calculates data for Report Table 1
  countfunc <- function(DV,dat) {
    DVcount <- summaryBy(as.formula(paste(DV,"~DOSEMGKG",sep="")), data=dat, FUN=lengthNA)
    names(DVcount) <- c("DoseGroup","DVcount")
    DVcount
  }
  PLAcount <- countfunc("PLA",ivdfall)
  BRAcount <- countfunc("BRA",ivdfall)
  LVRcount <- countfunc("LVR",ivdfall)
  MSCcount <- countfunc("MSC",ivdfall)
  HRTcount <- countfunc("HRT",ivdfall)
  SPLcount <- countfunc("SPL",ivdfall)
  LUNcount <- countfunc("LUN",ivdfall)
  KIDcount <- countfunc("KID",ivdfall)

  #Subject count by Dosegroup
  SUBcount <- ddply(ivdfall, .(DOSEMGKG), function(df) count.unique(df$UID))
  names(SUBcount) <- c("DoseGroup","SUBcount")
  SUBcount

  #Dose count by Dosegroup
  AMTcount <- ddply(ivdfall, .(DOSEMGKG), function(df) lengthNA(df$AMT))
  names(AMTcount) <- c("DoseGroup","AMTcount")
  AMTcount

  # Average DV and AMT per Subject
  DVsum <- cbind(rbind(PLAcount,BRAcount,LVRcount,MSCcount,HRTcount,SPLcount,LUNcount,KIDcount),SUBcount[-1],AMTcount[-1])
  DVsum$DVperSUB <- round(DVsum$DVcount/DVsum$SUBcount,0)
  DVsum$AMTperSUB <- round(DVsum$AMTcount/DVsum$SUBcount,0)
  DVsum$TISSUE <- rep(c("PLA","BRA","LVR","MSC","HRT","SPL","LUN","KID"),each=5)
  DVsum

  filename.out <- paste(output.dir,"DVsum.csv",sep="/")
  write.csv(DVsum, file=filename.out)


  # Count missing DV data

  # Missing DV by TADNOM
  DVdata <- subset(ivdfall2, select=c(UID,ID,TADNOM,DV,MDV))

  missingDVbyTADNOM <- ddply(DVdata, .(TADNOM), colwise(calculate.percent.missing))
  missingDVbyTADNOM

  filename.out <- paste(output.dir,"missingDVbyTADNOM.csv",sep="/")
  write.csv(missingDVbyTADNOM, file=filename.out)


  # DV data present by TISSUE
  DV.present <- function(x)  if (any(is.numeric(x))==T) result <- 1 else result <- 0
  DVcountdata <-  ddply(ivdfall2, .(TISSUE,UID), function(df) DV.present(df$DV))


  withDVbyTISSUE <- ddply(DVcountdata, .(TISSUE), function(df) sum(df$V1))  #GOLD
  withDVbyTISSUE

  # Two subjects missing from all Tissue groups
  filename.out <- paste(output.dir,"DVwith_tissue.csv",sep="/")
  write.csv(withDVbyTISSUE, file=filename.out)

#-------------------------------------------------------------------------------------
  # Basic PK plot
  ivdfall3 <- ivdfall2
  BINnumber <- 3

  ivdfall3$DOSEMGKGf <- as.factor(ivdfall3$DOSEMGKG)
  ivdfall3$TISSUEf <- as.factor(ivdfall3$TISSUE)
  ivdfall3$WT_bin <- cut2(ivdfall3$WT, g=BINnumber)
  ivdfall3$DOSE_bin <- cut2(ivdfall3$DOSEMG, g=BINnumber)

  # Conc vs TAFD
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(ivdfall3,aes(colour=DOSEMGKGf))
  plotobj <- plotobj + geom_point(aes(TIME,DV),size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj + scale_x_continuous("Time (minutes)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("DoseGroup")
  plotobj <- plotobj + facet_wrap(~TISSUE,ncol=4)
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  plotobj

  filename.out <- paste(output.dir,"Overview_ConcObs_vs_TIME_by_TISSUE",sep="/")
  to.png(plotobj,filename.out)

#-------------------------------------------------------------------------------------
  # Basic PK plot

	dosefacetplot <- function(dfall,dose) {
	  dfsub <- dfall[dfall$DOSEMGKG==dose,]
  # Conc vs TAFD
    plotobj <- NULL
    titletext <- paste("Observed Concentrations\n")
    plotobj <- ggplot(dfsub,aes(colour=DOSEMGKGf))
    plotobj <- plotobj + geom_point(aes(TIME,DV),size=3, alpha=0.5)
    plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
    plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
    plotobj <- plotobj + scale_x_continuous("Time (minutes)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
    plotobj <- plotobj + scale_colour_discrete("DoseGroup")
    plotobj <- plotobj + facet_wrap(~TISSUE,ncol=4)
	  plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
    plotobj

    filename.out <- paste(output.dir,paste0("Overview_ConcObs_vs_TIME_by_TISSUE_",dose),sep="/")
    to.png(plotobj,filename.out)
  }

	 dosefacetplot(ivdfall3,0.5)
	 dosefacetplot(ivdfall3,1.5)
	 dosefacetplot(ivdfall3,5)
	 dosefacetplot(ivdfall3,10)
#-------------------------------------------------------------------------------------
  # Influence of Covariates

  # Function to plot by factor

  plotByFactor <- function(factorColname,factorText,dat) {
    spanfactor <- 1

  #Concentration plots
    plotobj <- NULL
    titletext <- paste("All Tissue Concentrations\n")
    plotobj <- ggplot(data=dat)
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y="DV", colour=factorColname), size=2, alpha=0.5)
    #plotobj <- plotobj + geom_smooth(aes_string(x="TIME", y="DV"), method=loess, span=spanfactor, se=F, size=1, colour="black")
    plotobj <- plotobj + ggtitle(titletext)
    plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
    plotobj <- plotobj + scale_x_continuous("Time after dose (mins)")
    plotobj <- plotobj + scale_colour_brewer(factorText, palette="Set1")
    plotobj <- plotobj + facet_wrap(as.formula(paste("~",factorColname,sep="")),ncol=2)
    plotobj

    filename.out <- paste(output.dir,"/",factorText,"_ConcObs_vs_TAD_facet",sep="")
    to.png.wx1(plotobj,filename.out)

  #Concentration plots
    plotobj <- NULL
    titletext <- paste("All Tissue Concentrations\n")
    plotobj <- ggplot(data=dat)
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y="DV", colour=factorColname), size=2, alpha=0.5)
    #plotobj <- plotobj + geom_smooth(aes_string(x="TIME", y="DV", colour=factorColname), method=loess, span=spanfactor, se=F, size=1)
    plotobj <- plotobj + ggtitle(titletext)
    plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
    plotobj <- plotobj + scale_x_continuous("Time after dose (mins)")
    plotobj <- plotobj + scale_colour_brewer(factorText, palette="Set1")
    plotobj

    filename.out <- paste(output.dir,"/",factorText,"_ConcObs_vs_TAD",sep="")
    to.png.wx1(plotobj,filename.out)
  }

  plotByFactor("DOSE_bin","Binned Dose (ug)",ivdfall3)
  plotByFactor("WT_bin","Binned Weight (g)",ivdfall3)



### ---------------------------------------------PO---------------------------------------------- ###
  output.dir <- paste(working.dir,"/",scriptname,"_OutputPO",sep="")
  if(!file.exists(output.dir)) {
    dir.create(output.dir)
  }

#-------------------------------------------------------------------------------------
 ### Column names
  # As presented
  names(datapo)
  # Sorted
  sort(names(datapo))
  # Structure
  str(datapo)
#-------------------------------------------------------------------------------------
 ### Plot PK data
  with(datapo, table(Dose..mg.kg., useNA = "always"))				#NA's from Repeat measures

  CvTplot("Plasma","Plasma.DV..ng.mL.",datapo)

  # Number of samples (includes repeats)
  nsub <- length(unique(datapo$Sample.ID))
  nsub

#-------------------------------------------------------------------------------------
 ### Convert datapo to standard format
  datapo2 <- data.frame("SAMP" = datapo$Sample.ID,"DOSEMGKG" = datapo$Dose..mg.kg.)

  with(datapo, table(Avg.time..min., useNA = "always"))
  datapo2$TIME <- datapo$Time..min.

  datapo2$DV <- datapo$Plasma.DV..ng.mL.

#-------------------------------------------------------------------------------------
 ### Obtain standardised ID
  IDpo <- end.splitter(datapo2$SAMP)
  datapo3 <- cbind(IDpo,datapo2[-1])
  podfall <- orderBy(~DOSEMGKG+TIME, data=datapo3)

#-------------------------------------------------------------------------------------
 ### Collate average values in standard format

  averpo2 <- data.frame("DOSEMGKG" = averpo$Dose..mg.kg., "TIME" = averpo$Avg.time..min.)

  averpo2$TADNOM <- IDpo$TADNOM[!is.na(datapo$Plasma.Avg..ng.mL.)]

  averpo2$DV <- averpo$Plasma.Avg..ng.mL.

  podfave <- orderBy(~DOSEMGKG+TIME, data=averpo2)

#-------------------------------------------------------------------------------------
  # Check subject numbers
  with(podfall, table(ID))
  with(podfall, table(UID))

  # Check the PK sample times
  with(podfall, table(DOSEMGKG,TADNOM))

  # Check the dose columns
  with(podfall, table(DOSEMGKG))	#dose by dose group

#-------------------------------------------------------------------------------------
 ### Calculate dose normalized concentrations
  # Calculate dose normalised DV
  # Units are ng/ml per mg

  podfall2 <- podfall
  podfave2 <- podfave
  # Check distribution of DV
  distplot <- function(input,output) {
    plotobj <- ggplot(input, aes(DV))
    plotobj <- plotobj + geom_histogram()
    plotobj <- plotobj + facet_wrap(~DOSEMGKG)
    filename.out <- paste(output.dir,paste("Histogram_DV",output,sep="_"),sep="/")
    to.png(plotobj,filename.out)
    plotobj

    plotobj <- ggplot(input, aes(log(DV)))
    plotobj <- plotobj + geom_histogram()
    plotobj <- plotobj + facet_wrap(~DOSEMGKG)
    filename.out <- paste(output.dir,paste("Histogram_DVlog",output,sep="_"),sep="/")
    to.png(plotobj,filename.out)
    plotobj
  }
  distplot(podfall2,"ALL")
  distplot(podfave2,"AVE")

  podfall2$DVNORM <- podfall2$DV/podfall2$DOSEMG

  podfall2$MDV <- ifelse(is.na(podfall2$DV),1,0)

#-------------------------------------------------------------------------------------
  # Count missing covariate data
  # Missing by Study
  covdata <- subset(podfall, select=c("UID","ID","DOSEMGKG"))

  # Reassign missing
  covdata[covdata==-1] <- NA

  # Finish off
  missingbysample <- ddply(covdata, .(ID), colwise(calculate.percent.missing))

  filename.out <- paste(output.dir,"Missing_by_Sample.csv",sep="/")
  write.csv(missingbysample, file=filename.out, row.names=F)

  # Missing by Subject
  # Requires classify.missing function before it works
  missingbysubject <- ddply(covdata, .(ID,UID), colwise(calculate.percent.missing))
  missingsummary <- missingbysubject[missingbysubject$DOSEMGKG != 0,]
  missingsummary <- orderBy(~UID,missingsummary)

  filename.out <- paste(output.dir,"Missing_by_Subject.csv",sep="/")
  write.csv(missingsummary, file=filename.out, row.names=F)

#-------------------------------------------------------------------------------------
  # Subset covariates

  #Keeps missing as -1 - use for categorical summary
  dataallone.po <- lapplyBy(~UID, data=podfall,  oneperID)
  dataallone.po <- bind.list(dataallone.po)
  dim(dataallone.po)

  dataallone.po$UIDf <- as.factor(dataallone.po$UID)
  dataallone.po$IDf <- as.factor(dataallone.po$ID)
  dataallone.po$DOSEMGKGf <- as.factor(dataallone.po$DOSEMGKG)

  dataallone.po2 <- dataallone.po

#-------------------------------------------------------------------------------------
  # Summary of study characteristics
  # Do all subjects have PK data
  testfunc("DV",podfall)

  # DV count by Dosegroup
  # Calculates data for Report Table 1
  DVcount <- countfunc("DV",podfall)

  #Subject count by Dosegroup
  SUBcount <- ddply(podfall, .(DOSEMGKG), function(df) count.unique(df$UID))
  names(SUBcount) <- c("DoseGroup","SUBcount")
  SUBcount


  # Average DV and AMT per Subject
  DVsum <- cbind(DVcount,SUBcount[-1])
  DVsum$DVperSUB <- round(DVsum$DVcount/DVsum$SUBcount,0)
  DVsum

  filename.out <- paste(output.dir,"DVsum.csv",sep="/")
  write.csv(DVsum, file=filename.out)


  # Count missing DV data

  # Missing DV by TADNOM
  DVdata <- subset(podfall2, select=c(UID,ID,TADNOM,DV,MDV))

  missingDVbyTADNOM <- ddply(DVdata, .(TADNOM), colwise(calculate.percent.missing))
  missingDVbyTADNOM

  filename.out <- paste(output.dir,"missingDVbyTADNOM.csv",sep="/")
  write.csv(missingDVbyTADNOM, file=filename.out)


  # DV data present by DOSEMGKG
  DV.present <- function(x)  if (any(is.numeric(x))==T) result <- 1 else result <- 0
  DVcountdata <-  ddply(podfall2, .(DOSEMGKG,UID), function(df) DV.present(df$DV))


  withDVbyDOSEMGKG <- ddply(DVcountdata, .(DOSEMGKG), function(df) sum(df$V1))  #GOLD
  withDVbyDOSEMGKG

  # Two subjects missing from all DOSEMGKG groups
  filename.out <- paste(output.dir,"DVwith_dosemgkg.csv",sep="/")
  write.csv(withDVbyTISSUE, file=filename.out)

#-------------------------------------------------------------------------------------
  # Basic PK plot
  podfall3 <- podfall2

  podfall3$DOSEMGKGf <- as.factor(podfall3$DOSEMGKG)

  # Conc vs TAFD
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(podfall3,aes(colour=DOSEMGKGf))
  plotobj <- plotobj + geom_point(aes(TIME,DV),size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj + scale_x_continuous("Time (minutes)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("DoseGroup")
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  plotobj

  filename.out <- paste(output.dir,"Overview_ConcObs_vs_TIME",sep="/")
  to.png(plotobj,filename.out)

### ---------------------------------------------PR---------------------------------------------- ###
  output.dir <- paste(working.dir,"/",scriptname,"_OutputPR",sep="")
  if(!file.exists(output.dir))
  {
    dir.create(output.dir)
  }

#-------------------------------------------------------------------------------------
 ### Column names
  # As presented
  names(datapr)
  # Sorted
  sort(names(datapr))
  # Structure
  str(datapr)
#-------------------------------------------------------------------------------------
 ### Plot PK data
  with(datapr, table(Dose..mg.kg., useNA = "always"))				#NA's from Repeat measures
  with(datapr, table(Mouse.Wt..g., useNA = "always"))				#NA's from Repeat measures

  # Plot Conc vs TAFD
  CvTplot("DV","Plasma.DV..ng.mL.",datapr)

  # Number of samples (includes repeats)
  nsub <- length(unique(datapr$Sample.ID))
  nsub

#-------------------------------------------------------------------------------------
 ### Convert datapr to standard format
  datapr2 <- data.frame("SAMP" = datapr$Sample.ID,"DOSEMGKG" = datapr$Dose..mg.kg.)

  datapr2$DOSEMG <- datapr2$DOSEMGKG*datapr$Mouse.Wt..g./1000

  datapr2$AMT <- datapr2$DOSEMG*1000	#dose in ug

  datapr2$WT <- datapr$Mouse.Wt..g.

  with(datapr, table(Avg.time..min., useNA = "always"))
  datapr2$TIME <- datapr$Time..min.

  datapr2$DV <- datapr$Plasma.DV..ng.mL.

#-------------------------------------------------------------------------------------
 ### Obtain standardised ID
  IDpr <- end.splitter(datapr2$SAMP)
  datapr3 <- cbind(IDpr,datapr2[-1])
  prdfall <- orderBy(~DOSEMGKG+TIME, data=datapr3)

#-------------------------------------------------------------------------------------
 ### Collate average values in standard format

  averpr2 <- data.frame("DOSEMGKG" = averpr$Dose..mg.kg., "TIME" = averpr$Avg.time..min.)

  averpr2$TADNOM <- IDpr$TADNOM[!is.na(datapr$Plasma.Avg..ng.mL.)]

  averpr2$DV <- averpr$Plasma.Avg..ng.mL.

  prdfave <- orderBy(~DOSEMGKG+TIME, data=averpr2)

#-------------------------------------------------------------------------------------
  # Check subject numbers
  with(prdfall, table(ID))
  with(prdfall, table(UID))

  # Check the PK sample times
  with(prdfall, table(DOSEMGKG,TADNOM))
  #!!! 25min TADNOM outlier !!!#

  # Check the dose columns
  with(prdfall, table(DOSEMG))
  with(prdfall, table(DOSEMG,DOSEMGKG))	#dose by dose group

#-------------------------------------------------------------------------------------
 ### Calculate dose normalized concentrations
  # Calculate dose normalised DV
  # Units are ng/ml per mg

  prdfall2 <- prdfall
  prdfave2 <- prdfave

  # Check distribution of DV
  distplot(prdfall2,"ALL")
  distplot(prdfave2,"AVE")

  prdfall2$DVNORM <- prdfall2$DV/prdfall2$DOSEMG

  prdfall2$MDV <- ifelse(is.na(prdfall2$DV),1,0)

#-------------------------------------------------------------------------------------
  # Count missing covariate data
  # Missing by Study
  covdata <- subset(prdfall, select=c("UID","ID","WT","DOSEMGKG"))

  # Reassign missing
  covdata[covdata==-1] <- NA

  # Finish off
  missingbysample <- ddply(covdata, .(ID), colwise(calculate.percent.missing))

  filename.out <- paste(output.dir,"Missing_by_Sample.csv",sep="/")
  write.csv(missingbysample, file=filename.out, row.names=F)

  # Missing by Subject
  # Requires classify.missing function before it works
  missingbysubject <- ddply(covdata, .(ID,UID), colwise(calculate.percent.missing))
  missingsummary <- missingbysubject[missingbysubject$DOSEMGKG != 0,]
  missingsummary <- orderBy(~UID,missingsummary)

  filename.out <- paste(output.dir,"Missing_by_Subject.csv",sep="/")
  write.csv(missingsummary, file=filename.out, row.names=F)

#-------------------------------------------------------------------------------------
  # Subset covariates

  #Keeps missing as -1 - use for categorical summary
  dataallone.pr <- lapplyBy(~UID, data=prdfall,  oneperID)
  dataallone.pr <- bind.list(dataallone.pr)
  dim(dataallone.pr)

  #Sets missing to NA - use for continuous summary
  covdataone <- lapplyBy(~UID, data=covdata,  oneperID)
  covdataone <- bind.list(covdataone)
  dim(covdataone)

  dataallone.pr$UIDf <- as.factor(dataallone.pr$UID)
  dataallone.pr$IDf <- as.factor(dataallone.pr$ID)
  dataallone.pr$DOSEMGKGf <- as.factor(dataallone.pr$DOSEMGKG)

  dataallone.pr2 <- dataallone.pr

#-------------------------------------------------------------------------------------
  # Summary of study characteristics
  # Do all subjects have PK data
  testfunc("DV",prdfall)
  # Do all subjects have dose data
  testfunc("AMT",prdfall)

  # DV count by Dosegroup
  # Calculates data for Report Table 1
  DVcount <- countfunc("DV",prdfall)

  #Subject count by Dosegroup
  SUBcount <- ddply(prdfall, .(DOSEMGKG), function(df) count.unique(df$UID))
  names(SUBcount) <- c("DoseGroup","SUBcount")
  SUBcount

  #Dose count by Dosegroup
  AMTcount <- ddply(prdfall, .(DOSEMGKG), function(df) lengthNA(df$AMT))
  names(AMTcount) <- c("DoseGroup","AMTcount")
  AMTcount

  # Average DV and AMT per Subject
  DVsum <- cbind(DVcount,SUBcount[-1],AMTcount[-1])
  DVsum$DVperSUB <- round(DVsum$DVcount/DVsum$SUBcount,0)
  DVsum$AMTperSUB <- round(DVsum$AMTcount/DVsum$SUBcount,0)
  DVsum

  filename.out <- paste(output.dir,"DVsum.csv",sep="/")
  write.csv(DVsum, file=filename.out)

  # Count missing DV data

  # Missing DV by TADNOM
  DVdata <- subset(prdfall2, select=c(UID,ID,TADNOM,DV,MDV))

  missingDVbyTADNOM <- ddply(DVdata, .(TADNOM), colwise(calculate.percent.missing))
  missingDVbyTADNOM

  filename.out <- paste(output.dir,"missingDVbyTADNOM.csv",sep="/")
  write.csv(missingDVbyTADNOM, file=filename.out)


  # DV data present by DOSEMGKG
  DV.present <- function(x)  if (any(is.numeric(x))==T) result <- 1 else result <- 0
  DVcountdata <-  ddply(prdfall2, .(DOSEMGKG,UID), function(df) DV.present(df$DV))


  withDVbyDOSEMGKG <- ddply(DVcountdata, .(DOSEMGKG), function(df) sum(df$V1))  #GOLD
  withDVbyDOSEMGKG

  # Two subjects missing from all Tissue groups
  filename.out <- paste(output.dir,"DVwith_dosemgkg.csv",sep="/")
  write.csv(withDVbyTISSUE, file=filename.out)

#-------------------------------------------------------------------------------------
  # Basic PK plot
  prdfall3 <- prdfall2
  BINnumber <- 3

  prdfall3$DOSEMGKGf <- as.factor(prdfall3$DOSEMGKG)
  prdfall3$WT_bin <- cut2(prdfall3$WT, g=BINnumber)
  prdfall3$DOSE_bin <- cut2(prdfall3$DOSEMG, g=BINnumber)

  # Conc vs TAFD
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(prdfall3,aes(colour=DOSEMGKGf))
  plotobj <- plotobj + geom_point(aes(TIME,DV),size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj + scale_x_continuous("Time (minutes)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("DoseGroup")
	 plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  plotobj

  filename.out <- paste(output.dir,"Overview_ConcObs_vs_TIME",sep="/")
  to.png(plotobj,filename.out)

#-------------------------------------------------------------------------------------
  # Influence of Covariates

  # Function to plot by factor
  plotByFactor("DOSE_bin","Binned Dose (ug)",prdfall3)
  plotByFactor("WT_bin","Binned Weight (g)",prdfall3)

#----------------------------------------------------------------------------
#Save workspace

save(list = ls(all = TRUE), file=workspacefilename)
