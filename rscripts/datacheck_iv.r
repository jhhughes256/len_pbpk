# Data Check Script for IV Data
# -----------------------------------------------------------------------------
# Data already published:
# Rozewski DM, Herman SEM, Towns WH, Mahoney E, Stefanovski MR, Shin JD, et al.
# Pharmacokinetics and Tissue Disposition of Lenalidomide in Mice. AAPS J.
# 2012;14(4):872-82.
# -----------------------------------------------------------------------------
# IV Data - All Tissues
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory
  git.dir <- "E:/Hughes/Git"
  reponame <- "len_pbpk"

# Setup directory
  source(paste(git.dir, reponame, "rscripts",
    "datacheck_load.R", sep = "/"))

  scriptname <- "datacheck_iv"
  output.dir <- paste(working.dir, "plot", scriptname, sep = "/")
  if(!file.exists(output.dir)) {
    dir.create(output.dir)
  }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explore Data for Cleaning and Data Extraction
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# View data names and structure
  names(rawiv)
  sort(names(rawiv))
  str(rawiv)

# Check distribution of doses and weights
  with(rawiv, table(Dose..mg.kg., useNA = "always"))
  with(rawiv, table(Mouse.Wt..g., useNA = "always"))

# Check NA's for non-binned data
  any(is.na(rawiv$Sample.ID))
  any(is.na(rawiv$Time..min.))

# Isolate NA's
  rawiv[which(is.na(rawiv$Dose..mg.kg.)), ]
  rawiv[which(is.na(rawiv$Mouse.Wt..g.)), ]
  rawiv[which(is.na(rawiv$Time..min.)), ]
  # NA's appear to be as a result of a repeat sample, it was ascertained from
  # Mitch that these are in fact for the same mouse. Should use the values from
  # each entry to inform the other to form one row.

# Number of samples (includes repeats)
  length(unique(rawiv$Sample.ID))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Clean & Extract Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Convert data for all data
  dataiv <- data.frame("ID" = rawiv$Sample.ID,"DOSEMGKG" = rawiv$Dose..mg.kg.)
  dataiv$DOSEMG <- dataiv$DOSEMGKG*rawiv$Mouse.Wt..g./1000
  dataiv$AMT <- dataiv$DOSEMG*1000	#dose in ug
  dataiv$WT <- rawiv$Mouse.Wt..g.
  dataiv$TIME <- rawiv$Time..min.
  dataiv$PLA <- rawiv$Plasma.DV..ng.mL.
  dataiv$BRA <- rawiv$Brain..ng.mL.
  dataiv$LVR <- rawiv$Liver..ng.mL.
  dataiv$MSC <- rawiv$Mscl..ng.mL.
  dataiv$HRT <- rawiv$Hrt..ng.mL.
  dataiv$SPL <- rawiv$Spln..ng.mL.
  dataiv$LUN <- rawiv$Lung..ng.mL.
  dataiv$KID <- rawiv$Kidney..ng.mL.

# Clean sample IDs
  IDiv <- end.splitter(dataiv$ID)
  dataiv <- dataiv %>%
    `[`(-1) %>%
    cbind(IDiv, .)

# Check the TADNOM against DOSEMGKG
  with(dataiv, table(DOSEMGKG, TADNOM))
  dataiv[which(dataiv$TADNOM == "25"), ]
  # There are one set of values stated to have TADNOM as 25, but values were
  # measured at 20mins, one at 13 which is odd
  # This TADNOM is likely meant to be 20
  dataiv$TADNOM[which(dataiv$TADNOM == "25")] <- "20"

# Fix repeated UID
# The original sample has NA WT
# The repeated sample has NA DOSEMGKG
  subiv <- dataiv[which(with(dataiv, TADNOM == "20" & DOSEMGKG == 5 & !is.na(dataiv$WT))), ]
  IDori <- dataiv[which(is.na(dataiv$WT)), ]
  IDrep <- dataiv[which(is.na(dataiv$DOSEMGKG)), ] %>% arrange(UID)
  tissue.mean <- function(x) {
    colwise(mean)(x[c("PLA", "BRA", "LVR", "MSC", "HRT", "SPL", "LUN", "KID")])
  }
  rbind(tissue.mean(subiv), tissue.mean(IDori), tissue.mean(IDrep))
  # it looks like the repeat row is no good, and was done due to the plasma
  # concs not working in the original test.

# Combine into one row, use DV and DOSEMGKG from orig, use WT & TIME from repeat
  IDori$TIME <- IDrep$TIME
  IDori$WT <- IDrep$WT
  IDori$DOSEMG <- IDori$DOSEMGKG*IDori$WT/1000
  IDori$AMT <- IDori$DOSEMG*1000  # units: ug
  IDori$PLA <- NA

# Then replace both samples with this combined sample and order the data.frame
  IDrem1 <- dataiv %>%
    `[`("WT") %>%
    is.na() %>%
    which()
  IDrem2 <- dataiv %>%
    `[`("DOSEMGKG") %>%
    is.na() %>%
    which()
  dataiv <- dataiv %>%
    `[`(-c(IDrem1, IDrem2), ) %>%
    rbind(IDori) %>%
    arrange(DOSEMGKG, TIME)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data Check
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Check subject numbers
  with(dataiv, table(ID))
  with(dataiv, table(UID))
  any(with(dataiv, table(UID)) > 2)
  # Successfully removed repeat sample

# Check PK dose data
  with(dataiv, table(DOSEMGKG, TADNOM))

# Check the dose columns
  with(dataiv, table(DOSEMG))
  hist(dataiv$DOSEMG)

# Check distribution of DV
  meltiv <- melt(dataiv, id = c("UID", "ID", "TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME"))
  colnames(meltiv) <- c(head(colnames(meltiv), 8), "TISSUE", "DV")
  iv.distplot(meltiv, "ALL.png")
  iv.distplot(ivdfave2, "AVE.png")

# Calculate dose normalized concentrations and mark missing DV
# Units are ng/ml per mg
   meltiv$DVNORM <- meltiv$DV/meltiv$DOSEMG
   meltiv$MDV <- ifelse(is.na(meltiv$DV), 1, 0)
   dataiv$DVNORM <- dataiv$DV/dataiv$DOSEMG
   dataiv$MDV <- ifelse(is.na(dataiv$DV), 1, 0)

# Plot PK data
  iv.CvTplot(meltiv[which(!is.na(meltiv$DOSEMGKG)), ])
  iv.CvTplot(meltiv[which(!is.na(meltiv$DOSEMGKG)), ], dosenorm = T)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Covariate Data Check
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Count missing covariate data
  missingbysample <- ddply(dataiv, .(DOSEMGKG), colwise(calculate.percent.missing))

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
  dataallone.iv <- lapplyBy(~UID, data=dataiv,  oneperID)
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
  testfunc("PLA",dataiv)
  testfunc("BRA",dataiv)
  testfunc("LVR",dataiv)
  testfunc("MSC",dataiv)
  testfunc("HRT",dataiv)
  testfunc("SPL",dataiv)
  testfunc("LUN",dataiv)
  testfunc("KID",dataiv)
  # Do all subjects have dose data
  testfunc("AMT",dataiv)	#repeat samples

  # DV count by Dosegroup
  # Calculates data for Report Table 1
  countfunc <- function(DV,dat) {
    DVcount <- summaryBy(as.formula(paste(DV,"~DOSEMGKG",sep="")), data=dat, FUN=lengthNA)
    names(DVcount) <- c("DoseGroup","DVcount")
    DVcount
  }
  PLAcount <- countfunc("PLA",dataiv)
  BRAcount <- countfunc("BRA",dataiv)
  LVRcount <- countfunc("LVR",dataiv)
  MSCcount <- countfunc("MSC",dataiv)
  HRTcount <- countfunc("HRT",dataiv)
  SPLcount <- countfunc("SPL",dataiv)
  LUNcount <- countfunc("LUN",dataiv)
  KIDcount <- countfunc("KID",dataiv)

  #Subject count by Dosegroup
  SUBcount <- ddply(dataiv, .(DOSEMGKG), function(df) count.unique(df$UID))
  names(SUBcount) <- c("DoseGroup","SUBcount")
  SUBcount

  #Dose count by Dosegroup
  AMTcount <- ddply(dataiv, .(DOSEMGKG), function(df) lengthNA(df$AMT))
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
  DVdata <- subset(dataiv2, select=c(UID,ID,TADNOM,DV,MDV))

  missingDVbyTADNOM <- ddply(DVdata, .(TADNOM), colwise(calculate.percent.missing))
  missingDVbyTADNOM

  filename.out <- paste(output.dir,"missingDVbyTADNOM.csv",sep="/")
  write.csv(missingDVbyTADNOM, file=filename.out)


  # DV data present by TISSUE
  DV.present <- function(x)  if (any(is.numeric(x))==T) result <- 1 else result <- 0
  DVcountdata <-  ddply(dataiv2, .(TISSUE,UID), function(df) DV.present(df$DV))


  withDVbyTISSUE <- ddply(DVcountdata, .(TISSUE), function(df) sum(df$V1))  #GOLD
  withDVbyTISSUE

  # Two subjects missing from all Tissue groups
  filename.out <- paste(output.dir,"DVwith_tissue.csv",sep="/")
  write.csv(withDVbyTISSUE, file=filename.out)

#-------------------------------------------------------------------------------------
  # Basic PK plot
  dataiv3 <- dataiv2
  BINnumber <- 3

  dataiv3$DOSEMGKGf <- as.factor(dataiv3$DOSEMGKG)
  dataiv3$TISSUEf <- as.factor(dataiv3$TISSUE)
  dataiv3$WT_bin <- cut2(dataiv3$WT, g=BINnumber)
  dataiv3$DOSE_bin <- cut2(dataiv3$DOSEMG, g=BINnumber)

  # Conc vs TAFD
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(dataiv3,aes(colour=DOSEMGKGf))
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

	 dosefacetplot(dataiv3,0.5)
	 dosefacetplot(dataiv3,1.5)
	 dosefacetplot(dataiv3,5)
	 dosefacetplot(dataiv3,10)
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

  plotByFactor("DOSE_bin","Binned Dose (ug)",dataiv3)
  plotByFactor("WT_bin","Binned Weight (g)",dataiv3)



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
