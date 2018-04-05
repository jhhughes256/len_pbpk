###model_summary.r
##Goal: To observe models and their number of parameters, AIC etc.
##Note: Based heavily off of Model_Summary_aa_df_11July2016.r -> DF code

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()

# Set the working directory
   work.dir <- "E:/Hughes/Data/PK/FLAG"
   scriptname <- "model_summary"
   setwd(work.dir)

# Load libraries
  library(doBy)
  library(plyr)

# Source utility functions file
   source("E:/Hughes/functions_utility.r")

#Set NONMEM directory extension
  nm.dir <- "nm7"
  search.term <- paste("*",nm.dir, sep="")

#Get a list of subdirectory names
  dir.names <- dir(pattern=glob2rx(search.term))

###---------------------------------------------------------------
#The following need to be commented out or not depending if you have those types or runs or not
#For example, if you don't have a VPC run, then you must comment this code out.
#If you don't, then the code will crash.

#Find the ones with VPC in the name indicating VPC run
 indexoutVPC <- which(grepl(dir.names, pattern=glob2rx(c("*VPC*", "*CELGENE*", "*LOPEZ*"))))
 dir.names <- dir.names[-indexoutVPC]

#Find the ones with VPC in the name indicating VPC run
 #indexoutBS <- which(grepl(dir.names, pattern=glob2rx("*.bs.*")))
 #dir.names <- dir.names[-indexoutBS]

 #Find the ones with IMP in the name indicating IMPortance sampling run
 #indexoutIMP <- which(grepl(dir.names, pattern=glob2rx("*IMP*")))
 #dir.names <- dir.names[-indexoutIMP]


#Find the ones with MIX in the name indicating MIXture model run
 #indexoutMIX <- which(grepl(dir.names, pattern=glob2rx("*MIX*")))
 #dir.names <- dir.names[-indexoutMIX]
###-----------------------------------------------------------------

#Function to count number of parameters given a directory name
collate_SHK_AIC <- function(dir.name)
 {

  print(dir.name)

#Debug
  #work.dir <- master.dir
  #setwd(work.dir)
  #lst.file.path <-  paste(master.dir,"RUN001_1COMP_1ABS_PPV_CL.nm7/RUN001_1COMP_1ABS_PPV_CL.lst",sep="/")
  #ctl.file.name <- "RUN001_1COMP_1ABS_PPV_CL.ctl"
  #lst.file.name <- "RUN001_1COMP_1ABS_PPV_CL.lst"

  result1 <- NA
  #Makes use of the fact the *.lst file has a line starting "PARAMETER:" that then lists the initial values of all parameters.
  #This may be over more than 1 line, but the parameters end before the line starting GRADIENT:

  lst.file.name <- gsub(nm.dir,"lst",dir.name)
  ctl.file.name <- gsub(nm.dir,"ctl",dir.name)
  lst.file.path <- paste(work.dir,dir.name,lst.file.name, sep="/")

#Scrape data from the *.lst file
  if (file.exists(lst.file.path)==T)  #to screen for missing file
  {
   lst.lines <- readLines(lst.file.path)           #read all the lines of the lst file
    if (length(lst.lines) > 2)                         #to screen for running files where only 2 lines are present
     {

      #Calculate the number of parameters
      first.line <- grep("PARAMETER:",lst.lines)      #find the index number of the first parameter line
      first.line <- first.line[1]

      last.line <- grep("GRADIENT:",lst.lines)         #find the index number of the last parameter line
      last.line <- last.line[1] - 1

      param.lines <- lst.lines[first.line:last.line]   #subset all the parameter lines

      param.lines  <- scan(textConnection(param.lines), what="char", quiet=T)  #turn into vector
      closeAllConnections()

      npar <- length(param.lines)-1
      print(npar)

      #Get the termination status
      #1 is successful, 0 is terminated
      term.line <- grep("#TERM",lst.lines) + 1
      term.text <- lst.lines[term.line]
      term.code <- grep("SUCCESSFUL",term.text)
      if (length(term.code) > 0) term.code <- "Successful"
      if (length(term.code) == 0) term.code <- "Terminated"
      print(term.code)

      #Get the covariance status - this will be absent if covariance step failed
      #1 is successful, 0 is failed
      cov.line <- grep("STANDARD ERROR OF ESTIMATE",lst.lines) #cannot just rely upon "STANDARD ERROR" because Mixture models use that term in their output unrelated to $COV status
      if (length(cov.line) > 0) cov.code <- "OK"
      if (length(cov.line) == 0) cov.code <- "FAILED"
      print(cov.code)

      #Get the objective function value
           obj.line <- grep("OBJV:",lst.lines)      #find the index number of the first parameter line
           OBJ <- lst.lines[obj.line]
           OBJ <- scan(textConnection(OBJ), what="char", quiet=T)
           OBJ <- as.numeric(OBJ[2])
           print(OBJ)

      #Get the theta parameters
      first.line <- grep("FINAL PARAMETER ESTIMATE",lst.lines)      #find the index number of the first theta line
      first.line <- first.line[1]+6

      last.line <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS",lst.lines)         #find the index number of the last theta line
      last.line <- last.line[1] - 1      #first occurrence of this string!

      theta.lines <- lst.lines[first.line:last.line]   #subset all the theta lines
      theta.lines <- gsub("TH ","TH",theta.lines) #stop TH 1 being read as "TH","1"

      theta.lines  <- scan(textConnection(theta.lines), what="char", quiet=T)  #turn into vector
      closeAllConnections()

      theta.vals <- as.numeric(theta.lines)



      #Get the standard errors of the theta parameters
      if (cov.code=="OK")
      {
      first.line <- grep("STANDARD ERROR OF ESTIMATE",lst.lines)      #find the index number of the first se line
      first.line <- first.line[1]+6

      last.line <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS",lst.lines)         #find the index number of the last se line
      last.line <- last.line[2] - 1      #second occurrence of this string!

      theta.se.lines <- lst.lines[first.line:last.line]   #subset all the se lines
      theta.se.lines <- gsub("TH ","TH",theta.se.lines)

      theta.se.lines  <- scan(textConnection(theta.se.lines), what="char", quiet=T)  #turn into vector
      closeAllConnections()

      theta.se.vals <- as.numeric(theta.se.lines)

      theta.se.percent <- 100*(theta.se.vals/theta.vals)

      theta.se.max <- max(theta.se.percent, na.rm=T)

      } else
      {
        theta.se.max <- NA
      }
      print(theta.se.max)


      #Get the omega parameters
      first.line <- grep("OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS",lst.lines)      #find the index number of the first omega line
      first.line <- first.line[1]+3

      last.line <- grep("SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS",lst.lines)         #find the index number of the last omega line
      last.line <- last.line[1] - 1      #first occurrence of this string!

      omega.lines <- lst.lines[first.line:last.line]   #subset all the omega lines
      omega.lines <- gsub("ETA ","ETA",omega.lines) #stop ETA 1 being read as "ETA","1"

      omega.lines  <- scan(textConnection(omega.lines), what="char", quiet=T)  #turn into vector
      closeAllConnections()

      omega.vals <- as.numeric(omega.lines)


      #Get the standard errors of the omega parameters
      if (cov.code=="OK")
      {
      first.line <- grep("OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS",lst.lines)      #find the index number of the first se line
      first.line <- first.line[2]+3  # second occurrence of this string!

      last.line <- grep("SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS",lst.lines)         #find the index number of the last se line
      last.line <- last.line[2] - 1      #second occurrence of this string!

      omega.se.lines <- lst.lines[first.line:last.line]   #subset all the se lines
      omega.se.lines <- gsub("ETA ","ETA",omega.se.lines)

      omega.se.lines  <- scan(textConnection(omega.se.lines), what="char", quiet=T)  #turn into vector
      closeAllConnections()

      omega.se.vals <- as.numeric(omega.se.lines)

      omega.se.percent <- 100*(omega.se.vals/(omega.vals))

#This bit of code allows for off diagonal COV matrices.
#It only collects the SE associated with the diagonal elements otherwise it would also report
#on SE's for off diagonal elements (which can be very large).
#Comment out if you want the overall largest SE including off-diagonal elements


	  	  	if (length(omega.se.percent) == 4)				# 1 ETA
			{
			omega.se.percent <- omega.se.percent[c(4)]
			}

	  	  	if (length(omega.se.percent) == 9)				# 2 ETA's
			{
			omega.se.percent <- omega.se.percent[c(5,9)]
			}

	  		if (length(omega.se.percent) == 15)				# 3 ETA's
			{
			omega.se.percent <- omega.se.percent[c(6,10,15)]
			}

			if (length(omega.se.percent) == 22)				# 4 ETA's, there is a nice pattern to it ;-)
			{
			omega.se.percent <- omega.se.percent[c(7,11,16,22)]
			}

	  #calculates OMEGA %SE correctly.  Uses SD from the CORRelation matrix, and matches WFN output.
	  #see http://www.mail-archive.com/nmusers@globomaxnm.com/msg00845.html

      omega.se.max <- max(omega.se.percent, na.rm=T)

      } else
      {
        omega.se.max <- NA
      }
      print(omega.se.max)

      result1 <- data.frame(ctl.file.name,term.code,cov.code,OBJ,npar,theta.se.max,omega.se.max)  #turn into dataframe

      #Calculate the AIC
      result1$AIC <- result1$OBJ+2*result1$npar

      result1

     }
  }


  result2 <- NA
  shk.file.name <- gsub(nm.dir,"shk",dir.name)
  shk.file.path <- paste(work.dir,dir.name,shk.file.name, sep="/")

 #Debug
  #shk.file.path <- "D:/Foster/Canine_Pimobendan/PK_modelling/All_data/RUN7_PIMO_1COMP_TR1_FEDF.ctl.shk"
  #shk.file.name <- "RUN7_PIMO_1COMP_TR1_FEDF.ctl.shk"




#Scrape data from the *.shk file
  if (file.exists(shk.file.path)==T) {  #to screen for missing file
    shk.data <- read.table(shk.file.path, skip=1, header=T)           #read all the lines of the shk file
	  shk.data <- subset(shk.data, SUBPOP==1) #Only take SUBPOP 1 for mixture models

	  shrink.vals <- subset(shk.data, TYPE==4) #ETA shrinkage
	  lowestSHK <- min(shrink.vals[-c(1,2)])  #remove first 2 columns
    highestSHK <- max(shrink.vals[-c(1,2)])
    print(highestSHK)

	  pvalue.vals <- subset(shk.data, TYPE==3) #P values
	  lowestPVAL <- min(pvalue.vals[-c(1,2)])  #remove first 2 columns
    highestPVAL <- max(pvalue.vals[-c(1,2)])
    print(lowestPVAL)

    #result2 <- data.frame(ctl.file.name,lowestSHK,highestSHK,lowestPVAL,highestPVAL)  #turn into dataframe
    result2 <- data.frame(ctl.file.name,highestSHK,lowestPVAL)  #turn into dataframe
    result2
  }
  result <- cbind(result1,result2)
  }


#Run the function for every directory
 #Look at test list to find failed model runs - comment out extension (e.g. from nm72Win64 to nm72Win64off)
 #test <- mlply(dir.names, collate_SHK_AIC)

 #Run without errors to return as dataframe rather than list
  allrundata  <- mdply(dir.names, collate_SHK_AIC, .inform = TRUE)

#Order by AIC
  allrundata <- orderBy(~AIC, allrundata)

#Paste the name of the model in allAIC
  allrundata$ctl.file.name <- gsub(".ctl","",allrundata$ctl.file.name)

#Rename Columns
  allrundata <- rename(allrundata, c("ctl.file.name"="Run","term.code"="Minimization Step","cov.code"="Covariance Step"))

#Write to file
  write.csv(allrundata, file="TERM_SE_AIC_SHK_PVAL.csv", row.names=F)

	View(read.csv("TERM_SE_AIC_SHK_PVAL.csv"))
