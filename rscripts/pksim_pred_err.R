# Simulation processing script for Lenalidomide PKSim IV Model
# -----------------------------------------------------------------------------
# Incorporates observed data which was previously published:
# Rozewski DM, Herman SEM, Towns WH, Mahoney E, Stefanovski MR, Shin JD, et al.
# Pharmacokinetics and Tissue Disposition of Lenalidomide in Mice. AAPS J.
# 2012;14(4):872-82.
# -----------------------------------------------------------------------------
# Ready workspace
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32", 
      "C:/Users/hugjh001/Documents/len_pbpk", "C:/Users/Jim Hughes/Documents/GitRepos/len_pbpk")

    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2] | getwd() == wd[5]) {
      git.dir <- "C:/Users/hugjh001/Documents"
      reponame <- "len_pbpk"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "len_pbpk"
    } else if (getwd() == wd[6]) {
      git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
      reponame <- "len_pbpk"
    }
    rm("wd")
  }

# Set up environment with observed data present
  source(paste(git.dir, reponame, "rscripts",
    "data_iv.R", sep = "/"))

# Load additional libraries
  library(plyr)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Simulation Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Simulation data is output from PKSim
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
# Read in simulation data
  file.dir <- "raw_data/PKSim_paper/"
  model.names <- c("M.BBB", "M.BBB.HydroAll", "M.BBB.HydroBrain", "M.LitPGP", 
    "M.NoBBB", "M.NoBBB.HydroAll", "M.NoBBB.HydroPlas", "M.BBB.HydroComb")

  filename <- list(
    BBB = list(
      Dose0_5 = paste0(file.dir, model.names[1], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[1], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[1], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[1], ".10.xlsx")
    ),
    BBB.HydroAll = list(
      Dose0_5 = paste0(file.dir, model.names[2], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[2], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[2], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[2], ".10.xlsx")
    ),
    BBB.HydroBrain = list(
      Dose0_5 = paste0(file.dir, model.names[3], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[3], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[3], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[3], ".10.xlsx")
    ),
    LitPGP = list(
      Dose0_5 = paste0(file.dir, model.names[4], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[4], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[4], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[4], ".10.xlsx")
    ),
    NoBBB = list(
      Dose0_5 = paste0(file.dir, model.names[5], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[5], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[5], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[5], ".10.xlsx")
    ),
    NoBBB.HydroAll = list(
      Dose0_5 = paste0(file.dir, model.names[6], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[6], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[6], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[6], ".10.xlsx")
    ),
    NoBBB.HydroPlas = list(
      Dose0_5 = paste0(file.dir, model.names[7], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[7], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[7], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[7], ".10.xlsx")
    ),
    BBB.HydroComb = list(
      Dose0_5 = paste0(file.dir, model.names[8], ".0.5.xlsx"),
      Dose1_5 = paste0(file.dir, model.names[8], ".1.5.xlsx"),
      Dose5 = paste0(file.dir, model.names[8
        ], ".5.xlsx"),
      Dose10 = paste0(file.dir, model.names[8], ".10.xlsx")
    )
  )
  
  allsims <- ldply(filename, function(x) {
  # Start up read excel loop
    read.in <- lapply(x, function(y) {
      df <- as.data.frame(read_excel(y))  # read in excel as data.frame
      mgkg <- as.numeric(substr(y, nchar(y) - 7, nchar(y) - 5))  # determine dose
      if (is.na(mgkg)) mgkg <- 5  # if dose is na then it is 5
      if (mgkg == 0.1) mgkg <- 10
      df$DOSEMGKG <- mgkg  # fill column in the data with this info
    # change names so rbind can work before binding the data.frames together in a list
      names(df) <- str_replace(names(df), substr(y, 22, nchar(y) - 5), "")
      df
    })
    
  # List of data.frames is easily bound together with do.call
    simdata_raw <- do.call("rbind", read.in)
    
  # Remove hydrolysis metabolite columns
    simdata <- simdata_raw[!str_detect(names(simdata_raw), "Hydrolysis.Metabolite")]
    
  # Rename columns
    names(simdata) <- c(
      "TIME", "Venous Blood Plasma", "bone", "brain_pl", "brain_bc", "brain_is", 
      "brain_ic", "Brain Tissue", "Heart Tissue", "urine", "Kidney Tissue", 
      "Liver Tissue", "Lung Tissue", "Muscle Tissue", "Spleen Tissue", "DOSEMGKG" 
    )
    names(dataiv) <- c(
      "UID", "ID", "TADNOM", "DOSEMGKG", "DOSEMG", "AMT", "WT", "TIME", 
      "Venous Blood Plasma", "Brain Tissue", "Liver Tissue", "Muscle Tissue", 
      "Heart Tissue", "Spleen Tissue", "Lung Tissue", "Kidney Tissue"
    )
    
  # Update time units to minutes and concentration units to ng/mL (ug/L)
    simdata$TIME <- simdata$TIME*60
    umol2ug <- function(x) {
      x*259.26
    }
    simdata[-c(1, 16)] <- colwise(umol2ug)(simdata[-c(1, 16)])
    
  # Melt the data for relevant concentration columns
    plotdata <- melt(simdata[c(1,2,8,9,11:16)], id.vars = c("TIME", "DOSEMGKG"),
      variable.name = "Tissue", value.name = "Concentration")
    plotdata$TIME <- round(plotdata$TIME, 1)
    obsdata <- melt(dataiv[c(3:4, 9:16)], id.vars = c("TADNOM", "DOSEMGKG"),
      variable.name = "Tissue", value.name = "Concentration")
    names(obsdata) <- c("TIME", "DOSEMGKG", "TISSUE", "OBS")
    
    predata <- plotdata[plotdata$TIME %in% unique(obsdata$TIME), ]
    names(predata) <- c("TIME", "DOSEMGKG", "TISSUE", "PRED")
    alldata <- merge(obsdata, predata)
    cleandata <- alldata[!is.na(alldata$OBS),]
    cleandata
  })
  
  allsims$PRED_ERR <- (allsims$PRED - allsims$OBS)/allsims$OBS
  
  dlply(allsims, .(.id, TISSUE), function(df) {
    c(sqrt(mean(df$PRED_ERR^2, na.rm = T)), quantile(df$PRED_ERR, prob = c(0.05, 0.5, 0.95)))
  })
  
  newsims <- ddply(allsims, .(.id, DOSEMGKG, TISSUE), function(df) {
    df$relCMAX <- max(df$PRED, na.rm = T)/max(df$OBS, na.rm = T) 
    df
  })
  
  dlply(newsims, .(.id, TISSUE), function(df) {
    cmax <- df[!duplicated(df$DOSEMGKG), "relCMAX"]
    c(
      meanCmax = mean(cmax),
      sdCmax = sd(cmax),
      seCmax = sd(cmax)/length(cmax),
      CI90loCmax = quantile(cmax, probs = 0.05),
      medCmax = median(cmax),
      CI90hiCmax = quantile(cmax, probs = 0.95)
    )
  })
  
  meanpool <- ddply(allsims, .(.id, TIME, DOSEMGKG, TISSUE), function(df) {
    df$OBS <- mean(df$OBS, na.rm = T)
    unique(df[-7])
  })
  
  auc_fn <- function(conc, time) {
    auc <- c(NULL)
    for (i in 2:length(conc)) {
      h <- time[i] - time[i-1]
      auc [i-1] <- (conc[i-1] + conc[i])*h/2
    }
    return(sum(auc))
  }
  
  meanpool <- ddply(meanpool, .(.id, DOSEMGKG, TISSUE), function(df) {
    df$relAUC <- auc_fn(df$PRED, df$TIME)/auc_fn(df$OBS, df$TIME)
    df
  })
  
  dlply(meanpool, .(.id, TISSUE), function(df) {
    auc <- df[!duplicated(df$DOSEMGKG), "relAUC"]
    c(
      meanAUC = mean(auc),
      sdAUC = sd(auc),
      seAUC = sd(auc)/length(auc),
      CI90loAUC = quantile(auc, probs = 0.05),
      medAUC = median(auc),
      CI90hiAUC = quantile(auc, probs = 0.95)
    )
  })
  
  t_test_fn <- function(mod1, mod2, tissue) {
    df1_cmax <- newsims[with(newsims, .id == mod1 & TISSUE == tissue), ]
    df1_cmax <- df1_cmax[with(df1_cmax, !duplicated(DOSEMGKG)), "relCMAX"]
    df2_cmax <- newsims[with(newsims, .id == mod2 & TISSUE == tissue), ]
    df2_cmax <- df2_cmax[with(df2_cmax, !duplicated(DOSEMGKG)), "relCMAX"]
    df1_auc <- meanpool[with(meanpool, .id == mod1 & TISSUE == tissue), ]
    df1_auc <- df1_auc[with(df1_auc, !duplicated(DOSEMGKG)), "relAUC"]
    df2_auc <- meanpool[with(meanpool, .id == mod2 & TISSUE == tissue), ]
    df2_auc <- df2_auc[with(df2_auc, !duplicated(DOSEMGKG)), "relAUC"]
    list(
      cmax = t.test(df1_cmax, df2_cmax, mu = 0),
      auc = t.test(df1_auc, df2_auc, mu = 0)
    )
  }
  t_test_fn("LitPGP", "BBB", "Heart Tissue")
  
  t_test_prederr_fn <- function(mod1, mod2, tissue) {
    df1_prer <- allsims[with(allsims, .id == mod1 & TISSUE == tissue), "PRED_ERR"]
    df2_prer <- allsims[with(allsims, .id == mod2 & TISSUE == tissue), "PRED_ERR"]
    list(
      medPREDERR = c(median(df1_prer, na.rm = T), Mod2 = median(df2_prer, na.rm = T)),
      sePREDERR = c(sd(df1_prer, na.rm = T)/length(na.omit(df1_prer)), Mod2 = sd(df2_prer, na.rm = T)/length(na.omit(df2_prer))),
      ttestPREDERR = t.test(df1_prer, df2_prer, mu = 0)
    )
  }
  t_test_prederr_fn("LitPGP", "BBB", "Kidney Tissue")
  