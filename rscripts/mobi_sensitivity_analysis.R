# -----------------------------------------------------------------------------
# Ready workspace
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32", "C:/Users/hugjh001/Documents/len_pbpk")

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
    }
    rm("wd")
  }

# Load libraries
  library(readxl)
  library(stringr)
  library(plyr)
  library(reshape2)
  library(ggplot2)

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Sensitivity Analysis Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# SA data is output from PKSim
# Value column is dimensionless
  sapre <- "raw_data/MoBi_paper/"
  sasuf <- "Sumexp52.xlsx"
  mbra <- as.data.frame(read_excel(paste0(sapre, "SA_Brain", sasuf)))
  mhrt <- as.data.frame(read_excel(paste0(sapre, "SA_Heart", sasuf)))
  mkid <- as.data.frame(read_excel(paste0(sapre, "SA_Kidney", sasuf)))
  mlng <- as.data.frame(read_excel(paste0(sapre, "SA_Lung", sasuf)))
  hbra <- as.data.frame(read_excel(paste0(sapre, "SA_HumanBrain", sasuf)))
  hhrt <- as.data.frame(read_excel(paste0(sapre, "SA_HumanHeart", sasuf)))
  hkid <- as.data.frame(read_excel(paste0(sapre, "SA_HumanKidney", sasuf)))
  hlng <- as.data.frame(read_excel(paste0(sapre, "SA_HumanLung", sasuf)))
  
  mbra$Tissue <- "MouseBrain"
  mhrt$Tissue <- "MouseHeart"
  mkid$Tissue <- "MouseKidney"
  mlng$Tissue <- "MouseLung"
  hbra$Tissue <- "HumanBrain"
  hhrt$Tissue <- "HumanHeart"
  hkid$Tissue <- "HumanKidney"
  hlng$Tissue <- "HumanLung"
  
  sadata_art <- rbind(mbra, mhrt, mkid, hbra, hhrt, hkid)
  sadata_ven <- rbind(mlng, hlng)
  
# Subset data
# Only interested in arterial and venous blood
  sadata_art <- sadata_art[
    str_detect(sadata_art$Output, "ArterialBlood") & 
    str_detect(sadata_art$Output, "Plasma"),
  ]
  sadata_ven <- sadata_ven[
    str_detect(sadata_ven$Output, "VenousBlood") &
    str_detect(sadata_ven$Output, "Plasma"),
  ]
  
# Now can combine lung with other tissues
  sadata_raw <- rbind(sadata_art, sadata_ven)
  
# PK Parameter of interest: AUC_inf and C_max
  sadata_raw <- sadata_raw[sadata_raw$`PK-Parameter` %in% c("AUC_inf", "C_max"), ]
  
# Change id variables to factors and change values
  sadata_raw$PKParamf <- factor(sadata_raw$`PK-Parameter`)
  levels(sadata_raw$PKParamf) <- c("AUCinf", "Cmax")
  
# Remove sum of exponential parameters
  sadata_raw <- sadata_raw[!str_detect(sadata_raw$Parameter, "IV.0.5mgkg-Application_1-"),]
  
# Remove columns with no useful information
  sadata <- cbind(sadata_raw[8:9], sadata_raw[c(3, 7)])
  
# Remove all zeroes
  sadata$AbsValue <- abs(sadata$Value)
  sadata <- sadata[sadata$AbsValue > 0,]
  
# Now must determine most important parameters for each of the 8 Outputs
  ranksa <- ddply(sadata, .(Tissue, PKParamf), function(x) {
    sub <- x
    out <- sub[order(sub$AbsValue, decreasing = T),]
    cbind(out[1:2,], 1:2)
  })
 
  ranksa[ranksa$PKParamf == "Cmax",]
  