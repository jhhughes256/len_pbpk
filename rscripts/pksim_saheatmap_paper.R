# Sensitivity analysis processing script for Lenalidomide IV PKSim Model Paper
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
  library(ggplot2)
  library(cowplot)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy Sensitivity Analysis Data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# SA data is output from PKSim
# Value column is dimensionless
  filename0_5 <- "raw_data/PKSim_data/Sensitivity.Analysis.Tong.0.5.xlsx"
  # filename1_5 <- "raw_data/PKSim_data/Sensitivity.Analysis.Tong.1.5.xlsx"
  # filename5 <- "raw_data/PKSim_data/Sensitivity.Analysis.Tong.5.xlsx"
  # filename10 <- "raw_data/PKSim_data/Sensitivity.Analysis.Tong.10.xlsx"
  
  sadata0_5 <- as.data.frame(read_excel(filename0_5))
  # sadata1_5 <- as.data.frame(read_excel(filename1_5))
  # sadata5 <- as.data.frame(read_excel(filename5))
  # sadata10 <- as.data.frame(read_excel(filename10))
  sadata_raw <- sadata0_5
  
# Subset data
# PK Parameter of interest: AUC_inf and C_max
  sadata_raw <- sadata_raw[sadata_raw$`PK-Parameter` %in% c("AUC_inf", "C_max"), ]
  
# Not interested in hydrolysis metabolite or bone
  sadata_raw <- sadata_raw[!(str_detect(sadata_raw$Output, "Hydrolysis.Metabolite") | str_detect(sadata_raw$Output, "Bone")), ]
  
# Only interested in tissues and venous blood
  sadata_raw <- sadata_raw[str_detect(sadata_raw$Output, "Tissue") | str_detect(sadata_raw$Output, "VenousBlood"), ]
  
# Change id variables to factors and change values
  sadata_raw$Outputf <- factor(sadata_raw$Output)
  sadata_raw$PKParamf <- factor(sadata_raw$`PK-Parameter`)
  
  levels(sadata_raw$Outputf) <- c("Brain", "Heart", "Kidney", "Liver", "Lung", "Muscle", "Spleen", "Venous Blood")
  levels(sadata_raw$PKParamf) <- c("AUCinf", "Cmax")
  
# Remove columns with no useful information
  sadata <- cbind(sadata_raw[8:9], sadata_raw[c(3, 7)])
  
# Now must determine most important parameters for each of the 8 Outputs
  sadata$AbsValue <- abs(sadata$Value)
  ignoreSA <- c("IV.0.5mgkg-Dose", "Organism-Plasma protein scale factor",
    "Lenalidomide-ABCB1-Theoretical-Transporter concentration",
    "Lenalidomide-Hydrolysis.Tissues-Theoretical-Enzyme concentration",
    "Hydrolysis.Brain-Brain-Intracellular-Relative expression (normalized)")
  ranksa <- ddply(sadata, .(Outputf, PKParamf), function(x) {
    sub <- x[!x$Parameter %in% ignoreSA,]
    out <- sub[order(sub$AbsValue, decreasing = T),]
    cbind(out[1:10,], 1:10)
  })
  ranksa_list <- dlply(ranksa, .(PKParamf), function(x) {
    out <- x[order(x$`1:10`),]
    unique(out$Parameter)
  })
  
# Subset data for plotting and then clean up parameter names
  bestsa_auc <- sadata[sadata$Parameter %in% ranksa_list$AUCinf[1:10] & 
    sadata$PKParamf == "AUCinf", ]
  bestsa_cmax <- sadata[sadata$Parameter %in% ranksa_list$Cmax[1:10] &
    sadata$PKParamf == "Cmax", ]
  
  # factorise parameter names
  bestsa_auc$Paramf <- factor(bestsa_auc$Parameter)
  bestsa_cmax$Paramf <- factor(bestsa_cmax$Parameter)
  
  # change level names
  levels(bestsa_auc$Paramf) <- c("P-gp Liver Expr.", "P-gp Total Expr.", 
    "Hydro. Total Expr.", "Kidney Volume", "P-gp Km", "P-gp Vmax", 
    "Fraction Unbound", "Hydro. Metabolism", "Lipophilicity", "Muscle Volume")
  levels(bestsa_cmax$Paramf) <- c("P-gp Liver Expr.", "P-gp Total Expr.", 
    "GFR", "Kidney Volume", "P-gp Km", "P-gp Vmax", "Fraction Unbound", 
    "Lipophilicity", "Muscle Volume", "Hematocrit")
  
  # reorder factor levels
  bestsa_auc$Paramf <- factor(bestsa_auc$Paramf, 
    levels(bestsa_auc$Paramf)[c(5:6, 1:2, 8, 3:4, 10:9, 7)])
  bestsa_cmax$Paramf <- factor(bestsa_cmax$Paramf, 
    levels(bestsa_cmax$Paramf)[c(5:6, 1:2, 3, 10, 4, 9:8, 7)])
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot data
# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# Create heatmap for each plot, then cowplot them together
  max_auc <- max(abs(bestsa_auc$Value))
  max_cmax <- max(abs(bestsa_cmax$Value))
  
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_tile(aes(x = Outputf, y = Paramf, fill = Value), 
    data = bestsa_auc, colour = "white")
  p1 <- p1 + scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
    space = "Lab",  breaks = round(max_auc*c(-1, -0.5, 0, 0.5, 1), 1),
    name = "Sensitivity Index\n        AUCinf", midpoint = 0, limit = round(max_auc*c(-1.05, 1.05), 1))
  # p1 <- p1 + xlab("\nTissue AUCinf")
  p1 <- p1 + xlab("")
  p1 <- p1 + ylab("Model Parameter")
  p1 <- p1 + theme_minimal()
  p1
  
  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_tile(aes(x = Outputf, y = Paramf, fill = Value), 
    data = bestsa_cmax, colour = "white")
  p2 <- p2 + scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
    space = "Lab", breaks = round(max_cmax*c(-1, -0.5, 0, 0.5, 1), 1),
    name = "Sensitivity Index\n        Cmax", midpoint = 0, limit = round(max_cmax*c(-1.05, 1.05), 1))
  p2 <- p2 + xlab("\nTissue")
  p2 <- p2 + ylab("Model Parameter")
  p2 <- p2 + theme_minimal()
  p2
  
  plot_grid(p1, p2, ncol = 1)
  
  setwd("produced_data")
  ggsave("Figure 3.png", width = 23.2, height = 17.2, units = "cm")