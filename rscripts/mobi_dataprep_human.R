# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  library(reshape2)

# Customize ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy PKSim Simulation Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation data is output from MoBi
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
  
# Read in PKSim data
  file.dir <- "raw_data/MoBi_paper/"
  model.name <- "PKSim_withArterial_Human.xlsx"
  pksim.raw <- as.data.frame(read_excel(paste0(file.dir, model.name)))
  
# Rename columns
  pksim <- pksim.raw
  names(pksim) <- c(
    "TIME", 
    paste0(
      rep(c("pl", "bc"), times = 2), 
      rep(c("Venous", "Arterial"), each = 2)
    ),
    paste0(
      rep(c("pl", "bc", "is", "ic", "ti"), times = 4), 
      rep(c("Brain", "Heart", "Kidney", "Lung"), each = 5)
    )
  )
  
# Convert time to minutes
  pksim$TIME <- pksim$TIME*60
  
# Melt pksim data
  pksim.id <- c("TIME", "plVenous", "bcVenous", "plArterial", "bcArterial")
  pksim.melt <- melt(pksim, 
    id.vars = pksim.id, variable.name = "var", value.name = "val"
  )
  
# Split variable into tissue and compartment type
  pksim.melt$var <- as.character(pksim.melt$var)
  pksim.melt$comp <- with(pksim.melt, paste0(substr(var, 1, 2), "Tissue"))
  pksim.melt$TISSUE <- with(pksim.melt, substr(var, 3, nchar(var)))
  
# Label blood vessel compartments as afferent and efferent for all but lung
  pksim.melt$plAfferent <- pksim.melt$plArterial
  pksim.melt$bcAfferent <- pksim.melt$bcArterial
  pksim.melt$plEfferent <- pksim.melt$plVenous
  pksim.melt$bcEfferent <- pksim.melt$bcVenous
    
# Redefine values for lung so they are swapped
  cond <- pksim.melt$TISSUE == "Lung"
  pksim.melt$plAfferent[cond] <- pksim.melt$plVenous[cond]
  pksim.melt$bcAfferent[cond] <- pksim.melt$bcVenous[cond]
  pksim.melt$plEfferent[cond] <- pksim.melt$plArterial[cond]
  pksim.melt$bcEfferent[cond] <- pksim.melt$bcArterial[cond]
  
# Cast data
  obsdata <- dcast(pksim.melt[, !names(pksim.melt) == "var"],
    TIME + plAfferent + bcAfferent + plEfferent + bcEfferent + TISSUE ~ comp,
    value.var = "val"
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tidy MoBi Simulation Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulation data is output from MoBi
# Time units are hours
# Concentration units are umol/L
# Molar mass of lenalidomide = 259.26 g/mol
    
# Read in simulation data and rename columns
  model.name <- "PKSimBrain_Human.xlsx"
  brain_raw <- read_excel(paste0(file.dir, model.name))
  brain <- as.data.frame(brain_raw)
  names(brain) <- c(
    "TIME", 
    paste0(
      c(rep(c("pl", "bc"), times = 3), "is", "ic", "ti"),
      c(rep(c("Efferent", "Afferent"), each = 2), rep("Tissue", times = 5))
    )
  )
  
  model.name <- "PKSimHeart_Human.xlsx"
  heart_raw <- read_excel(paste0(file.dir, model.name))
  heart <- as.data.frame(heart_raw)
  names(heart) <- c(
    "TIME", 
    paste0(
      c(rep(c("pl", "bc"), times = 3), "is", "ic", "ti"),
      c(rep(c("Efferent", "Afferent"), each = 2), rep("Tissue", times = 5))
    )
  )
  
  model.name <- "PKSimKidney_Human.xlsx"
  kidney_raw <- read_excel(paste0(file.dir, model.name))
  kidney <- as.data.frame(kidney_raw)
  names(kidney) <- c(
    "TIME", 
    paste0(
      c(rep(c("pl", "bc"), times = 3), "is", "ic", "ti"),
      c(rep(c("Efferent", "Afferent"), each = 2), rep("Tissue", times = 5))
    )
  )
  
  model.name <- "PKSimLung_Human.xlsx"
  lung_raw <- read_excel(paste0(file.dir, model.name))
  lung <- as.data.frame(lung_raw)
  names(lung) <- c(
    "TIME", 
    paste0(
      c(rep(c("pl", "bc"), times = 3), "is", "ic", "ti"),
      c(rep(c("Afferent", "Efferent"), each = 2), rep("Tissue", times = 5))
    )
  )
  
# Add tissue type column
  brain$TISSUE <- "Brain"
  heart$TISSUE <- "Heart"
  kidney$TISSUE <- "Kidney"
  lung$TISSUE <- "Lung"
  
# Remove kidney urine concentration column
  kidney <- kidney[, !names(kidney) == "uriTissue"]
  
# Join datasets together and change units of time to minutes
  simdata <- rbind(brain, heart, kidney, lung)
  simdata$TIME <- simdata$TIME*60
  