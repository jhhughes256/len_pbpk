# Calculation of protein binding from equilibrium dialysis results
# -----------------------------------------------------------------------------
# Set up directories
# Put the directory that all the data is located inside here
  data_folder <- "E:/Hughes/Data/RAW_NonClinical/protein_binding/20171009"

# Load libraries
  library(readxl)  # read in excel spreadsheets without conversion to csv
  library(stringr)  # for data extraction from sample names
  library(plyr)  # for easy application of functions on split datasets

# -----------------------------------------------------------------------------
# Load in all data
# 10/02 Data
  filename_in <- paste0(data_folder, "/mouse plasma protein binding results_9282017_std passing_Short.xls")
  mouse_plas1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/Plasma protein binding_plasma_human_rerun_finalresults_9282017_Short.xls")
  human_plas1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  # No Mitch edits
  # filename_in <- paste0(data_folder, "/Plasmaproteinbinding_PBS_results_10022017_Short_171006133211.xls")
  # both_pbs1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  # Mitch edits
  filename_in <- paste0(data_folder, "/Plasmaproteinbinding_PBS_mitchedits_results_10022017_Short.xls")
  both_pbs1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

# 10/18 Data
  filename_in <- paste0(data_folder, "/mouse plasma 4h_results_10182017_Short.xls")
  mouse_plas2 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/human plasma protein binding_4h_LenaAPCI_results_10182017_Short.xls")
  human_plas2 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/Plasmaproteinbinding_4h_PBS_results_10182017_Short.xls")
  both_pbs2 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)
  
# -----------------------------------------------------------------------------
# Combine data
  # any(names(mouse_plas) != names(human_plas) & names(mouse_plas) != names(both_pbs))
  alldata <- rbind(mouse_plas1, human_plas1, both_pbs1, mouse_plas2, human_plas2, both_pbs2)

# Fix column names
  name_correction <- c("Specified.Amount", "Calculated.Amount")
  names(alldata) <- str_replace_all(names(alldata),"[ ()#]",".")
  names(alldata)[names(alldata) == "Amount"] <- name_correction

# Take subset of sample data
  subdata <- alldata[alldata$Sample.Type == "Unknown Sample" & !is.na(alldata$Sample.Type),]

# Add covariate data using Filename string
  subdata$Species <- "mouse"
  subdata$Species[str_detect(subdata$Filename, "h")] <- "human"

  subdata$Vehicle <- "plasma"
  subdata$Vehicle[str_detect(subdata$Filename, "b")] <- "buffer"

  subdata$id <- 4
  subdata$id[str_detect(subdata$Filename, "5")] <- 5
  subdata$id[str_detect(subdata$Filename, "6")] <- 6
  subdata$id[str_detect(subdata$Filename, "7")] <- 7
  subdata$id[str_detect(subdata$Filename, "8")] <- 8
  subdata$id[str_detect(subdata$Filename, "9")] <- 9

  subdata$Conc <- 3
  subdata$Conc[str_detect(subdata$Filename, "0_3")] <- 0.3
  subdata$Conc[str_detect(subdata$Filename, "_1_")] <- 1
  subdata$Conc[str_detect(subdata$Filename, "0_03")] <- 0.03
  subdata$Conc[str_detect(subdata$Filename, "10_")] <- 10

# Clean up subdata with clean column names
  subdata <- with(subdata, data.frame(
    ID = id,
    name = Filename,
    dv = suppressWarnings(as.numeric(Calculated.Amount)),
    species = Species,
    vehicle = Vehicle,
    conc = Conc,
    peak_status = Peak.Status,
    area = Area,
    area_istd = ISTD.Area,
    area_ratio = Area.Ratio,
    rt = RT
  ))

# Handle negatives then split into human and mouse datasets
  subdata$dv[subdata$dv <= 0.3] <- NA

# Calculate fraction bound
# ddply splits the data according to the second argument
# e.g. it splits it into separate ID, conc and species terms ID being sample 
#   numbers, conc being concentrations
# The resulting split data consists of 2 rows, one for the plasma and one for 
#   the buffer 
  out <- ddply(subdata, .(species, conc, ID), function(x) {
    Cp <- x$dv[x$vehicle == "plasma"]  # total plasma conc
    Cd <- x$dv[x$vehicle == "buffer"]  # free dialysate conc
    Ci <- unique(x$conc)*1000  # initial plasma conc; convert from uM to nM
    fb <- (Cp-Cd)/Cp*100  # fraction bound
    dr <- Cd/Ci*100  # drug recovered in dialysate
    pr <- Cp/Ci*100  # drug recovered in plasma
    c(fraction_bound = fb, plas_recovered = pr, dial_recovered = dr, total_recovered = dr+pr)
  })

# Save to spreadsheet
  filename_out <- paste(git.dir, "protein_binding.csv", sep = "/")
  write.csv(out, filename_out, row.names = F, quote = F)
