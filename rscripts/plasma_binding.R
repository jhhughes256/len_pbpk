# Calculation of protein binding from equilibrium dialysis results
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "len_pbpk"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
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

# -----------------------------------------------------------------------------
# Load in data
  # data_folder <- "E:/Hughes/Data/RAW_NonClinical/protein_binding/20171009"
  data_folder <- paste(git.dir, reponame, "raw_data", sep = "/")

  filename_in <- paste0(data_folder, "/mouse plasma protein binding results_9282017_std passing_Short.xls")
  mouse_plas1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/Plasma protein binding_plasma_human_rerun_finalresults_9282017_Short.xls")
  human_plas1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  # filename_in <- paste0(data_folder, "/Plasmaproteinbinding_PBS_results_10022017_Short_171006133211.xls")
  # both_pbs1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/Plasmaproteinbinding_PBS_mitchedits_results_10022017_Short.xls")
  both_pbs1 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/mouse plasma 4h_results_10182017_Short.xls")
  mouse_plas2 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/human plasma protein binding_4h_LenaAPCI_results_10182017_Short.xls")
  human_plas2 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

  filename_in <- paste0(data_folder, "/Plasmaproteinbinding_4h_PBS_results_10182017_Short.xls")
  both_pbs2 <- read_excel(filename_in, sheet = "Lenalidomide", skip = 4)

# Combine data
  # any(names(mouse_plas) != names(human_plas) & names(mouse_plas) != names(both_pbs))
  alldata <- rbind(mouse_plas1, human_plas1, both_pbs1, mouse_plas2, human_plas2, both_pbs2)

# Fix column names
  name_correction <- c("Specified.Amount", "Calculated.Amount")
  names(alldata) <- str_replace_all(names(alldata),"[ ()#]",".")
  names(alldata)[names(alldata) == "Amount"] <- name_correction

# Take subset of sample data and add covariate info
  subdata <- alldata[alldata$Sample.Type == "Unknown Sample" & !is.na(alldata$Sample.Type),]

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

  out <- ddply(subdata, .(species, conc, ID), function(x) {
    Cp <- x$dv[x$vehicle == "plasma"]  # total plasma conc
    Cd <- x$dv[x$vehicle == "buffer"]  # free dialysate conc
    Ci <- unique(x$conc)*1000  # initial plasma conc
    # (Cp-Cd)*1.1/((Cp-Cd)*1.1+Cd)
    # browser()
    fb <- (Cp-Cd)/Cp*100
    dr <- Cd/Ci*100
    pr <- Cp/Ci*100
    c(fraction_bound = fb, plas_recovered = pr, dial_recovered = dr, total_recovered = dr+pr)
  })

  fb_fun <- function(x, z, dv) {
    dc <- x - z
    dc*dv/(dc*dv+z)
  }

  filename_out <- paste(git.dir, "protein_binding.csv", sep = "/")
  write.csv(out, filename_out, row.names = F, quote = F)

  statsum <- ddply(out, .(species, conc), function(x) {
    fb_na <- is.na(x$fraction_unbound)
    n <- 3
    if (any(fb_na)) {
      n <- n - length(which(fb_na))
    }
    data.frame(
      var = c("fraction_bound", "plas_recovered", "dial_recovered", "total_recovered"),
      median = c(
        mean(x$fraction_bound, na.rm = T), mean(x$fraction_bound, na.rm = T),
        mean(x$fraction_bound, na.rm = T), mean(x$fraction_bound, na.rm = T)
      ),
      range = c(
        mean(x$fraction_bound, na.rm = T ), mean(x$fraction_bound, na.rm = T),
        mean(x$fraction_bound, na.rm = T), mean(x$fraction_bound, na.rm = T)
      ),
    )
  })
