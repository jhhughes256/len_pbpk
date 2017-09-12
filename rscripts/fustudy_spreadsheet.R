# Script for creation of a spreadsheet to catalogue samples
# -----------------------------------------------------------------------------
# Protein binding dialysis experiment conducted on 25/07/2017 - 03/08/2017
# Completed at Riffe Building Level 6, The Ohio State University, Columbus
# -----------------------------------------------------------------------------
# Load libraries
  library(plyr)  # arrange()

# Set working directory
  work.dir <- "C:/Users/Jim Hughes/Documents/GitRepos/len_pbpk/produced_data/"
  setwd(work.dir)

# Set up variables
  var.names <- c("Sample", "Incubation_(hours)", "Box", "Species",
    "Concentration_(uM)", "Tube", "Vehicle", "Note")
  var.hours <- c(2, 4, 6)
  var.species <- c("m", "h")  # m = mouse, h = human
  var.conc <- c(10, 3, 1, 0.3, 0.1, 0.03)  # units = uM
  var.rep <- 1:9
  var.solv <- c("T", "B")  # T = plasma, B = buffer

  var.df <- expand.grid(var.hours, var.species, var.conc, var.rep, var.solv)
  var.box <- ceiling(var.df$Var4/3)

# Create sample IDs
  var.id <- apply(var.df, 1, paste, collapse = "")
  var.id <- paste(substr(var.id, 1, 7), substr(var.id, 8, 9), sep = "-")
  which.id3 <- which(substr(var.id, 3, 3) == "1")  # 10uM
  which.id6 <- which(substr(var.id, 6, 7) == "00")  # 10, 3, 1 uM
  which.id7 <- which(substr(var.id, 7, 7) == "0")  # 10, 3, 1, 0.3, 0.1 uM

  var.id[which.id7] <- paste0(
    substr(var.id[which.id7], 1, 6),
    substr(var.id[which.id7], 8, 10)
  )
  var.id[which.id6] <- paste0(
    substr(var.id[which.id6], 1, 4),
    substr(var.id[which.id6], 7, 9)
  )
  var.id[which.id3] <- paste(
    substr(var.id[which.id3], 1, 2),
    substr(var.id[which.id3], 3, 7)
  )

# Combine columns
  sampdf <- arrange(data.frame(var.id, var.df[1], var.box, var.df[,-1]),
    Var1, var.box, Var2, -Var3, Var4)

# Add note for samples that don't exist
# This is Box 3 for 2 hours and 6 hours Incubation
# Ran out of mouse plasma, so there is no 0.1 or 0.03 tubes
  sampdf$Note <- ""
  which.missing <- which(
    sampdf$var.box == 3 &  # choose replicate box 3
    sampdf$Var1 != 4 &  # choose incubation times 2 & 6
    sampdf$Var2 == "m" &  # choose mouse samples
    sampdf$Var3 %in% c(0.1, 0.03)  # choose concentrations
  )
  sampdf$Note[which.missing] <- "no_sample"

# Add note for samples that came from an old membrane
  which.old <- which(
    sampdf$var.box == 1 &  # choose replicate box 1
    sampdf$Var1 == 6 &  # choose incubation times 6
    sampdf$Var2 == "h" &  # choose human samples
    sampdf$Var3 %in% c(0.3, 0.1, 0.03)  # choose concentrations
  )
  sampdf$Note[which.old] <- "old_membrane"

# Add note for uncertain quality of certain samples
  which.poor <- which(
    sampdf$var.box == 3 &  # choose replicate box 3
    sampdf$Var1 != 4 &  # choose incubation times 2 & 6
    sampdf$Var2 == "h" &  # choose human samples
    sampdf$Var3 %in% c(0.1, 0.03)  # choose concentrations
  )
  sampdf$Note[which.poor] <- "potentially_poor_quality"

# Add note for first priority samples
  which.first <- which(
    sampdf$Var1 == 4 &  # choose incubation times 4
    sampdf$Var3 %in% c(10, 0.3, 0.03)  # choose concentrations
  )
  sampdf$Note[which.first] <- "priority"

# Rename columns, change categorical levels, add columns for data entry and save
  names(sampdf) <- var.names
  levels(sampdf$Species) <- c("mouse", "human")
  levels(sampdf$Vehicle) <- c("plasma", "buffer")
  sampdf$Result <- ""
  filename.out <- paste0(work.dir, "Protein_Binding_Results.csv")
  #write.csv(sampdf, filename.out, row.names = F)

# This spreadsheet is then:
# - Saved as Excel workbook
# - Bottom border added to first row
# - View > Freeze Panes > Freeze Top Row
# - Alt+D, Alt+F, Alt+F (all cells)
# - Home > Format > AutoFit Column Width (all cells)
