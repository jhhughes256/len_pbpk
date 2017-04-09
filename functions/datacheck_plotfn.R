# Plot functions used in datacheck_nonclin.R

iv.distplot <- function(input, output, output.dir) {
  plotobj <- ggplot(input, aes(log(DV)))
  plotobj <- plotobj + geom_histogram()
  plotobj <- plotobj + facet_wrap(~TISSUE, ncol = 4)
  filename.out <- paste(
    output.dir,
    paste("Histogram_DVlog", paste0(output, ".png"), sep = "_"),
    sep = "/"
  )  # paste
  plotobj
  ggsave(filename.out)
  plotobj
}

po.distplot <- function(input, output, output.dir) {
  plotobj <- ggplot(input, aes(log(DV)))
  plotobj <- plotobj + geom_histogram()
  plotobj <- plotobj + facet_wrap(~DOSEMGKG)
  filename.out <- paste(output.dir, paste0("Histogram_DVlog_", output, ".png"), sep = "/")
  plotobj
  ggsave(filename.out)
  plotobj
}

po.CvTplot <- function(dat, name, output.dir, dosenorm = F) {
  plotobj <- NULL
  plotobj <- ggplot(data = dat[which(!is.na(dat$DOSEMGKG)), ])
  if (!dosenorm) {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DV"),
      size = 2, alpha = 0.5, colour = "blue")
    plotobj <- plotobj + facet_wrap(~DOSEMGKG, ncol = 4)
    info <- c("Observed Concentrations\n", "DV", 7.2, 9.2)
  } else {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DVNORM"),
      size = 2, alpha = 0.5, colour = "blue")
    info <- c("Dose-Normalised Concentrations\n", "DVNORM", 7.2, 9.2)
  }
  plotobj <- plotobj + ggtitle(info[1])
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)\n")
  plotobj <- plotobj + scale_x_continuous("\nTime after first dose (minutes)")
  filename.out <- paste(
    output.dir,
    paste0("CvTplot_", name, "_", info[2], ".png"),
    sep = "/"
  )  # paste
  plotobj
  ggsave(filename.out, width = as.numeric(info[3]), height = as.numeric(info[4]))
  plotobj
}

po.CvTplot <- function(dat, name, output.dir, dosenorm = F) {
  plotobj <- NULL
  plotobj <- ggplot(data = dat[which(!is.na(dat$DOSEMGKG)), ])
  if (!dosenorm) {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DV"),
      size = 2, alpha = 0.5, colour = "blue")
    plotobj <- plotobj + facet_wrap(~TISSUE+DOSEMGKG, ncol = 4)
    info <- c("Observed Concentrations\n", "DV", 7.2, 18.4)
  } else {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DVNORM"),
      size = 2, alpha = 0.5, colour = "blue")
    plotobj <- plotobj + facet_wrap(~TISSUE, ncol = 3)
    info <- c("Dose-Normalised Concentrations\n", "DVNORM", 7.2, 9.2)
  }
  plotobj <- plotobj + ggtitle(info[1])
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)\n")
  plotobj <- plotobj + scale_x_continuous("\nTime after first dose (minutes)")
  filename.out <- paste(
    output.dir,
    paste0("CvTplot_", name, "_", info[2], ".png"),
    sep = "/"
  )  # paste
  plotobj
  ggsave(filename.out, width = as.numeric(info[3]), height = as.numeric(info[4]))
  plotobj
}


plotByFactor <- function(factorColname, factorText, dat, output.dir) {
  spanfactor <- 1

#Concentration plots
  plotobj <- NULL
  titletext <- paste("All Tissue Concentrations\n")
  plotobj <- ggplot(data = dat)
  plotobj <- plotobj + geom_point(aes_string(x = "TIME", y = "DV",
    colour = factorColname), size = 2, alpha = 0.5)
  #plotobj <- plotobj + geom_smooth(aes_string(x="TIME", y="DV"),
  #  method=loess, span=spanfactor, se=F, size=1, colour="black")
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj + scale_x_continuous("Time after dose (mins)")
  plotobj <- plotobj + scale_colour_brewer(factorText, palette = "Set1")
  plotobj <- plotobj + facet_wrap(~TISSUEf, ncol = 2)
  plotobj

  filename.out <- paste(output.dir, "/",
    factorText, "_ConcObs_vs_TAD_facet.png", sep = "")
  ggsave(filename.out)

#Concentration plots
  plotobj <- NULL
  titletext <- paste("All Tissue Concentrations\n")
  plotobj <- ggplot(data = dat)
  plotobj <- plotobj + geom_point(aes_string(x="TIME", y="DV",
    colour = factorColname), size = 2, alpha = 0.5)
  #plotobj <- plotobj + geom_smooth(aes_string(x="TIME", y="DV", colour=factorColname), method=loess, span=spanfactor, se=F, size=1)
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj + scale_x_continuous("Time after dose (mins)")
  plotobj <- plotobj + scale_colour_brewer(factorText, palette="Set1")
  plotobj

  filename.out <- paste(output.dir, "/",
    factorText, "_ConcObs_vs_TAD.png", sep = "")
  ggsave(filename.out)
}
