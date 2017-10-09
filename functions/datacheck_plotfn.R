# Plot functions used in datacheck_nonclin.R

iv.distplot <- function(input, output) {
  plotobj <- ggplot(input, aes(log(DV)))
  plotobj <- plotobj + geom_histogram()
  plotobj <- plotobj + facet_wrap(~TISSUE, ncol = 4)
  plotobj
}

po.distplot <- function(input, output) {
  plotobj <- ggplot(input, aes(log(DV)))
  plotobj <- plotobj + geom_histogram()
  plotobj <- plotobj + facet_wrap(~DOSEMGKG)
  plotobj
}

iv.CvTplot <- function(dat, dosenorm = F) {
  plotobj <- NULL
  plotobj <- ggplot(data = dat[which(!is.na(dat$DOSEMGKG)), ])
  if (!dosenorm) {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DV"),
      size = 2, alpha = 0.5, colour = "blue")
    plotobj <- plotobj + facet_wrap(DOSEMGKG~TISSUE, ncol = 4)
    info <- c("Observed Concentrations\n", "DV", 7.2, 9.2)
  } else {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DVNORM"),
      size = 2, alpha = 0.5, colour = "blue")
    plotobj <- plotobj + facet_wrap(~TISSUE, ncol = 4)
    info <- c("Dose-Normalised Concentrations\n", "DVNORM", 7.2, 9.2)
  }
  plotobj <- plotobj + ggtitle(info[1])
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)\n")
  plotobj <- plotobj + scale_x_continuous("\nTime after first dose (minutes)")
  plotobj
}

po.CvTplot <- function(dat, dosenorm = F) {
  plotobj <- NULL
  plotobj <- ggplot(data = dat[which(!is.na(dat$DOSEMGKG)), ])
  if (!dosenorm) {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DV"),
      size = 2, alpha = 0.5, colour = "blue")
    plotobj <- plotobj + facet_wrap(~DOSEMGKG, ncol = 4)
    info <- c("Observed Concentrations\n", "DV", 7.2, 18.4)
  } else {
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y = "DVNORM"),
      size = 2, alpha = 0.5, colour = "blue")
    info <- c("Dose-Normalised Concentrations\n", "DVNORM", 7.2, 9.2)
  }
  plotobj <- plotobj + ggtitle(info[1])
  plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)\n")
  plotobj <- plotobj + scale_x_continuous("\nTime after first dose (minutes)")
  plotobj
}

plotByFactor <- function(factorColname, factorText, dat) {
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
}
