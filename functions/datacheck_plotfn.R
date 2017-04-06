# Plot functions used in datacheck_nonclin.R

iv.distplot <- function(input, output) {
  plotobj <- ggplot(input, aes(log(DV)))
  plotobj <- plotobj + geom_histogram()
  plotobj <- plotobj + facet_wrap(~TISSUE, ncol = 4)
  filename.out <- paste(
    output.dir,
    paste("Histogram_DVlog", output, sep = "_"),
    sep = "/"
  )  # paste
  plotobj
  ggsave(filename.out)
}

iv.CvTplot <- function(dat, dosenorm = F) {
  plotobj <- NULL
  plotobj <- ggplot(data = dat)
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
    paste0("CvTplot_", info[2], ".png"),
    sep = "/"
  )  # paste
  plotobj
  ggsave(filename.out, width = as.numeric(info[3]), height = as.numeric(info[4]))
}
