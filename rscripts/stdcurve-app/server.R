# Server script for application for choosing standard curves
# -----------------------------------------------------------------------------

shinyServer(function(input, output, session) {
# Reactive UI
  Rui <- reactive({
    init <- initVals[[as.numeric(input$numrun)]][[as.numeric(input$datatype)]]
    div(
      fluidRowCol(1, init),
      fluidRowCol(2, init),
      fluidRowCol(3, init),
      fluidRowCol(4, init),
      fluidRowCol(5, init),
      fluidRowCol(6, init),
      fluidRowCol(7, init),
      fluidRowCol(8, init),
      fluidRowCol(9, init)
    )
  })
  
  output$pointui <- renderUI({
    Rui()
  })

# -----------------------------------------------------------------------------
# Reactive Plot
  
  Rdata <- reactive({
    subdata[subdata$Spreadsheet == input$datatype & subdata$Run == input$numrun,]
  })  # reactive_Rdata
  
  RdataStd <- reactive({
    dataStd <- Rdata()[Rdata()$Sample.Type == "Std Bracket Sample",]
    ddply(dataStd, .(Specified.Amount), function(x) {
      amt <- unique(x$Specified.Amount)
      if (get("input")[[paste0("status", which(concVec == amt))]]) {
        x[x$batch == get("input")[[paste0("set", which(concVec == amt))]], ]
      } else {
        x[x$batch == 0, ]
      }
    })
  })  # reactive_RdataStd
  
  Rlm <- reactive({
    lmres <- lm(Area.Ratio ~ Level, data = RdataStd(), weight = 1/Level**2)
  })  # reactive_Rlm
  
  Rline <- reactive({
    lmcoef <- Rlm()$coefficients
    level <- exp(seq(log(0.1), log(5000), by = 0.1))
    data.frame(
      Level = level,
      Area.Ratio = lmcoef["Level"]*level + lmcoef["(Intercept)"]
    )
  })  # reactive_Rline
  
  Rplot <- reactive({
    p <- NULL
    p <- ggplot()
    p <- p + geom_point(aes(x = Level, y = Area.Ratio), data = RdataStd(),
      colour = "blue", shape = 1, size = 3)
    p <- p + geom_line(aes(x = Level, y = Area.Ratio), data = Rline(),
      colour = "red", size = 0.8)
    p <- p + geom_point(aes(x = Level, y = Area.Ratio), 
      data = Rdata()[Rdata()$Sample.Type == "QC Sample",],
      colour = "black", size = 2)
    if (input$logx) p <- p + scale_x_log10()
    if (input$logy) p <- p + scale_y_log10()
    p
  })  # Rplot
  
  output$stdcurve <- renderPlot({
    Rplot()
  })  # renderPlot_stdcurve
  
  output$r2 <- renderText({
    paste("R-Squared =", summary(Rlm())$r.squared)
  })  # renderText_r2
  
# -----------------------------------------------------------------------------
# Reactive Table
  Rsamp <- reactive({
    lmcoef <- Rlm()$coefficients
    samp <- subdata[subdata$Sample.Type == "Unknown Sample" & subdata$Run == input$numrun,]
    
    samp$Species <- "mouse"
    samp$Species[str_detect(samp$Filename, "h")] <- "human"
    samp$Vehicle <- "plasma"
    samp$Vehicle[str_detect(samp$Filename, "b")] <- "buffer"
    samp$ID <- 4
    samp$ID[str_detect(samp$Filename, "5")] <- 5
    samp$ID[str_detect(samp$Filename, "6")] <- 6
    samp$ID[str_detect(samp$Filename, "7")] <- 7
    samp$ID[str_detect(samp$Filename, "8")] <- 8
    samp$ID[str_detect(samp$Filename, "9")] <- 9
    samp$Conc <- 3
    samp$Conc[str_detect(samp$Filename, "0_3")] <- 0.3
    samp$Conc[str_detect(samp$Filename, "_1_")] <- 1
    samp$Conc[str_detect(samp$Filename, "0_03")] <- 0.03
    samp$Conc[str_detect(samp$Filename, "10_")] <- 10
    
    samp$dv <- (samp$Area.Ratio - lmcoef["(Intercept)"])/lmcoef["Level"]
    samp$dv[samp$dv <= 0] <- NA
    samp
  })  # reative_Rsamp
  
  output$fbdf <- renderTable({
    ddply(Rsamp(), .(Species, Conc, ID), function(x) {
      Cp <- x$dv[x$Vehicle == "plasma"]  # total plasma conc
      Cd <- x$dv[x$Vehicle == "buffer"]  # free dialysate conc
      Ci <- unique(x$Conc)*1000  # initial plasma conc; convert from uM to nM
      fb <- (Cp-Cd)/Cp*100  # fraction bound
      dr <- Cd/Ci*100  # drug recovered in dialysate
      pr <- Cp/Ci*100  # drug recovered in plasma
      data.frame(
        fraction_bound = fb, 
        plas_recovered = pr, 
        dial_recovered = dr, 
        total_recovered = dr+pr
      )
    })
  })  # renderDataTable_fbdf

# -----------------------------------------------------------------------------
# Open console for R session
  # observe(label = "console", {
  #   if(input$console != 0) {
  #     options(browserNLdisabled = TRUE)
  #     isolate(browser())
  #   }
  # })
  
# Close the R session when browser closes
  session$onSessionEnded(function(){
   stopApp()
  })  # session.onSessionEnded
})  # shinyServer