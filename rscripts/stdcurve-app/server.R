# Server script for application for choosing standard curves
# -----------------------------------------------------------------------------

shinyServer(function(input, output, session) {
  Rdata <- reactive({
    subdata[subdata$Spreadsheet == input$datatype,]
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
  })  # reactive_Rdata
  
  Rlm <- reactive({
    lmres <- lm(Area.Ratio ~ Level, data = RdataStd(), weight = 1/Level**2)
  })
  
  Rline <- reactive({
    lmcoef <- Rlm()$coefficients
    level <- exp(seq(log(0.1), log(5000), by = 0.1))
    data.frame(
      Level = level,
      Area.Ratio = lmcoef["Level"]*level + lmcoef["(Intercept)"]
    )
  })
  
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
  })
  
  output$stdcurve <- renderPlot({
    Rplot()
  })  # renderPlot_stdcurve
  
  output$r2 <- renderText({
    paste("R-Squared =", summary(Rlm())$r.squared)
  })
  
  # Open console for R session
  # observe(label = "console", {
  #   if(input$console != 0) {
  #     options(browserNLdisabled = TRUE)
  #     saved_console <- ".RDuetConsole"
  #     if (file.exists(saved_console)) load(saved_console)
  #     isolate(browser())
  #     save(file = saved_console, list = ls(environment()))
  #   }
  # })
  
# Close the R session when browser closes
  session$onSessionEnded(function(){
   stopApp()
  })  # session.onSessionEnded
})  # shinyServer