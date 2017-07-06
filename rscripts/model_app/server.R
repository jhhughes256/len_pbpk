# server.R for len_pbpk model exploration
# ----------------------------------------------------------------------------
shinyServer(function(input, output, session) {
  Rconc <- reactive({
    # Set up sample times
    ID <- 1:4  # each ID represents a dose level
    ID2 <- sort(c(rep(ID, times = length(TIME))))
    times <- rep(TIME, times = length(ID))
    input.simdata <- data.frame(
      ID = ID2,
      time = times,
      amt = 0,
      evid = 0,
      rate = 0,
      cmt = 1,
      WT = 28,  # input here
      fu = input$fu,
      fuT = input$fuT,
      PSsplstd = input$PSspl,
      PSdiffstd = input$PSdiff,
      Vmaxt = input$Vmax,
      kmt = input$km
    )
    # Set up dose times
    dose.times <- 0
    dosedata <- input.simdata[input.simdata$time %in% dose.times, ]
    dosedata$amt <- c(0.5, 1.5, 5, 10)*unique(input.simdata$WT)*10^3
    # /10^3 for weight to kg, *10^6 for dose to ng -> *10^3
    dosedata$evid <- 1
    dosedata$rate <- dosedata$amt*60  # dose administered in 1 second
    dosedata$cmt <- 1
    # Combine dose times and sample times
    input.simdata <- rbind(input.simdata, dosedata)
    input.simdata <- input.simdata[with(input.simdata, order(ID, time)), ]
    # Simulate data
    output.simdata <- as.data.frame(mrgsim(
      data_set(mouse.mod, input.simdata)
    ))  # mrgsim
    output.simdata
  })  # Rconc

  Rmelt <- reactive({
    # Set up for data melting
    input.data <- Rconc()
    init.str <- names(mouse.mod@init)  # determine compartment names of model
    init.n <- length(init.str)  # determine number of compartments
    n.cols <- dim(input.data)[2]  # determine number of columns in Rsim()
    # Create dosemgkg column
    dosemgkg <- factor(input.data$ID)
    levels(dosemgkg) <- c(0.5, 1.5, 5, 10)
    input.data$dosemgkg <- as.numeric(levels(dosemgkg))[dosemgkg]
    # Melt data for plotting
    output.data <- cbind(
      melt(
        input.data[c(1:2, (init.n+3):(2+init.n*2), n.cols+1)],  # conc columns
        c("ID", "time", "dosemgkg")  # id.vars - ?melt.data.frame
      )
    )
    # Clean up melted data
    names(output.data) <- c("ID", "TIME", "DOSEMGKG", "COMP", "C")
    output.data$COMP <- as.factor(output.data$COMP)
    levels(output.data$COMP) <- c(
      toupper(substr(init.str, 2, nchar(init.str)))
    )  # cleans up factor levels
    output.data <- output.data[output.data$TIME != 0, ]
    output.data
  })  # Rmelt

  output$meltPlotTitle <- renderUI({
    if (input$comp == "PA") {
      h2("Plasma Concentration Time Profile")
    } else if (input$comp == "ART") {
      h2("Lung Concentration Time Profile")
    } else if (input$comp == "BRA") {
      h2("Brain Concentration Time Profile")
    } else if (input$comp == "LVR") {
      h2("Liver Concentration Time Profile")
    } else if (input$comp == "SPS") {
      h2("Spleen Concentration Time Profile")
    } else if (input$comp == "KID") {
      h2("Kidney Concentration Time Profile")
    } else if (input$comp == "HRT") {
      h2("Heart Concentration Time Profile")
    } else if (input$comp == "MSC") {
      h2("Muscle Concentration Time Profile")
    }
  })

  output$meltPlot <- renderPlot({
    # Subset data according to input
    Rplotdata <- Rmelt()[Rmelt()$COMP == input$comp,]
    NRplotdata <- NRmelt[NRmelt$COMP == input$comp,]
    # Add mg/kg to each facet label
    Rplotdata$DOSEMGKG <- factor(Rplotdata$DOSEMGKG)
    levels(Rplotdata$DOSEMGKG) <- paste(levels(Rplotdata$DOSEMGKG), "mg/kg")
    NRplotdata$DOSEMGKG <- factor(NRplotdata$DOSEMGKG)
    levels(NRplotdata$DOSEMGKG) <- paste(levels(NRplotdata$DOSEMGKG), "mg/kg")
    # Ready plot object
    plotobj <- NULL
    plotobj <- ggplot()
    # Plot data
    plotobj <- plotobj + geom_line(aes(x = TIME, y = C), data = Rplotdata,
      colour = "red", size = 1)
    plotobj <- plotobj + geom_point(aes(x = TIME, y = C), data = NRplotdata[NRplotdata$TYPE == 0,],
      colour = "blue", alpha = 0.5)
    plotobj <- plotobj + geom_line(aes(x = TIME, y = C), data = NRplotdata[NRplotdata$TYPE == 1,],
      colour = "blue", linetype = "dashed")
    # Set facet, axis titles and limit
    plotobj <- plotobj + facet_wrap(~DOSEMGKG, ncol = 4)
    plotobj <- plotobj + scale_x_continuous("\nTime (mins)")
    plotobj <- plotobj + scale_y_log10("Concentrations (ng/mL)\n",
      labels = scales::comma)
    plotobj <- plotobj + coord_cartesian(xlim = c(0, 100))
    plotobj
  })  # output.meltPlot

  # Close the R session when browser closes
  session$onSessionEnded(function(){
   stopApp()
  })  # onSessionEnded

  # Open debug console for R session
  observe(label = "console", {
    if(input$console != 0) {
      options(browserNLdisabled = TRUE)
      # saved_console <- ".RDuetConsole"
      # if (file.exists(saved_console)) load(saved_console)
      isolate(browser())
      # save(file = saved_console, list = ls(environment()))
    }
  })  # observe.console
})  # shinyServer
