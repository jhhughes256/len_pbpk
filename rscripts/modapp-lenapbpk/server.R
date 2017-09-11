# server.R for len_pbpk model exploration
# ----------------------------------------------------------------------------
shinyServer(function(input, output, session) {
  Rconc <- reactive({
    # Set up sample times
    ID <- 1:11  # each ID represents a dose level
    ID2 <- sort(c(rep(ID, times = length(TIME))))
    times <- rep(TIME, times = length(ID))
    input.simdata <- data.frame(
      ID = ID2,
      time = times,
      amt = 0,
      evid = 0,
      rate = 0,
      cmt = 1,
      WT = input$wt,  # input here
      fu = input$fu,
      fuT = input$fuT,
      #PSsplstd = input$PSspl,
      khyd = input$khyd,
      PSdiffstd = input$PSdiff,
      Vmaxt = input$Vmax,
      kmt = input$km
    )  # input.simdata
    # Set up parameters dependent on input$comp
    if (!input$comp %in% c("PA", "ART")) {
      input.simdata$Q <- input[[paste0("Q", Rcomp()[1])]]
      input.simdata$V <- input[[paste0("V", Rcomp()[1])]]
      ncols <- dim(input.simdata)[2]
      names(input.simdata)[c(ncols-1, ncols)] <- c(
        paste0("Q", Rcomp()[1], "std"), paste0("V", Rcomp()[1], "std")
      )  # input.simdata.names
    }
    # Set up dose times
    dose.times <- 0
    dosedata <- input.simdata[input.simdata$time %in% dose.times, ]
    dosedata$amt <- c(
      c(0.5, 1.5, 5, 10)*unique(input.simdata$WT)*10^3,
      c(0.5, 10)*unique(input.simdata$WT)*10^3,
      c(0.5, 10)*unique(input.simdata$WT)*10^3,
      rep(1, 3)
    )
    # /10^3 for weight to kg, *10^6 for dose to ng -> *10^3
    # 1 for DOSENORM models
    dosedata$evid <- 1
    dosedata$rate <- dosedata$amt*120  # dose administered in 0.5 second
    dosedata$cmt <- c(rep(1, 4), rep(2, 2), rep(4, 2), c(1, 2, 4))
    # Combine dose times and sample times
    input.simdata <- rbind(input.simdata, dosedata)
    input.simdata <- input.simdata[with(input.simdata, order(ID, time)), ]
    # Simulate data
    output.simdata <- as.data.frame(mrgsim(
      data_set(mouse.mod, input.simdata)
    ))  # mrgsim
    output.simdata
  })  # reactive.Rconc

  Rmelt <- reactive({
    # Set up for data melting
    input.data <- Rconc()
    init.str <- names(mouse.mod@init)  # determine compartment names of model
    init.n <- length(init.str)  # determine number of compartments
    n.cols <- dim(input.data)[2]  # determine number of columns in Rconc()
    # Create dosemgkg column
    input.data$dosemgkg <- rep(c(0.5, 1.5, 5, 10, 0.5, 10, 0.5, 10, 1, 1, 1), each = length(TIME)+1)
    # Create type column (0 = PO, 1 = IV)
    input.data$data <- rep(c(1, 1, 1, 1, 0, 0, 2, 2, 1, 0, 2), each = length(TIME)+1)
    # Melt data for plotting
    output.data <- cbind(
      melt(  # ID, time, concentration columns, dosemgkg, iv
        input.data[c(1:2, (init.n+3):(2+init.n*2), n.cols+1, n.cols+2)],
        c("ID", "time", "dosemgkg", "data")  # id.vars - ?melt.data.frame
      )  # melt.output.data
    )  # output.data
    # Clean up melted data
    names(output.data) <- c("ID", "TIME", "DOSEMGKG", "DATA", "COMP", "C")
    output.data$COMP <- as.factor(output.data$COMP)
    levels(output.data$COMP) <- c(
      toupper(substr(init.str, 2, nchar(init.str)))
    )  # cleans up factor levels
    output.data <- output.data[output.data$TIME != 0, ]
    output.data
  })  # reactive.Rmelt

  Rcomp <- reactive({
    if (input$comp == "PA") {
      c("mix", "Plasma", 14, 1.2)
    } else if (input$comp == "ART") {
      c("lng", "Lung", 14, 0.18)
    } else if (input$comp == "BRA") {
      c("bra", "Brain", 0.46, 0.43)
    } else if (input$comp == "LVR") {
      c("lvr", "Liver", 2.3, 1.4)
    } else if (input$comp == "SPS") {
      c("spl", "Spleen", 0.16, 0.09)
    } else if (input$comp == "KID") {
      c("kid", "Kidney", 1.3, 0.43)
    } else if (input$comp == "HRT") {
      c("hrt", "Heart", 0.92, 0.13)
    } else if (input$comp == "MSC") {
      c("msc", "Muscle", 2.2, 9.6)
    }  # ifelse end
  })  # reactive.Rcomp

  output$meltHeaderUI <- renderUI({
    h2(paste(Rcomp()[2], "Concentration Time Profile"))
  })  # output.meltHeaderUI

  output$meltPlot <- renderPlot({
    # Subset data according to input
    Rplotdata <- Rmelt()[Rmelt()$DATA == input$route,]
    NRplotdata <- NRmelt[NRmelt$DATA == input$route,]
    if (input$dosenorm == F) {
      Rplotdata <- Rplotdata[Rplotdata$ID %in% c(1:8),]
    } else {
      Rplotdata <- Rplotdata[Rplotdata$ID %in% c(9:11),]  # dose normalised
    }
    if (input$route == 1) {
      Rplotdata <- Rplotdata[Rplotdata$COMP == input$comp,]
      NRplotdata <- NRplotdata[NRplotdata$COMP == input$comp,]
    } else {
      Rplotdata <- Rplotdata[Rplotdata$COMP == "PA",]
    }
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
    if (input$dosenorm == F) {
      plotobj <- plotobj + geom_point(aes(x = TIME, y = C), data = NRplotdata[NRplotdata$TYPE == 0,],
        colour = "blue", alpha = 0.5)
      plotobj <- plotobj + geom_line(aes(x = TIME, y = C), data = NRplotdata[NRplotdata$TYPE == 1,],
        colour = "blue", linetype = "dashed")
      # Set facet, axis titles and limit
      plotobj <- plotobj + facet_wrap(~DOSEMGKG, ncol = 4)
      yaxis.title <- "Concentrations (ng/mL)\n"
    } else {
      plotobj <- plotobj + geom_point(aes(x = TIME, y = C/(DOSEMG*10^6)),
        data = NRplotdata[NRplotdata$TYPE == 0,], colour = "blue", alpha = 0.5)
      plotobj <- plotobj + geom_smooth(aes(x = TIME, y = C/(DOSEMG*10^6)),
        data = NRplotdata[NRplotdata$TYPE == 1,], colour = "blue", linetype = "dashed")
      yaxis.title <- "Dose Normalised Concentrations (ng/mL/ng)\n"
    }
    plotobj <- plotobj + scale_x_continuous("\nTime (mins)")
    plotobj <- plotobj + scale_y_log10("Concentrations (ng/mL)\n",
      labels = scales::comma)
    plotobj <- plotobj + coord_cartesian(xlim = c(0, 300))
    plotobj
  })  # output.meltPlot

  output$numericInputUI <- renderUI({
    if (!input$comp %in% c("PA", "ART")) {
      div(
        numericInput(paste0("Q", Rcomp()[1]),
          paste(Rcomp()[2], "Blood Flow (ml/min):"),
          value = Rcomp()[3]
        ),  # numericInput.reactiveQ
        numericInput(paste0("V", Rcomp()[1]),
          paste(Rcomp()[2], "Volume (ml):"),
          value = Rcomp()[4]
        )  # numericInput.reactiveV
      )  # div.reactiveInput
    } else {  # if (input$comp %in% c("PA", "ART"))
      div(
        p("")
      )  # div.empty
    }  # ifelse end
  })  # output.numericInputUI

  # Close the R session when browser closes
  session$onSessionEnded(function(){
   stopApp()
  })  # session.onSessionEnded

  # Open debug console for R session
  # observe(label = "console", {
  #   if(input$console != 0) {
  #     options(browserNLdisabled = TRUE)
  #     # saved_console <- ".RDuetConsole"
  #     # if (file.exists(saved_console)) load(saved_console)
  #     isolate(browser())
  #     # save(file = saved_console, list = ls(environment()))
  #   }
  # })  # observe.console
})  # shinyServer
