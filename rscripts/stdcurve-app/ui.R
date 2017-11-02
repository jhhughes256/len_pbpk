# UI script for application for choosing standard curves
# -----------------------------------------------------------------------------

fluidPage(
  titlePanel("Standard Curve App"),
  sidebarLayout(
    sidebarPanel(
      # actionButton("console","Debug Console"),
      selectInput("datatype",
        "Select dataset:",
        choices = list(
          "Mouse Plasma" = 1,
          "Human Plasma" = 2,
          "PBS Buffer" = 3
        )  # choices_datatype
      ),  # selectInput
      checkboxInput("logx",
        "x-axis log-scale",
        value = TRUE
      ),  # checkboxInput_logx
      checkboxInput("logy",
        "y-axis log-scale",
        value = TRUE
      ),  # checkboxInput_logy
      fluidRowCol(1),
      fluidRowCol(2),
      fluidRowCol(3),
      fluidRowCol(4),
      fluidRowCol(5),
      fluidRowCol(6),
      fluidRowCol(7),
      fluidRowCol(8),
      fluidRowCol(9)
      # uiOutput("")
    ),  # sidebarPanel
    mainPanel(
      plotOutput("stdcurve",
        height = "550px"
      ),  # plotOutput
      textOutput("r2")
    )  # mainPanel
  )  # sidebarLayout
)  # fluidPage