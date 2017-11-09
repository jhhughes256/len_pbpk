# UI script for application for choosing standard curves
# -----------------------------------------------------------------------------

fluidPage(
  titlePanel("Standard Curve App"),
  sidebarLayout(
    sidebarPanel(
      # actionButton("console","Debug Console"),
      selectInput("numrun",
        "Select run date:",
        choices = list(
          "2nd October" = 1,
          "18th October" = 2
        )  # choices_datatype
      ),  # selectInput
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
      uiOutput("pointui")
    ),  # sidebarPanel
    mainPanel(
      plotOutput("stdcurve",
        height = "550px"
      ),  # plotOutput
      textOutput("r2"),
      br(),
      tableOutput("fbdf")
    )  # mainPanel
  )  # sidebarLayout
)  # fluidPage