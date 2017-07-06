# ui.R for len_pbpk model exploration
# -----------------------------------------------------------------------------
header <- dashboardHeader(
  title = "Lenalidomide PBPK"
)  # dashboardHeader

sidebar <- dashboardSidebar(
  actionButton("console", "Debug"),
  selectInput("comp",
    "Compartment:",
    choices = list(
      "Plasma" = "PA",
      "Lungs" = "ART",
      "Brain" = "BRA",
      "Liver" = "LVR",
      "Spleen" = "SPS",
      "Kidney" = "KID",
      "Heart" = "HRT",
      "Muscle" = "MSC"
    ),  # choices.comp
    selected = "PA"
  )  # selectInput.comp
)  # dashboardSidebar

body <- dashboardBody(
  # Concentration-time plot
  fluidRow(
    box(width = 12, align = "center",
      uiOutput("meltPlotTitle"),
      h4("Simulated: Red | Observed: Blue"),
      plotOutput("meltPlot")
    )  # box.plot
  ),  # fluidRow.plot
  fluidRow(
    # Inputs
    box(width = 12, align = "center",
      h3("Model Parameters"),
      column(width = 6,
        numericInput("fu",
          "Fraction unbound in Plasma:",
          value = 0.6,
          step = 0.04
        ),  # numericInput.fu
        numericInput("PSspl",
          "PS Spleen (ml/min):",
          value = 0.18,
          step = 0.01
        ),  # numericInput.PSspl
        numericInput("Vmax",
          "Vmax of renal secretion:",
          value = 100
        )  # numericInput.Vmax
      ),  # leftcolumn.inputs
      column(width = 6,
        numericInput("fuT",
          "Fraction unbound in Tissue:",
          value = 0.48,
          step = 0.04
        ),  # numericInput.fuT
        numericInput("PSdiff",
          "PS Renal Tubular Cells (ml/min):",
          value = 0.5,
          step = 0.01
        ),  # numericInput.PSdiff
        numericInput("km",
          "km of renal secretion:",
          value = 50
        )  # numericInput.km
      )  # rightcolumn.inputs
    )  # box.inputs
  )  #fluidRow.inputs
)  # dashboardBody

dashboardPage(
  header, sidebar, body
)  # dashboardPage
