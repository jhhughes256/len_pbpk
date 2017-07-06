# ui.R for len_pbpk model exploration
# -----------------------------------------------------------------------------
header <- dashboardHeader(
  title = "Lenalidomide PBPK"
)  # dashboardHeader

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Model Simulation",
      tabName = "simtab"
    ),
    menuItem("Model Diagram",
      tabName = "diagtab"
    ),
    menuItem("Model Specifications",
      tabName = "spectab"
    )
  ),
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

simtab <- tabItem(tabName = "simtab",
  # Concentration-time plot
  fluidRow(
    box(width = 12, align = "center",
      uiOutput("meltHeaderUI"),
      h4("Simulated: Red | Observed: Blue"),
      plotOutput("meltPlot")
    )  # box.plot
  ),  # fluidRow.plot
  fluidRow(
    # Inputs
    box(width = 12, align = "center",
      h3("Model Parameters"),
      column(width = 4,
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
      column(width = 4,
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
      ),  # middlecolumn.inputs
      column(width = 4,
        numericInput("wt",
          "Body Weight (g):",
          value = 28
        ),  # numericInput.wt
        uiOutput("numericInputUI")
      )  # rightcolumn.inputs
    )  # box.inputs
  )  #fluidRow.inputs
)  # simtab

diagtab <- tabItem(tabName = "diagtab",
  box(width = 6, align = "center",
    img(src = "diag.png",
      width = "100%", height = "100%"  # so it scales with the window
    )  # img.diag
  ),  # box.img
  box(width = 6,
    withMathJax(),
    helpText("Some math here $$\\alpha+\\beta$$")
  )  # box.mathjax
)  # diagtab

spectab <- tabItem(tabName = "spectab",
  box(width = 12,
    h2("Model Specification File for Lenalidomide Mouse Model"),
    h4("For use with R package mrgsolve"),
    p("Sources for reference values:"),
    p(em("Brown RP, Delp MD, Lindstedt SL, Rhomberg LR, Beliles RP. Physiological
    Parameter Values for Physiologically Based Pharmacokinetic Models.
    Toxicology and Industrial Health. 1997;13(4):407-84.")),
    p(em("Davies B, Morris T. Physiological Parameters in Laboratory Animals and
    Humans. Pharmaceutical Research. 1993;10(7):1093-5.")),
    includeMarkdown("www/model.md")
  )  # box.md
)  # spectab

body <- dashboardBody(
  tabItems(simtab, diagtab, spectab)
)  # dashboardBody

dashboardPage(
  header, sidebar, body
)  # dashboardPage
