# ui.R for len_pbpk model exploration
# -----------------------------------------------------------------------------
header <- dashboardHeader(
  title = "Lenalidomide PBPK"
)  # dashboardHeader

sidebar <- dashboardSidebar(
  sidebarMenu(id = "tabs",
    menuItem("Model Simulation",
      tabName = "simtab"
    ),  # menuItem.simtab
    menuItem("Model Diagram",
      tabName = "diagtab"
    ),  # menuItem.diagtab
    menuItem("Model Specifications",
      tabName = "spectab"
    )  # menuItem.spectab
  ),  # sidebarMenu.tabs
  # actionButton("console", "Debug"),
  conditionalPanel(condition = "input.tabs == 'simtab'",
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
  )  # conditionalPanel.comp
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
          value = 6
        )  # numericInput.Vmax
      ),  # leftcolumn.inputs
      column(width = 4,
        numericInput("fuT",
          "Fraction unbound in Tissue:",
          value = 0.48,
          step = 0.04
        ),  # numericInput.fuT
        numericInput("PSdiff",
          "PS Kidney (ml/min):",
          value = 0.5,
          step = 0.01
        ),  # numericInput.PSdiff
        numericInput("km",
          "km of renal secretion:",
          value = 0.1
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
  box(width = 6, style = "font-size: 25px",
    withMathJax(),
    p("$$\\begin{align}
     \\ \\frac{dA_{VASC}}{dt}=&
        Q_{BRA}\\frac{A_{BRA}}{V_{BRA}}
        +Q_{LVR}\\frac{A_{LVR}}{V_{LVR}}
        +Q_{KID}\\frac{A_{KID}}{V_{KID}}
        +Q_{HRT}\\frac{A_{HRT}}{V_{HRT}} \\\\
     \\ &+Q_{MSC}\\frac{A_{MSC}}{V_{MSC}}
        +Q_{BOD}\\frac{A_{BOD}}{V_{BOD}}
        -Q_{CO}\\frac{A_{VASC}}{V_{VASC}} \\\\
     \\ \\frac{dA_{LNG}}{dt}=&
        Q_{CO}(\\frac{A_{VASC}}{V_{VASC}}-\\frac{A_{LNG}}{V_{LNG}}) \\\\
     \\ \\frac{dA_{BRA}}{dt}=&
        Q_{BRA}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{BRA}}{V_{BRA}}) \\\\
     \\ \\frac{dA_{SPL}}{dt}=&
        Q_{SPL}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{SPL}}{V_{SPL}})
        +PS_{SPL}(A_{SPLP}-A_{SPL}) \\\\
     \\ \\frac{dA_{SPLP}}{dt}=&
        PS_{SPL}(A_{SPL}-A_{SPLP}) \\\\
     \\ \\frac{dA_{LVR}}{dt}=&
        Q_{LVR}(Q{_{SPL}+Q_{GIT}})\\frac{A_{LNG}}{V_{LNG}}
        -Q_{LVR}\\frac{A_{LVR}}{V_{LVR}}
        +Q_{SPL}\\frac{A_{SPL}}{V_{SPL}}
        +Q_{GIT}\\frac{A_{GIT}}{V_{GIT}} \\\\
     \\ \\frac{dA_{GIT}}{dt}=&
        Q_{GIT}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{GIT}}{V_{GIT}}) \\\\
     \\ \\frac{dA_{KID}}{dt}=&
        Q_{KID}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{KID}}{V_{KID}})
        -CL_{GFR}\\frac{A_{LNG}}{V_{LNG}}
        +PS_{KID}(A_{TUBC}-A_{KID}) \\\\
     \\ CL_{TRAN}=&
        \\frac{V_{max}}{k_{m}+\\frac{A_{KID}}{V_{KID}}} \\\\
     \\ \\frac{dA_{TUBC}}{dt}=&
        PS_{KID}(A_{KID}-A_{TUBC})
        -CL_{TRAN}A_{TUBC} \\\\
     \\ \\frac{dA_{FILT}}{dt}=&
        CL_{GFR}\\frac{A_{LNG}}{V_{LNG}}
        +CL_{TRAN}A_{TUBC}
        -k_{URINE}A_{FILT} \\\\
     \\ \\frac{dA_{URINE}}{dt}=&
        k_{URINE}A_{FILT} \\\\
     \\ \\frac{dA_{HRT}}{dt}=&
        Q_{HRT}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{HRT}}{V_{HRT}}) \\\\
     \\ \\frac{dA_{MSC}}{dt}=&
        Q_{MSC}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{MSC}}{V_{MSC}}) \\\\
     \\ \\frac{dA_{BOD}}{dt}=&
        Q_{BOD}(\\frac{A_{LNG}}{V_{LNG}}-\\frac{A_{BOD}}{V_{BOD}}) \\\\
     \\end{align}$$"
    )  # p.differentials
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
