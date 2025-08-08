## pdac model shiny page – single adjusted survival curve
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  library(p, character.only = TRUE)
}
usePackage("shiny")
usePackage("shinydashboard")
usePackage("ggplot2")
usePackage("plotly")
usePackage("shinycssloaders")
usePackage("psc")

# grab internal helper from the psc package
modp <- psc:::modp

# load fitted counterfactual model
load("M:/RJ Fellowship/PDAC model/CFM.Rds")
smod <- CFM

# survival function from model + covariates
spline_surv_est <- function(lam, kn, k, haz_co, cov_co, cov = NULL, maxTime = 60, beta = 0) {
  if (is.null(cov)) cov <- rep(0, length(cov_co))
  tm   <- seq(0.01, maxTime, length = 500)
  logt <- log(tm)
  lp   <- as.numeric(cov %*% cov_co)

  z   <- NULL
  z_h <- NULL
  for (i in seq_len(k)) {
    z <- cbind(z,
      modp(logt - kn[i+1])^3
        - lam[i+1] * modp(logt - kn[1])^3
        - (1 - lam[i+1]) * modp(logt - kn[length(kn)])^3
    )
    z_h <- cbind(z_h,
      modp(logt - kn[i+1])^2
        - lam[i+1] * modp(logt - kn[1])^2
        - (1 - lam[i+1]) * modp(logt - kn[length(kn)])^2
    )
  }

  H0 <- exp(haz_co[1] + haz_co[2] * logt + z %*% haz_co[3:(2 + k)])
  h0 <- (H0 / tm) * (haz_co[2] + 3 * (z_h %*% haz_co[3:(2 + k)]))

  H <- H0 * exp(lp + beta)
  h <- h0 * exp(lp + beta)
  S <- exp(-H)
  f <- S * h

  data.frame(time = tm, S = S, f = f)
}

# wrapper to make calling the model easier
cfmEst <- function(maxTime = 60, beta = 0, cov = NULL) {
  spline_surv_est(
    lam    = smod$lam,
    kn     = smod$kn,
    k      = smod$k,
    haz_co = smod$haz_co,
    cov_co = smod$cov_co,
    cov    = cov,
    maxTime = maxTime,
    beta    = beta
  )
}

# ui
ui <- dashboardPage(
  dashboardHeader(title = "cfm survival explorer – PDAC model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Survival curve", tabName = "surv_tab", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "surv_tab",
        fluidRow(
          box(
            title = "Inputs", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 4,
            selectInput("t_stage", "T-stage:", choices = c(2,3,4), selected = 2),
            selectInput("grade",   "Grade:",   choices = c(1,2,3), selected = 2),
            selectInput("nodes",   "Nodes:",   choices = c(1,2),   selected = 1),
            numericInput("lca199", "LCA199:", value = 1, min = 0, step = 1),
            sliderInput("maxTime", "Max time (months):", min = 12, max = 120, value = 60, step = 12)
          ),
          box(
            title = "Survival plot", status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 8,
            withSpinner(plotlyOutput("survPlot", height = "500px"))
          )
        )
      )
    )
  )
)

# server
server <- function(input, output, session) {
  # build covariate matrix from inputs
  newCov <- reactive({
    df <- data.frame(
      t      = factor(input$t_stage, levels = c(2,3,4)),
      grade  = factor(input$grade,   levels = c(1,2,3)),
      nodes  = factor(input$nodes,   levels = c(1,2)),
      lca199 = input$lca199
    )
    model.matrix(~ t + grade + nodes + lca199, data = df)[, -1, drop = FALSE]
  })

  # get adjusted survival curve only
  survData <- reactive({
    cfmEst(maxTime = input$maxTime, cov = newCov())
  })

  # render Plotly
  output$survPlot <- renderPlotly({
    sd <- survData()
    df <- data.frame(time = sd$time, S = sd$S)

    p <- ggplot(df, aes(x = time, y = S)) +
      geom_line(size = 1.2) +
      labs(x = "Time (months)", y = "Survival probability") +
      theme_minimal()

    ggplotly(p) %>%
    layout(
      yaxis = list(
        range      = c(0, 1),      # lock to [0,1]
        autorange  = FALSE,        # disable auto-scaling
        fixedrange = TRUE          # prevent any y-zoom/pan
      ),
      xaxis = list(
        fixedrange = FALSE         # still allow zooming on x if you like
      )
    )
  })
}

# launch app
shinyApp(ui, server)
