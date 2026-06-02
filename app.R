# =========================================================
# Shiny App: Visualizing Distribution Behaviour & Inference
# Base R version (no car package required)
# Features:
# - Paired differences (X - Y)
# - Many experiments (sampling distributions)
# - Normal overlay (sample mean/sd vs theoretical)
# - Histogram + Q-Q plot with 95% pointwise CI envelope
# - Shapiro-Wilk
# - Confidence intervals and coverage
# - Wide selection of distributions
# =========================================================


# https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html
# Run the shinylive::export line to populate the docs folder 
# so that shinylive works from github
#shinylive::export(appdir = "../NormalDistributionApproximation/", destdir = "docs")
#httpuv::runStaticServer("docs/", port = 8008)


# -----------------------------
# Q-Q plot with 95% CI envelope
# (replicates car::qqPlot behaviour using base R)
# Uses the pointwise Beta-distribution-based envelope
# as described in Fox & Weisberg (car package methodology)
# -----------------------------
qqPlot_base <- function(x,
                        main  = "Q-Q Plot with 95% CI",
                        col   = "#007CB7",
                        col.lines = "red",
                        conf  = 0.95) {
  
  n    <- length(x)
  sx   <- sort(x)
  
  # Sample mean and sd used for theoretical quantiles
  mu   <- mean(x)
  sig  <- sd(x)
  
  # Theoretical normal quantiles for plotting positions
  probs <- ppoints(n)          # (1:n - 0.375) / (n + 0.25) — Blom formula
  qth   <- qnorm(probs, mu, sig)
  
  # Standardise to N(0,1) for the envelope calculation
  z     <- (sx - mu) / sig
  zq    <- qnorm(probs)        # standard normal quantiles
  
  # 95% pointwise envelope using order-statistic variance
  # Var(z_(i)) ≈ (p*(1-p)) / (n * phi(zq)^2)
  alpha <- 1 - conf
  zz    <- qnorm(1 - alpha / 2)
  se    <- (1 / dnorm(zq)) * sqrt(probs * (1 - probs) / n)
  lo    <- zq - zz * se
  hi    <- zq + zz * se
  
  # Back-transform envelope to original scale
  lo_orig <- lo * sig + mu
  hi_orig <- hi * sig + mu
  
  # Plot limits with a little padding
  xlim <- range(qth)
  ylim <- range(c(sx, lo_orig, hi_orig))
  pad  <- diff(ylim) * 0.05
  ylim <- ylim + c(-pad, pad)
  
  # Base plot
  plot(qth, sx,
       main  = main,
       xlab  = "Theoretical Quantiles",
       ylab  = "Sample Quantiles",
       pch   = 16,
       col   = col,
       xlim  = xlim,
       ylim  = ylim,
       las   = 1)
  
  # 45-degree reference line (through 1st & 3rd quartiles, like car)
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  tQ <- qnorm(c(0.25, 0.75), mu, sig)
  slope     <- (Q3 - Q1) / (tQ[2] - tQ[1])
  intercept <- Q1 - slope * tQ[1]
  abline(a = intercept, b = slope, col = col.lines, lwd = 2)
  
  # Confidence envelope
  lines(qth, lo_orig, lty = 2, col = col.lines, lwd = 1.5)
  lines(qth, hi_orig, lty = 2, col = col.lines, lwd = 1.5)
  
  invisible(list(x = qth, y = sx, lower = lo_orig, upper = hi_orig))
}


# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  
  # ── Button styling ───────────────────────────────────────────
  tags$head(
    tags$style(HTML("
      #resample {
        background-color: #007CB7;
        color: white;
        border: none;
        padding: 6px 12px;
        margin-top: 4px;
        margin-bottom: 10px;
        font-size: 14px;
      }
      #resample:hover {
        background-color: #005f8e;
        color: white;
      }
      #resample:active {
        transform: scale(0.97);
      }
    "))
  ),
  # ─────────────────────────────────────────────────────────────
  
  titlePanel("Distribution Behaviour & Sampling Inference"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Distribution selector
      selectInput("dist", "Underlying distribution:",
                  choices = c("Normal",
                              "t (df=3)",
                              "Log-normal",
                              "Uniform",
                              "Exponential (skewed)",
                              "Chi-square (df=3)",
                              "Gamma (shape=2)",
                              "Beta (2,5)",
                              "Poisson (lambda=3)",
                              "Binomial (n=20, p=0.3)"),
                  selected = "Normal"),
      
      # ── Resample button ──────────────────────────────────────
      actionButton("resample", "Resample"),
      # ─────────────────────────────────────────────────────────
      
      # Sample size
      sliderInput("n", "Sample size per experiment:",
                  min = 5, max = 200, value = 20),
      
      # Paired differences toggle
      checkboxInput("paired_diff",
                    "Use paired differences (X - Y)", TRUE),
      
      # Many experiments toggle
      checkboxInput("many_exp",
                    "Simulate many experiments (to demonstrate CLT)",
                    TRUE),
      
      # Number of experiments
      sliderInput("n_exp",
                  "Number of experiments (only used with Simulation):",
                  min = 100, max = 5000, value = 1000, step = 100),
      
      # Standardize output
      checkboxInput("standardize",
                    "Standardize output (click to recalculate)", TRUE),
      
      # Normal overlay type
      checkboxInput("use_theoretical",
                    "Use theoretical normal (instead of sample mean/sd)",
                    FALSE),
      
      helpText("_______________________________"),
      helpText("Glenn Tattersall, PhD"),
      helpText("For use in BIOL 3P96 - Biostatistics")
    ),
    
    mainPanel(
      plotOutput("plot", height = "500px"),
      verbatimTextOutput("tests"),
      verbatimTextOutput("explain")
    )
  )
)


# -----------------------------
# SERVER
# -----------------------------
server <- function(input, output) {
 
   # -----------------------------
  # Function to generate a sample
  # -----------------------------
  generate_sample <- function(n, dist) {
    switch(dist,
           "Normal"                  = rnorm(n),
           "Uniform"                 = runif(n),
           "t (df=3)"                = rt(n, df = 3),
           "Exponential (skewed)"    = rexp(n, 1),
           "Log-normal"              = rlnorm(n, 0, 1),
           "Chi-square (df=3)"       = rchisq(n, 3),
           "Gamma (shape=2)"         = rgamma(n, shape = 2, rate = 1),
           "Beta (2,5)"              = rbeta(n, 2, 5),
           "Poisson (lambda=3)"      = rpois(n, 3),
           "Binomial (n=20, p=0.3)"  = rbinom(n, 20, 0.3)
    )
  }
  
  # -----------------------------
  # Reactive: Generate values
  # Depends on input$resample so clicking the button
  # invalidates this reactive and triggers a fresh draw.
  # -----------------------------
  get_values <- reactive({
    input$resample   # take a dependency on the button
    
    n     <- input$n
    n_exp <- input$n_exp
    
    if (!input$many_exp) {
      
      # Single experiment
      X      <- generate_sample(n, input$dist)
      Y      <- generate_sample(n, input$dist)
      values <- if (input$paired_diff) X - Y else X
      
    } else {
      
      # Many experiments → sampling distribution
      values <- replicate(n_exp, {
        X <- generate_sample(n, input$dist)
        Y <- generate_sample(n, input$dist)
        if (input$paired_diff) mean(X - Y) else mean(X)
      })
    }
    
    # Standardize if requested
    if (input$standardize) {
      values <- (values - mean(values)) / sd(values)
    }
    
    values
  })
  
  # -----------------------------
  # Reactive: Confidence intervals (many experiments only)
  # -----------------------------
  get_results <- reactive({
    if (!input$many_exp) return(NULL)
    
    n     <- input$n
    n_exp <- input$n_exp
    
    results <- replicate(n_exp, {
      
      X <- generate_sample(n, input$dist)
      Y <- generate_sample(n, input$dist)
      
      if (input$paired_diff) {
        vals      <- X - Y
        true_mean <- 0
      } else {
        vals      <- X
        true_mean <- switch(input$dist,
                            "Exponential (skewed)"   = 1,
                            "Uniform"                = 0.5,
                            "Log-normal"             = exp(0 + 0.5),
                            "Chi-square (df=3)"      = 3,
                            "Gamma (shape=2)"        = 2,
                            "Normal"                 = 0,
                            "t (df=3)"               = 0,
                            "Beta (2,5)"             = 2 / 7,
                            "Poisson (lambda=3)"     = 3,
                            "Binomial (n=20, p=0.3)" = 6,
                            NA)
      }
      
      m     <- mean(vals)
      s     <- sd(vals)
      se    <- s / sqrt(n)
      lower <- m - 1.96 * se
      upper <- m + 1.96 * se
      
      c(mean    = m,
        lower   = lower,
        upper   = upper,
        true    = true_mean,
        covered = as.numeric(true_mean >= lower & true_mean <= upper),
        width   = upper - lower)
    })
    
    as.data.frame(t(results))
  })
  
  # -----------------------------
  # Reactive: Shapiro-Wilk test (shared by tests + explain outputs
  # so both always refer to the exact same sample)
  # -----------------------------
  sw_result <- reactive({
    values <- get_values()
    if (length(values) <= 5000) shapiro.test(values) else NULL
  })
  
  # -----------------------------
  # Plot: Histogram + Q-Q plot
  # -----------------------------
  output$plot <- renderPlot({
    
    values <- get_values()
    n      <- input$n
    
    # Determine normal overlay parameters
    if (input$use_theoretical) {
      
      if (!input$many_exp) {
        if (input$paired_diff) {
          mu    <- 0
          sigma <- sqrt(2)
        } else {
          mu    <- switch(input$dist,
                          "Exponential (skewed)"   = 1,
                          "Log-normal"             = exp(0 + 0.5),
                          "Chi-square (df=3)"      = 3,
                          "Gamma (shape=2)"        = 2,
                          "Uniform"                = 0.5,
                          "Normal"                 = 0,
                          "t (df=3)"               = 0,
                          "Beta (2,5)"             = 2 / 7,
                          "Poisson (lambda=3)"     = 3,
                          "Binomial (n=20, p=0.3)" = 6)
          sigma <- sd(values)
        }
      } else {
        if (input$paired_diff) {
          mu    <- 0
          sigma <- sqrt(2 / n)
        } else {
          mu    <- switch(input$dist,
                          "Exponential (skewed)"   = 1,
                          "Log-normal"             = exp(0 + 0.5),
                          "Chi-square (df=3)"      = 3,
                          "Gamma (shape=2)"        = 2,
                          "Uniform"                = 0.5,
                          "Normal"                 = 0,
                          "t (df=3)"               = 0,
                          "Beta (2,5)"             = 2 / 7,
                          "Poisson (lambda=3)"     = 3,
                          "Binomial (n=20, p=0.3)" = 6)
          sigma <- sd(values)
        }
      }
      
    } else {
      mu    <- mean(values)
      sigma <- sd(values)
    }
    
    # Layout: 1 row, 2 columns
    par(mfrow = c(1, 2))
    
    # --- Histogram ---
    hist(values,
         probability = TRUE,
         breaks      = 30,
         col         = "#007CB7",
         border      = "#263056",
         main        = "Histogram with Normal Overlay",
         xlab        = "Value")
    curve(dnorm(x, mean = mu, sd = sigma),
          col = "red", lwd = 2, add = TRUE)
    
    # True mean line in many-experiments mode
    if (input$many_exp) {
      res <- get_results()
      abline(v = unique(res$true), col = "blue", lwd = 2, lty = 2)
    }
    
    # --- Q-Q plot with 95% envelope (base R, no car) ---
    qqPlot_base(values,
                main      = "Q-Q Plot with 95% CI",
                col       = "#007CB7",
                col.lines = "red",
                conf      = 0.95)
  })
  
  # -----------------------------
  # Normality tests
  output$tests <- renderPrint({
    result <- sw_result()
    
    subject <- if (!input$many_exp && !input$paired_diff) {
      "the raw sampled values (single experiment)"
    } else if (!input$many_exp && input$paired_diff) {
      "the paired differences X - Y (single experiment)"
    } else if (input$many_exp && !input$paired_diff) {
      paste0("the sampling distribution of means (", input$n_exp, " simulated experiments)")
    } else {
      paste0("the sampling distribution of mean differences (", input$n_exp, " simulated experiments)")
    }
    
    if (!is.null(result)) {
      cat("Shapiro-Wilk Test on", subject, ":\n")
      print(result)
    } else {
      cat("Shapiro-Wilk Test: skipped (n > 5000)\n")
    }
  })
  
  # -----------------------------
  # Teaching explanation text
  # -----------------------------
  output$explain <- renderText({
    result <- sw_result()
    
    wrap80 <- function(...) paste(strwrap(paste0(...), width = 80), collapse = "\n")
    
    if (is.null(result)) {
      return("Sample too large for Shapiro-Wilk; use the Q-Q plot to assess normality.")
    }
    
    pvalueSW   <- result$p.value
    sig        <- pvalueSW < 0.05
    dist_label <- input$dist
    
    if (!input$many_exp && !input$paired_diff) {
      base <- wrap80(
        "What you're seeing: A single sample of ", input$n, " values drawn ",
        "directly from the ", dist_label, " distribution. ",
        "The Shapiro-Wilk test is asking: do these raw values look like they ",
        "came from a normal distribution?"
      )
      note <- if (sig) {
        paste0("\n\n", wrap80(
          "Result: The test suggests non-normality (p < 0.05). This reflects ",
          "the shape of the underlying distribution itself."
        ))
      } else {
        paste0("\n\n", wrap80(
          "Result: No strong evidence of non-normality (p >= 0.05). ",
          "The sample is consistent with normality."
        ))
      }
      
    } else if (!input$many_exp && input$paired_diff) {
      base <- wrap80(
        "What you're seeing: Paired differences (X - Y) from a single ",
        "experiment (n = ", input$n, "), where both X and Y are drawn from ",
        "the ", dist_label, " distribution. ",
        "The Shapiro-Wilk test is asking: do these differences look normally ",
        "distributed?"
      )
      note <- if (sig) {
        paste0("\n\n", wrap80(
          "Result: The differences appear non-normal (p < 0.05). Even after ",
          "differencing, the shape of the original distribution leaves a ",
          "detectable skew or irregularity."
        ))
      } else {
        paste0("\n\n", wrap80(
          "Result: The differences look approximately normal (p >= 0.05). ",
          "Subtracting two samples from the same distribution often produces ",
          "a more symmetric, bell-shaped result."
        ))
      }
      
    } else if (input$many_exp && !input$paired_diff) {
      base <- wrap80(
        "What you're seeing: The sampling distribution of the mean - ",
        input$n_exp, " simulated experiments, each drawing n = ", input$n,
        " values from the ", dist_label, " distribution, with the mean ",
        "computed each time. The Shapiro-Wilk test is asking: does this ",
        "collection of means look normally distributed? This is a direct ",
        "demonstration of the Central Limit Theorem (CLT)."
      )
      note <- if (sig) {
        paste0("\n\n", wrap80(
          "Result: The distribution of means is still detectably non-normal ",
          "(p < 0.05). With a small n or a heavily skewed distribution, the ",
          "CLT may need a larger sample size to fully kick in."
        ))
      } else {
        paste0("\n\n", wrap80(
          "Result: The distribution of means looks normal (p >= 0.05). ",
          "This is the CLT in action: even if the underlying distribution is ",
          "not normal, means of repeated samples tend toward normality."
        ))
      }
      
    } else {
      base <- wrap80(
        "What you're seeing: The sampling distribution of mean paired ",
        "differences - ", input$n_exp, " simulated experiments, each ",
        "computing the mean of (X - Y) with n = ", input$n, " pairs drawn ",
        "from the ", dist_label, " distribution. The Shapiro-Wilk test is ",
        "asking: does this collection of mean differences look normally ",
        "distributed?"
      )
      note <- if (sig) {
        paste0("\n\n", wrap80(
          "Result: The mean differences are still detectably non-normal ",
          "(p < 0.05). Try increasing the sample size to see the CLT pull ",
          "the distribution toward normality."
        ))
      } else {
        paste0("\n\n", wrap80(
          "Result: The mean differences look normal (p >= 0.05). Averaging ",
          "over many paired observations smooths out the original ",
          "distribution's shape, consistent with the CLT."
        ))
      }
    }
    
    paste0(base, note)
  })
}

# -----------------------------
# Run the app
# -----------------------------
shinyApp(ui = ui, server = server)