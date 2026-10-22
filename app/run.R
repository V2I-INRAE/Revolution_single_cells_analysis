# Run from the project root: Rscript app/run.R
# The listener is local to this computer; no deployment configuration is needed.
shiny::runApp("app", host = "127.0.0.1", port = 3838,
  launch.browser = interactive())
