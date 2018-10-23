# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @export
ClusteringTool.run <- function()
{
    library(flowCore)
    library(microbenchmark)
    library(ncdfFlow)
    library(shiny)
    library(shinydashboard)
    library(shinyjs)
    library(FlowSOM)
    library(cluster)
    library(parallel)
    library(doSNOW)
    library(Biobase)


    appDir <- system.file("shinyApp", "app", package = "ClusteringTool")
    if (appDir == "")
    {
        stop("Could not find app directory. Try re-installing `ClusteringTool`.", call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal", launch.browser = T)
}
