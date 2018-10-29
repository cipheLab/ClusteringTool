# Clustering Tool
Shiny app developped to use different clustering algorithms on FCS files. The files are enriched and can then be downloaded and used with the Analysis Tool.
	 
>[User manual :](doc/manual_clusteringtool.pdf)

## Requirements
  * software: R(Version 3.4.3 to 3.5), Rstudio(optional)
  * R packages: flowcore, microbenchmark, ncdfFlow, shiny, shinydashboard, shinyjs, doSNOW, cluster, parallel, ggcyto
  
## Quick installation guide

  1. Run the following command in R/RStudio:
```
install.packages(c("microbenchmark,, "shiny", "shinyjs", "shinydashboard","cluster","parallel","doSNOW"))
source("https://bioconductor.org/biocLite.R")
biocLite("ggcyto")
biocLite("flowCore")
biocLite("FlowSOM")
biocLite("ncdfFlow")
```
  >You may be asked to reload your environment, if so, accept.
  
  2. Run the next commands:
```
library("devtools")
install_github("isambens/ClusteringTool")
```

  
## Launching the shiny application

  1. Run the following commands in R/RStudio:
```
library("ClusteringTool")
ClusteringTool.run()
```  