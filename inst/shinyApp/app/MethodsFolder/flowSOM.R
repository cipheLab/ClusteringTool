library(FlowSOM)
library(Biobase)

fct.parameters <- list("xgrid"=c(10,10,100,20),"ygrid"=c(10,10,100,20))

BRP_BM.flowSOM.execute <- function(fcs.file, params = list(10,10), markers_col)
{
    fcs.out.SOM <- ReadInput(fcs.file, transform = FALSE, scale = FALSE)
    
    
    x.d <- as.numeric(params[[1]])
    y.d <- as.numeric(params[[2]])
    if(x.d*y.d > nrow(fcs.file@exprs))
    {
        x.d <- as.integer(sqrt(nrow(fcs.file@exprs)))
        y.d <- as.integer(sqrt(nrow(fcs.file@exprs)))
    }
    
    fcs.out.SOM <- BuildSOM(fcs.out.SOM, colsToUse = as.numeric(markers_col), 
							xdim=x.d, 
							ydim=y.d)
    fcs.labels <- matrix(fcs.out.SOM$map$mapping[,1], ncol=1)
    colnames(fcs.labels) <- paste0("cluster_FlowSOM.",ncol(fcs.file@exprs)+1)

    return(fcs.labels)
}
