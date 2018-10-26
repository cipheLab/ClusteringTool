library(flowCore)

fct.parameters <- list("k"=c(1,10,100,10))

BRP_BM.hclust.execute <- function(fcs.file, params = list(10,40000), markers_col)
{
    nmb.events <- min(as.numeric(params[[2]]),nrow(fcs.file@exprs))
    k <- min(nrow(fcs.file@exprs), 40000)
    fcs.labels <- cutree(hclust(dist(fcs.file@exprs[,as.numeric(markers_col)])),k)
    
    fcs.labels <- matrix(fcs.labels, ncol=1)
    colnames(fcs.labels) <- paste0("cluster_hclust.",ncol(fcs.file@exprs)+1)
    
    return(fcs.labels)
}
