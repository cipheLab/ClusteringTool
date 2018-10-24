library(stats)

fct.parameters <- list("centers"=c(0,10,500,100),"iterations"=c(0,100,1000,100), "nstart"=c(0,1,50,6))

BRP_BM.kmeans.execute <- function(fcs.file, params = list(50,10,1), markers_col)
{
    k <- min(as.numeric(params[[1]]), nrow(fcs.file@exprs))
    fcs.out.kmeans <- kmeans(fcs.file@exprs[,as.numeric(markers_col)],
							centers=k,
							iter.max=as.numeric(params[[2]]),
							nstart=as.numeric(params[[3]]))
    fcs.labels <- matrix(fcs.out.kmeans$cluster, ncol=1)
    colnames(fcs.labels) <- paste0("cluster_K-Means.",ncol(fcs.file@exprs)+1)

    return(fcs.labels)
}
