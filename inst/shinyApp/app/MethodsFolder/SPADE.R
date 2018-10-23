library(cluster)
library(Biobase)
library(spadeCIPHE)

fct.parameters <- list("k"=c(0,100,1000,200), "apprx_mult"=c(0.1,0.1,3,1.5), "kernel_mult"=c(0.1,0.1,10,5), 
                       "exclude_pctile"=c(0,0.01,1,0.01), "target_pctile"=c(0,0.01,1,0.05))

BRP_BM.SPADE.execute <- function(fcs.file, params = list(200,1.5,5,0.01,0.05), markers_col)
{
    x <- fcs.file@exprs[,-cluster.col]
    
    #ADD DENSITY============================================================================================================================================
    fcs.dens <- SPADE.addDensityToFCSCIPHE(fcs.file, cols = markers_col, transforms = flowCore::linearTransform(), comp = F,
                                           apprx_mult = as.numeric(params[[2]]), kernel_mult = as.numeric(params[[3]]))
    
    #DOWNSAMPLE=============================================================================================================================================
    fcs.ds <- SPADE.downsampleFCSCIPHE(fcs.dens, exclude_pctile=as.numeric(params[[4]]), target_pctile = as.numeric(params[[5]]))
    
    #FCSToTree==============================================================================================================================================
    fcs.tree <- SPADE.FCSToTreeCIPHE(fcs.ds, cols = markers_col, k = as.numeric(params[[1]]), transforms = flowCore::linearTransform(), comp = F)
    
    #UPSAMPLE===============================================================================================================================================
    fcs.up <- SPADE.addClusterToFCSCIPHE(fcs.dens, fcs.tree[[1]], cols = markers_col, transforms = flowCore::linearTransform(), comp = F)
    
    #GENERATE FCS===========================================================================================================================================
    fcs.labels <- matrix(fcs.up@exprs[,ncol(fcs.up@exprs)],ncol=1)
    colnames(fcs.labels) <- paste0("cluster_SPADE",ncol(fcs.file@exprs)+1)
    
    return(fcs.labels)
}
