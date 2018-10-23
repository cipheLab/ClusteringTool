library(flowCore)

'%notin%' <- Negate('%in%')

get.params.list <- function(runs)
{
    params.names.list <- c()
    "%notin%" <- Negate("%in%")
    ordered.table <- c()
    
    #INIT LISTS
    lapply(1:length(runs), function(current.run.id)
    {
        lapply(1:length(runs[[current.run.id]]), function(p)
        {
            if(names(runs[[current.run.id]])[p]%notin%c(params.names.list))
            {
                params.names.list <<- c(params.names.list, names(runs[[current.run.id]])[p])
            }
        })
    })
    
    return(params.names.list)
}

get.ordered.table.from.runs <- function(runs, param)#Runs = list of list of parameters only (list 1 = run 1, ...)
{
    ordered.param.ids <- c()
    ordered.param.values <- c()
    irrelevant.param.ids <- c()
    params.names.list <- c()
    "%notin%" <- Negate("%in%")
    ordered.table <- c()
    
    #INIT LISTS
    lapply(1:length(runs), function(current.run.id)
    {
        par.id <- which(names(runs[[current.run.id]])%in%param)
        if(length(par.id)>0)
        {
            ordered.param.ids <<- c(ordered.param.ids, current.run.id)
            ordered.param.values <<- c(ordered.param.values, runs[[current.run.id]][[par.id[[1]]]])
            
        }
        else
        {
            irrelevant.param.ids <<- c(irrelevant.param.ids, current.run.id)
        }
        
        lapply(1:length(runs[[current.run.id]]), function(p)
        {
            if(names(runs[[current.run.id]])[p]%notin%c(params.names.list,param))
            {
                params.names.list <<- c(params.names.list, names(runs[[current.run.id]])[p])
            }
        })
    })
    
    
    #ORDER LISTS BY PARAM
    ordered.param.ids <- ordered.param.ids[order(ordered.param.values)]
    ordered.param.values <- ordered.param.values[order(ordered.param.values)]
    ordered.param.values <- split(ordered.param.values,factor(ordered.param.values))
    t <- lapply(1:length(ordered.param.values), function(i)
    {
        v <- ordered.param.ids[1:length(ordered.param.values[[i]])]
        ordered.param.ids <<- ordered.param.ids[-c(1:length(ordered.param.values[[i]]))]
        return(unlist(v))
    })
    ordered.param.ids <- t


    #CREATE OTHER PARAMS CODES
    temp.codes <- c()
    lapply(1:length(runs), function(current.run.id)
    {
        t <- c()
        lapply(1:length(params.names.list), function(p)
        {
            p.val.id <- which(names(runs[[current.run.id]])==params.names.list[[p]])
            if(length(p.val.id)>0)
            {
                val <- runs[[current.run.id]][[p.val.id[[1]]]]
                t <<- c(t, val)
            }
            else
            {
                t <<- c(t, "NULL")
            }
        })
        temp.codes <<- rbind(temp.codes, t)
    })
    #print(temp.codes)
    colnames(temp.codes) <- params.names.list

    #CREATE TABLE
    ordered.table.1 <- NULL
    ordered.table.2 <- NULL
    t <- unlist(ordered.param.ids)
    if(length(t)>0)
    {
        ordered.table.1 <- cbind(unlist(ordered.param.ids), unlist(ordered.param.values), t(temp.codes[unlist(ordered.param.ids),]))
    }
    t <- unlist(irrelevant.param.ids)
    if(length(t)>0)
    {
        ordered.table.2 <- cbind(unlist(irrelevant.param.ids), "NULL", t(temp.codes[unlist(irrelevant.param.ids),]))
    }
    ordered.table <- cbind(ordered.table.1, ordered.table.2)
    rownames(ordered.table) <- rep(NULL,nrow(ordered.table))
    colnames(ordered.table)[1:2] <- c("id",paste0(param))
    
    
    
    
    return(ordered.table)
    
}

get.IDs.identical.main.param <- function(ordered.table)
{
    values <- unique(ordered.table[,2])
    ids.list <- lapply(1:length(values), function(i)
    {
        t <- unlist(which(ordered.table[,2]==values[i]))
        t <- ordered.table[as.integer(t),1]
        return(as.integer(unlist(t)))
    })
    
    return(ids.list)
}

order.table.by.main.param <- function(ordered.table)
{
    values <- unique(ordered.table[,2])
    ids <- lapply(1:length(values), function(i)
    {
        t <- unlist(which(ordered.table[,2]==values[i]))
        return(as.integer(unlist(t)))
    })
    return( ordered.table[unlist(unlist(ids)), ] )
}

get.IDs.identical.sub.params <- function(ordered.table)
{
    values.mat <- as.matrix(t(unique(ordered.table[,-c(1,2)])))
    colnames(values.mat) <- colnames(ordered.table)[-c(1,2)]
    t.mat <- as.matrix(t(ordered.table[,-c(1,2)]))
    ids.list <- lapply(1:nrow(values.mat), function(n)
    {
        tmp <- sapply(1:nrow(t.mat), function(nt)
        {
            return(identical(t.mat[nt,],values.mat[n,]))
        })
        tmp <- unlist(tmp)
        tmp <- ordered.table[as.integer(tmp),1]
        return(as.integer(tmp))
    })
    return(ids.list)
}

order.table.by.sub.params <- function(ordered.table)
{
    values.mat <- as.matrix(t(unique(ordered.table[,-c(1,2)])))
    colnames(values.mat) <- colnames(ordered.table)[-c(1,2)]
    t.mat <- as.matrix(t(ordered.table[,-c(1,2)]))
    ids <- lapply(1:nrow(values.mat), function(n)
    {
        tmp <- sapply(1:nrow(t.mat), function(nt)
        {
            return(identical(t.mat[nt,],values.mat[n,]))
        })
        tmp <- unlist(tmp)
        return(as.integer(tmp))
    })
    return( ordered.table[unlist(unlist(ids)), ] )
}

compute.purity.points <- function(clusters, annotations, purity.matrix.clust)
{
    purity.values <- lapply(1:nrow(purity.matrix.clust), function(cl)
    {
        pur <- max(purity.matrix.clust[cl,])
        
        return(pur)
    })
    
    associated.annot <- lapply(1:length(clusters), function(cl)
    {
        annot.id <- 0
        lapply(1:length(annotations), function(an)
        {
            if(length(which(annotations[[an]][[2]]%in%clusters[[cl]][[2]]))==clusters[[cl]][[1]])
            {
                annot.id <<- an
            }
        })
        
        return(annot.id)
    })
    
    
    points.ids <- rep(NA,length(annotations))
    lapply(1:length(associated.annot), function(cl)
    {
        points.ids[[associated.annot[[cl]]]] <<- cl
    })
    
    return(list(points.ids, purity.values, associated.annot))
}