draw.cumulated.filled.plots <- function(points.list, max.height=2, x.values.range=c(0.15,1), x.lab="par", y.lab="val", nmb.pts)
{
    plot(x.values.range,c(0,max.height),xlab = x.lab, ylab = y.lab)
    y.vals <- as.numeric(unlist(points.list[[1]]))
    x.pts <- c(1:nmb.pts,rev(1:nmb.pts))
    lapply(1:(length(points.list)-1), function(i)
    {
        y.pts <- c(y.vals,rev(y.vals+as.numeric(unlist(points.list[[i+1]]))))
        polygon(x.pts, y.pts,col=colors()[i*8])
        
        y.vals <<- unlist(y.vals) + unlist(as.numeric(points.list[[i+1]]))
    })
}


draw.F.score.barplot <- function(F.score.matrix, populations.names, populations.sizes)
{
    scores.list.col <- sapply(1:nrow(F.score.matrix), function(r)
    {
        return(which(F.score.matrix[r,]==max(F.score.matrix[r,]))[[1]])
    })
    pop.order <- order(populations.sizes)
    
    scores.list <- rep(NA,length(scores.list.col))
    lapply(1:length(scores.list), function(i)
    {
        col.id <- as.integer(unlist(scores.list.col)[[i]])
        scores.list[[col.id]] <<- F.score.matrix[i,col.id]
    })
    
    names(scores.list) <- as.character(unlist(populations.names))
    scores.list <- scores.list[pop.order]
    
    barplot(scores.list, main="F scores by population (ordered by %events)", horiz = T, 
            names.arg = names(scores.list), cex.names = 0.8, xlim = c(0,1.05))
    
    
}


plot.selected.clusters <- function(val.mat, clusters, markers)
{
    highlighted = rep("gray7",nrow(val.mat))
    lapply(clusters, function(cl)
    {
        highlighted[unlist(as.integer(cl))] <<- "firebrick"
    })
    
    plot(val.mat[,markers], col=highlighted, xlim=c(-0.5,4.5), ylim=c(-0.5,4.5), pch=".")
    
    lapply(1:length(clusters), function(cl.id)
    {
        cl <- clusters[[cl.id]]
        xco <- mean(val.mat[unlist(as.integer(cl)), markers[1]])
        yco <- mean(val.mat[unlist(as.integer(cl)), markers[2]])
        text(xco,yco, names(clusters)[cl.id], col = "blue")
    })
}

plot.purity.by.annot <- function(annot.clusters, purity.val, annot.sizes, purity.threshold)
{
    annot.order <- order(annot.sizes)
    
    pts <- c()
    lapply(1:length(annot.clusters), function(an)
    {
        an.id <- annot.order[an]
        lapply(1:length(annot.clusters[[an.id]]), function(cl)
        {
            pts <<- c(pts, as.numeric(purity.val[as.integer(annot.clusters[[an.id]][[cl]])]))
            names(pts)[length(pts)] <<- as.integer(annot.clusters[[an.id]][[cl]])
        })
    })
    
    Y <- pts
    X <- as.integer(names(pts))
    below.points <- which(Y>=purity.threshold)
    Y <- Y[below.points]
    X <- X[below.points]
    
    plt.colors <- sample(colours(), length(annot.clusters), replace = F)
    plot(X, Y, ylim=c(0,1.05), xlim=c(0,max(X)+1))
    x0 <- 0
    y0 <- 0
    lapply(1:length(annot.clusters), function(an)
    {
        an.id <- annot.order[an]
        x1 <- x0+length(annot.clusters[[an.id]])
        y1 <- 1.05
        rect(x0,y0,x1,y1, col = plt.colors[an])
        x0 <<- x1
    })
    points(as.integer(names(pts)), pts, pch="x")
    abline(h=purity.threshold)
    
    return(X)
}