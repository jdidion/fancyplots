plot.pca.matrix <- function(m, group.df=NULL, PCs=c(1,2), returnData=FALSE, plotData=TRUE, returnFig=FALSE, ...) {
    pca <- prcomp(t(as.matrix(m)), center=TRUE, scale.=FALSE)
    pct.var <- percent.var(pca)
    if (plotData) {
        fig <- plot.pca(pca, group.df, PCs, pct.var, !returnFig, ...)
    }
    if (returnData || returnFig) {
        result <- list()
        if (returnData) {
            result$data <- pca
            result$pct.var <- pct.var
        }
        if (returnFig) {
            result$fig <- fig
        }
        invisible(result)
    }
}

percent.var <- function(pca) pca$sdev^2/sum(pca$sdev^2)

plot.pca <- function(pca, group.df, PCs=c(1,2), pct.var=NULL, plot=TRUE, ...) {
    if (is.null(pct.var)) {
        pct.var <- percent.var(pca)
    }
    x <- pca$x[,PCs[1]]
    y <- pca$x[,PCs[2]]
    plot.eig(x, y, group.df, PCs, pct.var=pct.var, plot=plot, dim.label="PC", ...)
}

#' Perform MDS analysis and plot results.
#' 
#' @param d distance matrix
#' @param group.df data.frame of grouping variable(s)
plot.mds <- function(d, group.df, dim.plot=c(1, 2), ndim=max(dim.plot), returnData=FALSE, 
                     plotData=TRUE, returnFig=FALSE, main=NULL) {
    mds <- suppressWarnings(cmdscale(as.dist(d), k=ndim, eig=TRUE))
    coords <- mds$points
    gof.str <- paste0("GoF: ", round(mds$GOF[1], 3))
    if (is.null(main)) {
        main <- gof.str
    }
    else {
        main <- paste(main, gof.str)
    }
    if (dim.plot[1] > ncol(coords)) {
        x <- rep.int(0, nsamples)
        warning(paste("dimension", dim.plot[1], "is degenerate or all zero"))
    }
    else x <- coords[, dim.plot[1]]
    if (dim.plot[2] > ncol(coords)) {
        y <- rep.int(0, nsamples)
        warning(paste("dimension", dim.plot[2], "is degenerate or all zero"))
    }
    else y <- coords[, dim.plot[2]]
    if (plotData) {
        fig <- plot.eig(x, y, group.df, dim.plot, plot=!returnFig, main=main, dim.label="Coordinate")
    }
    if (returnData || returnFig) {
        result <- list()
        if (returnData) {
            result$data <- mds
        }
        if (returnFig) {
            result$fig <- fig
        }
        invisible(result)
    }
}

plot.eig <- function(x, y, group.df, PCs, pt.size=2, pct.var=NULL, plot=TRUE, main=NULL, dim.label="PC") {
    d <- data.frame(PC1=x, PC2=y, group=1)
    if (!is.null(group.df)) {
        d$group <- factor(apply(group.df, 1, paste, collapse = " : "))
        #d <- cbind(d, group.df, names=rownames(group.df))
    }
    x.lab <- paste(dim.label, PCs[1])
    y.lab <- paste(dim.label, PCs[2])
    if (!is.null(pct.var)) {
        x.lab <- paste0(x.lab, ": ", round(pct.var[PCs[1]] * 100), "% variance")
        y.lab <- paste0(y.lab, ": ", round(pct.var[PCs[2]] * 100), "% variance")
    }
    g <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
        geom_point(size=pt.size) + xlab(x.lab) + ylab(y.lab)
    if (!is.null(main)) {
        g <- g + ggtitle(main)
    }
    if (plot) {
        print(g)
    }
    invisible(g)
}