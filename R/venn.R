#' Draw a Venn diagram of intersections between lists of items. Supports between 2-5 lists.
#
#' @param lists matrix with between 2-5 columns or list with 2-5 vectors.
#' @param plotData logical whether to draw the plot
#' @param returnFig logical whether to return the figure object
#' @param pal vector of colors to use for filling circles
#' 
#' TODO: allow for different colors in each region
plot.venn <- function(lists, plotData=TRUE, returnFig=FALSE, pal=c("lightgreen","orange","violet","yellow","brown"), ...) {
    if (is.data.frame(lists) || is.matrix(lists)) {
        N <- ncol(lists)
        stopifnot(N >= 2 && N <= 5)
        l <- list()
        for (i in 1:N) {
            l[[colnames(lists)[i]]] <- lists[,i]
        }
        lists <- l
    }
    else if (is.list(lists)) {
        N <- length(lists)
        stopifnot(N >= 1 && N <= 5)
    }
    else {
        stop("Invalid lists object")
    }
    fig <- venn.diagram(lists, NULL, na="stop", fill=pal[1:length(lists)], ...)
    if (plotData) {
        grid.draw(fig)
    }
    retval <- NULL
    if (returnFig) {
        retval <- fig
    }
    invisible(retval)
}
