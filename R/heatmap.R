#' Plot a matrix as a heatmap.
# source: http://www.phaget4.org/R/image_matrix.html
plot.matrix <- function(x, xlab='', ylab='', xLabels=NULL, yLabels=NULL, grid.on=TRUE, cex.axis=0.7, 
        grid.col="white", grid.lwd=1, ColorRamp=NULL, reverse.y=TRUE, zlim=NULL, pdf.file=NULL, 
        pdf.width=NULL, pdf.height=NULL, region=c("all", "upper.tri","lower.tri"), axis.ivl=1, 
        hide.x.names=FALSE, hide.y.names=FALSE, col.dendro=NULL, ...) {
    
    if (is.null(zlim)) {
        zlim <- c(min(x), max(x))
    }
    
    if (is.null(xLabels)) {
        xLabels <- colnames(x)
    }
    if (is.null(yLabels)) {
        yLabels <- rownames(x)
    }
    
    title <- c()

    region <- match.arg(region, c("all", "upper.tri","lower.tri"))
    if (region == "upper.tri") {
        x[lower.tri(x)] <- zlim[1]-1
    }
    else if (region == "lower.tri") {
        x[upper.tri(x)] <- zlim[1]-1
    }

    # check for additional function arguments
    if (length(list(...))) {
        Lst <- list(...)
        if (!is.null(Lst$zlim)) {
            zlim <- Lst$zlim[1]
        }
        if (!is.null(Lst$yLabels)) {
            yLabels <- c(Lst$yLabels)
        }
        if (!is.null(Lst$xLabels)) {
            xLabels <- c(Lst$xLabels)
        }
        if (!is.null(Lst$title)) {
            title <- Lst$title
        }
    }
    
    if (is.null(xLabels)) {
        xLabels <- c(1:ncol(x))
    }
    
    if (is.null(yLabels)) {
        yLabels <- c(1:nrow(x))
    }

    # Red and green range from 0 to 1 while Blue ranges from 1 to 0
    if (is.null(ColorRamp)) {
        ColorRamp <- rgb(seq(0, 1, length=256), seq(0, 1, length=256), seq(1, 0, length=256))
    }
    ColorLevels <- seq(zlim[1], zlim[2], length=length(ColorRamp))

    # Reverse Y axis
    if (reverse.y) {
        reverse <- nrow(x):1
        yLabels <- yLabels[reverse]
        x <- x[reverse, ]
    }

    if (!is.null(pdf.file)) {
        pdf(pdf.file, pdf.width, pdf.height)
    }

    if (is.null(col.dendro)) {
        layout(matrix(data=c(1, 2), nrow=1, ncol=2), widths=c(6, 1), heights=1)
    }
    else {
        layout(matrix(data=c(1, 2, 0, 3), nrow=2, ncol=2), widths=c(6, 1), heights=c(1, 6))
    }
    par(mar=c(0, 2, 0, 2), oma=c(0, 0, ifelse(is.null(title),0,2.5), 0))
    if (!is.null(col.dendro)) {
        plot(col.dendro, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    par(mar=c(8, 2, 0, 2))
    image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp,
        xlab=xlab, ylab=ylab, axes=FALSE, zlim=zlim)
    box()
    if (grid.on) {
        segments(1:(length(xLabels)+1)+0.5, 0, 1:(length(xLabels)+1)+0.5, length(yLabels)+1, col=grid.col, lwd=grid.lwd)
        segments(0, 1:(length(yLabels)+1)+0.5, length(xLabels)+3, 1:(length(yLabels)+1)+0.5, col=grid.col, lwd=grid.lwd)
    }
    xax <- seq(1, length(xLabels), axis.ivl)
    yax <- seq(1, length(yLabels), axis.ivl)
    if (!hide.x.names) {
        axis(1, at=xax, labels=xLabels[xax], cex.axis=cex.axis, las=2)
    }
    if (!hide.y.names) {
        axis(2, at=yax, labels=yLabels[yax], cex.axis=cex.axis, las=2)
    }
    par(mar=c(3, 2.5, 2.5, 2))
    image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
        col=ColorRamp, xlab="", ylab="", xaxt="n")
    if (!is.null(title)) {
        title(main=title, outer=TRUE)
    }
    layout(1)

    if (!is.null(pdf.file)) {
        dev.off()
    }
}
