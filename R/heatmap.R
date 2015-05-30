#' Plot a matrix as a heatmap.
# source: http://www.phaget4.org/R/image_matrix.html
plot.matrix <- function(x, xlab='', ylab='', xLabels=NULL, yLabels=NULL, grid.on=TRUE, cex.axis=0.7, 
        grid.col="white", grid.lwd=1, ColorRamp=NULL, reverse.y=TRUE, zlim=NULL, na=0, na.text=as.character(na),
        pdf.file=NULL, pdf.width=NULL, pdf.height=NULL, region=c("all", "upper.tri","lower.tri"), axis.ivl=1, 
        plot.values=FALSE, value.col="black", hide.x.names=FALSE, hide.y.names=FALSE, col.dendro=NULL, 
        bottom.mar=8, left.mar=2, x.las=2, ...) {
    
    x.val <- x
    x[is.na(x)] <- na
    
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
        x <- x[reverse,,drop=F]
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
    par(mar=c(0, left.mar, 0, 2), oma=c(0, 0, ifelse(is.null(title), 1, 2.5), 0))
    if (!is.null(col.dendro)) {
        plot(col.dendro, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    par(mar=c(bottom.mar, left.mar, 0, 2))
    image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp,
        xlab=xlab, ylab=ylab, axes=FALSE, zlim=zlim)
    box()
    if (grid.on) {
        segments(1:(length(xLabels)+1)+0.5, 0, 1:(length(xLabels)+1)+0.5, length(yLabels)+1, col=grid.col, lwd=grid.lwd)
        segments(0, 1:(length(yLabels)+1)+0.5, length(xLabels)+3, 1:(length(yLabels)+1)+0.5, col=grid.col, lwd=grid.lwd)
    }
    if (plot.values) {
        print(ncol(x))
        for (i in 1:ncol(x)) {
            txt <- ifelse(is.na(x.val[,i]), na.text, x.val[,i])
            if (reverse.y) {
                txt <- rev(txt)
            }
            text(i, 1:nrow(x), txt, col=value.col)
        }
    }
    xax <- seq(1, length(xLabels), axis.ivl)
    yax <- seq(1, length(yLabels), axis.ivl)
    if (!hide.x.names) {
        axis(1, at=xax, labels=xLabels[xax], cex.axis=cex.axis, las=x.las)
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

#' Plot pairs of values as linear heatmaps.
#' 
#' @param A vector or two-column matrix of numeric values (see details)
#' @param B vector of numeric values the same length as `A`, or NULL if `A` is a matrix
#' @param xlim range of possible A and B values
#' @param palettes palettes for chosing colors for values less than and greater than `center. Either a vector or
#' a list of two vectors, where each vector contains two or more colors, ordered from closest to furthest from
#' `center`.
#' @param nslices number of color slices to use in each gradient
#' @param xlab x-axis label
#' @param spacing amount of space [0-0.5) between each heatmap
#' @param bottom.to.top if TRUE, plot values from the bottom to the top of the plot, otherwise top-down
#' @param ... additional arguments to be passed to `plot`
#' 
#' @examples
#' A <- matrix(runif(20) / 2, nrow=10, ncol=2, dimnames=list(paste('Sample', 1:10), c("A","B")))
#' par(mar=c(4,5,1,1))
#' paired.linear.heatmap(A, xlab='Allelic Ratio (B/A)')
paired.linear.heatmap <- function(A, B=NULL, xlim=c(0, 0.5), palettes=brewer.pal(9, "Reds"), nslices=50, xlab="", 
                                  spacing=0.05, bottom.to.top=TRUE, ...) {
    if (is.null(B)) {
        stopifnot(is.matrix(A) || is.data.frame(A))
        B <- A[,2]
        A <- A[,1]
    }
    else {
        stopifnot(length(A) == length(B))
    }
    if (is.null(xlim)) {
        xlim <- range(min(c(A,B)), max(c(A,B)))
    }
    else {
        stopifnot(all(A >= xlim[1] & A <= xlim[2]))
        stopifnot(all(B >= xlim[1] & B <= xlim[2]))
    }
    
    if (!bottom.to.top) {
        A <- rev(A)
        B <- rev(B)
    }
    
    # shift A and B to range 0-1
    rng <- xlim[2] - xlim[1]
    A.norm <- (A - xlim[1]) / rng
    B.norm <- (B - xlim[1]) / rng
    
    if (is.list(palettes) && length(palettes) >= 2) {
        A.pal <- palettes[[1]]
        B.pal <- palettes[[2]]
    }
    else {
        A.pal <- B.pal <- unlist(palettes)
    }
    A.colors <- colorRampPalette(A.pal)(nslices)
    A.colidxs <- round(nslices * A.norm)
    
    B.colors <- colorRampPalette(B.pal)(nslices)
    B.colidxs <- round(nslices * B.norm)
    
    N <- length(A)
    plot(0:N, type="n", xaxt="n", yaxt="n", xlab=xlab, ylab="", xlim=c(-1, 1), ...)
    for (i in 1:N) {
        gradient.rect(-A[i], i-1+spacing, 0, i-spacing, col=rev(A.colors[1:A.colidxs[i]]), border='black')
        gradient.rect(0, i-1+spacing, B[i], i-spacing, col=B.colors[1:B.colidxs[i]], border='black')
    }
    xlabs <- seq(xlim[1], xlim[2], length.out=6)
    axis(1, seq(-1,1,0.2), c(rev(xlabs), xlabs[-1]))
    axis(2, c(1:N)-0.5, names(A), las=2)
}

#' Plot values as linear heatmaps.
#'
#' @param x vector or two-column matrix of numeric values (see details)
#' @param center central value; each value of x will be shown relative to its difference from `center`
#' @param min minimum value; must be <= `center`
#' @param max value; must be >= `center`
#' @param palettes palettes for chosing colors for values less than and greater than `center. Either a vector or
#' a list of two vectors, where each vector contains two or more colors, ordered from closest to furthest from
#' `center`.
#' @param nslices number of color slices to use in each gradient
#' @param xlab x-axis label
#' @param spacing amount of space [0-0.5) between each heatmap
#' @param bottom.to.top if TRUE, plot values from the bottom to the top of the plot, otherwise top-down
#' @param ... additional arguments to be passed to `plot`
#' 
#' @examples
#' vals <- runif(10)
#' names(vals) <- paste('Sample', 1:10)
#' par(mar=c(4,5,1,1))
#' linear.heatmap(vals, xlab='Allelic Ratio (B/A)')
linear.heatmap <- function(x, center=0.5, minval=0.0, maxval=1.0, palettes=brewer.pal(9, "Reds"), nslices=50, xlab="", 
                           spacing=0.05, bottom.to.top=TRUE, ...) {
    if (is.null(minval)) {
        minval <- min(x)
    }
    else {
        stopifnot(all(x >= minval))
    }
    if (is.null(maxval)) {
        maxval <- max(x)
    }
    else {
        stopifnot(all(x <= maxval))
    }
    stopifnot(minval <= center)
    stopifnot(maxval >= center)
    
    if (!bottom.to.top) {
        x <- rev(x)
    }
    
    if (is.list(palettes) && length(palettes) >= 2) {
        l.pal <- palettes[[1]]
        r.pal <- palettes[[2]]
    }
    else {
        l.pal <- r.pal <- unlist(palettes)
    }
    l.colors <- colorRampPalette(l.pal)(nslices)
    r.colors <- colorRampPalette(r.pal)(nslices)
    
    l.idxs <- which(x < center)
    l.colidxs <- round(nslices * ((center - x[l.idxs]) / (center - minval)))
    
    r.idxs <- which(x > center)
    r.colidxs <- round(nslices * ((x[r.idxs] - center) / (maxval - center)))
    
    N <- length(x)
    plot(0:N, type="n", yaxt="n", ylab="", xlab=xlab, xlim=c(minval, maxval), ...)
    for (ii in 1:length(l.idxs)) {
        i <- l.idxs[ii]
        gradient.rect(x[i], i-1+spacing, center, i-spacing, col=rev(l.colors[1:l.colidxs[ii]]), border='black')
    }
    for (ii in 1:length(r.idxs)) {
        i <- r.idxs[ii]
        gradient.rect(center, i-1+spacing, x[i], i-spacing, col=r.colors[1:r.colidxs[ii]], border='black')
    }
    axis(2, c(1:N)-0.5, names(x), las=2)
}