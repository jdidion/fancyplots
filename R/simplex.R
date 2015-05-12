DEFAULT_COLRAMP <- colorRampPalette(blues9[-(1:3)])
DEFAULT_PCH <- 21

#' A simplex is a way to plot a specific type of multi-variate information in two-dimensions.
#' Each data point is taken to be a list of proportions that sum to 1.0. An example of data
#' that can be plotted with a simplex is admixture assignments from three or more founder
#' populations: the higher the value of a specific variable (population), the more of the 
#' sample's admixture is derived from that population.
#'
#' Points are plotted in the interior of a shape with a number of vertices equal to the number 
#' of variables (i.e. a triangle for data points in three dimensions, a square for four dimensions,
#' etc). The higher the value of a variable, the closer the point is to the vertex representing
#' that variable. The color of the point is also determined by it's closeness to each vertex.
#'
#' If `series` is a list, then each element is a Nx3 matrix, where each row is translated 
#' into a data point. The sum of each row in the data matrix must be 1. If `series` is a 
#' list of lists, each element has a `data` member, which is the matrix, and must also 
#' have pch and colramp elements. If `show.legend` is TRUE, each list must also have a name
#' element. Other optional elements are: highlight.idx and highlight.col, where highlight.idx 
#' can be an index or list of indices to highight with the specified color (this works best if 
#' the colramp is grayscale).
#'
#' @param series either a list or a list of lists (see details).
#' @param main plot title.
#' @param corner.labels labels for each of the vertices.
#' @param cex point cex.
#' @param show.ticks show tick marks on the sides of the triangle.
#' @param show.legend whether to show a legend defining the colors used for each variable.
#' @param legend.pos where to put the legend.
#' @param plotwidth size of plot (in inches).
#' @param new.plot logical whether to create a new plot or overplot on the existing device.
simplex.triangle <- function(series, main=NULL, corner.labels=NULL, cex=NULL,
                             show.ticks=FALSE, show.legend=FALSE,
                             legend.pos='topleft', plotwidth=7, new.plot=TRUE) {

    nseries <- length(series)
    
    # draw the frame
    if (show.legend) {
        if (new.plot) {
            dev.new(width=plotwidth+nseries, height=plotwidth)
        }
        par(oma=c(0,0,0,0),mar=c(0.1,0.1,ifelse(is.null(main), 0, 1),0.1))
        layout(matrix(seq(1, nseries), nrow=1), 
               widths=c(plotwidth, rep(1, nseries)), 
               heights=rep(plotwidth, nseries+1))
    }
    
    plot(0:1, type='n', axes=FALSE, xlim=c(-1,1), ylim=c(-1,1), xlab='', ylab='')
    if (!is.null(main)) {
        title(main, font.main=2)
    }
    # draw the triangle
    draw.polygon(show.ticks=show.ticks, corner.labels=corner.labels)

    if (!is.list(series[[1]])) {
        l = list()
        l[[series$name]] <- series
        series <- l
    }

    nm <- names(series)
    legs <- list()
    for (n in nm) {
        legs[[n]] <- simplex.plot.series(series[[n]], cex=cex)
    }
    
    if (show.legend) {
        pch <- unlist(lapply(series, function(x) default(x, 'pch', DEFAULT_PCH)))
        col <- unlist(lapply(series, function(x) default(x, 'colramp', DEFAULT_COLRAMP)(1)))
        bg <- unlist(lapply(series, function(x) ifelse(default(x, 'filled', TRUE), col, NULL)))
        legend(legend.pos, bty="n", legend=nm, pch=pch, col=col, pt.bg=bg)

        for (n in names(legs)) {
            plot(NA, xlim=c(0,10), ylim=c(0,11), type="n", ann=FALSE, axes=FALSE)
            rect(0, 1:10, 2, 2:11, border=NA, col=legs[[n]]$col)
            text(4, (1:10)+0.5, signif(legs[[n]]$density, 2))
         }
    }
}

default <- function(list, key, value) {
    if (key %in% names(list)) {
        return(list[[key]])
    }
    else {
        return(value)
    }
}

as.cartesian <- function(data, sidelen=1) {
    sides <- ncol(data)
    alpha <- 2 * (pi / sides)
    r <- sidelen / (2 * sin(alpha / 2))
    dx <- matrix(0, nrow(data), sides)
    dy <- matrix(0, nrow(data), sides)
    
    for (i in 1:sides) {
        a <- alpha * (i-1)
        dpt <- sapply(data[,i], function(x) pol2cart(r * x, a))
        dx[,i] <- dpt[1,]
        dy[,i] <- dpt[2,]
    }
    
    return(cbind(apply(dx, 1, sum), apply(dy, 1, sum)))
}

# returns c(name, pch, col)
simplex.plot.series <- function(series, sidelen=1, cex=NULL) {
    name <- series$name
    data <- series$data
    pch <- default(series, 'pch', DEFAULT_PCH)
    colramp <- default(series, 'colramp', DEFAULT_COLRAMP)
    filled <- default(series, 'filled', TRUE)
    hidx <- default(series, 'highlight.idx', NULL)
    hcol <- default(series, 'highlihgt.col', NULL)
    
    xy <- as.cartesian(data, sidelen)
    
    dcols <- densCols(xy, colramp=colramp)
    
    dens <- grDevices:::.smoothScatterCalcDensity(xy, nbin=128)
    dens <- as.numeric(dens$fhat)
    dens <- dens[dens>0]
    dens.rng <- range(dens)
    colLegend <- data.frame(density=seq(min(dens), max(dens), len=10), color=I(colramp(10)))
    
    bg <- NULL
    if (filled) {
        bg <- dcols
    }
    points(xy[,1], xy[,2], cex=cex, pch=pch, col=dcols, bg=bg)

    if (!is.null(hidx)) {
        points(xy[hidx,1], xy[hidx,2], cex=cex, pch=pch, col=hcol)
    }
    
    return(colLegend)
}

# polygon is drawn with center at (0,0).
# cornerlabels start at top and go clockwise.
draw.polygon <- function(sides=3, sidelen=1, corner.labels=NULL, show.ticks=FALSE, tick.size=0.02) {
    # Draw the polygon
    qpi <- (pi / 4)
    alpha <- 2 * (pi / sides)
    r <- sidelen / (2 * sin(alpha / 2))
    corners <- sapply(0:sides, function(x) pol2cart(r, alpha * x))
    lines(corners[1,], corners[2,])

    # label corners
    if (!is.null(corner.labels)) {
        par(xpd=TRUE)
        for (i in 1:sides) {
            pos <- octant2pos(alpha * (i - 1) / qpi)
            text(corners[1,i], corners[2,i], corner.labels[i], pos=pos)
        }
    }

    # draw ticks
    if (show.ticks) {
        t.r <- r - tick.size
        t.alpha <- alpha / 10
        t.int <- sidelen / 10
        t.a <- c((pi-alpha) / 2)
        t.c <- c(0)
        t.r <- c(r)
        for (i in 1:9) {
            s <- sqrt(t.r[i]^2 + t.int^2 - (2 * t.r[i] * t.int * cos(t.a[i])))
            t.r <- append(t.r, s)
            c <- asin(t.int * sin(t.a[i]) / s)
            t.a <- append(t.a, t.a[i] + c)
            t.c <- append(t.c, c + t.c[i])
        }

        pos <- calc.pos(sides, alpha) 
        pct <- as.character(c(seq(10, 40, by=10), 50, seq(40, 10, by=-10)))
        for (i in 1:sides) {
            a <- alpha * (i-1)
            m <- -1 * (corners[1,i+1] - corners[1,i]) / (corners[2,i+1] - corners[2,i])
            dx <- sqrt(tick.size^2 / (m^2 + 1))
            dy <- m * dx
            mid <- pol2cart(t.r[6], a + t.c[6])
            dir <- ifelse(dist(c(0,0), mid) > dist(c(0,0), c(mid[1] + dx, mid[2] + dy)), 1, -1)
            
            for (j in 2:10) {
                xy1 <- pol2cart(t.r[j], a + t.c[j])
                x2 <- ifelse(dir > 0, xy1[1] + dx, xy1[1] - dx)
                y2 <- ifelse(dir > 0, xy1[2] + dy, xy1[2] - dy)
                lines(c(xy1[1], x2), c(xy1[2], y2))
                text(xy1[1], xy1[2], pct[j-1], pos=pos[i], cex=0.5)
            }
        }
    }
}

pol2cart <- function(r, theta) {
  c(r*sin(theta), r*cos(theta))
}

rad2deg <- function(rad) {
    return(rad * 180 / pi)
}

dist <- function(p1, p2) {
    return(sqrt((p2[1] - p1[1])^2 + (p2[2]-p1[2])^2))
}

octant2pos <- function(octant) {
    if (octant > 3 && octant <= 5) {
        return(1)
    }
    else if (octant > 5 && octant <= 7) {
        return(2)
    }
    else if (octant > 7 || octant <= 1) {
        return(3)
    }
    else {
        return(4)
    }
}

calc.pos <- function(sides, alpha) {
    h <- sides / 2
    mx <- ceiling(h)
    odd <- mx != h
    r <- c()
    l <- c()
    for (i in 1:mx) {
        if (odd && i == mx) {
            r <- append(r, 1)
        }
        else if (i * alpha <= (pi * 72 / 180)) {
            r <- append(r, 3)
            l <- append(r, 3, 0)
        }
        else if ((i-1) * alpha >= (pi * 120 / 180)) {
            r <- append(r, 1)
            l <- append(l, 1, 0)
        }
        else {
            r <- append(r, 4)
            l <- append(l, 2, 0)
        }
    }
    return(c(r, l))
}

norm.rows <- function(m) {
    return(t(apply(m, 1, function(x) x / sum(x))))
}
