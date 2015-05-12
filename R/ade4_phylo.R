#' Draw a radial tree for an ade4 phylog object.
radial.phylog <- function (phylog, circle=1, show.leaves=TRUE, show.leaf.nodes=TRUE,
        labels.leaves=names(phylog$leaves), cex.leaves=1, bg.leaves='black', pch.leaves=21,
        fg.text.leaves='black', ragged.leaves=FALSE, pad.leaves=0, cex.leaf.points=1,
        number.leaves=FALSE, cex.leaf.numbers=0.5,
        labelptr.show=TRUE, labelptr.col=grey(0.7), labelptr.pad=NULL,
        show.nodes=FALSE, labels.nodes=phylog$nodes, cex.nodes=1, bg.nodes='white', pch.nodes=21,
        bg.colramp.nodes=NULL, fg.nodes='black', adj.nodes=0,
        lwd.branch=1, lwd.branch.default=1, col.branch='black', col.branch.default='black',
        col.branch.colramp=NULL, weights.branch=NULL,
        col.center='red', cex.center=2,
        legend.nodes=FALSE, legend.title.nodes="Nodes",
        legend.branch=FALSE, legend.title.branch="Branches",
        main=NULL, pdf.file=NULL, pdf.width=7, pdf.height=7) {
    if (is.character(phylog)) {
        phylog <- newick2phylog(readLines(phylog))
    }
    if (!inherits(phylog, "phylog"))
        stop("Non convenient data")
    if (circle < 0)
        stop("'circle': non convenient value")

    retval <- list(phylog=phylog)

    if (!is.null(pdf.file)) {
        pdf(pdf.file, pdf.width, pdf.height)
    }

    leaves.number <- length(phylog$leaves)
    leaves.names <- names(phylog$leaves)
    nodes.number <- length(phylog$nodes)
    nodes.names <- names(phylog$nodes)
    if (length(labels.leaves) == 1) {
        name.map <- read.table(labels.leaves, sep='=', header=FALSE, colClasses='character', row.names=1)
        labels.leaves <- name.map[substr(leaves.names, 2, nchar(leaves.names)), 1]
    }
    if (length(labels.leaves) != leaves.number)
        labels.leaves <- leaves.names
    if (length(labels.nodes) != nodes.number)
        labels.nodes <- names(phylog$nodes)

    dis <- phylog$droot
    dis <- (dis / max(dis)) * circle
    dist.leaves <- dis[leaves.names]
    dist.nodes <- dis[nodes.names]
    if (nodes.number == 1) {
        d.rayon <- circle
    }
    else {
        d.rayon <- circle / (nodes.number - 1)
    }
    theta <- (2 * pi) / leaves.number
    alpha <- theta * (1:leaves.number)
    names(alpha) <- leaves.names
    x <- dist.leaves * cos(alpha)
    y <- dist.leaves * sin(alpha)
    ang <- rep(0, length(dist.nodes))
    names(ang) <- names(dist.nodes)
    ang <- c(alpha, ang)

    opar <- par(mar=c(0.1, 0.1, 0.1, 0.1), oma=c(1 + ifelse(legend.nodes || legend.branch, 1, 0), 1, 1, 1))
    on.exit(par(opar))
    plot.default(0, 0, type="n", asp=1, xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(-2, 2),
        ylim=c(-2, 2), xaxs="i", yaxs="i", frame.plot=FALSE)
    if (!is.null(main)) {
        title(main, line=-1)
    }

    if (!is.null(bg.colramp.nodes)) {
        if (length(bg.colramp.nodes) == 1) {
            bg.colramp.nodes <- brewer.pal(9, bg.colramp.nodes)
        }
        nodelabs <- as.numeric(labels.nodes)
        node.ramp.start <- as.integer(floor(min(nodelabs)))
        node.ramp.end <- as.integer(ceiling(max(nodelabs)))
        node.ramp <- make.ramp(bg.colramp.nodes, node.ramp.start, node.ramp.end)
        bg.nodes <- sapply(nodelabs, node.ramp)
    }

    if (!is.null(col.branch.colramp)) {
        if (length(col.branch.colramp) == 1) {
            col.branch.colramp <- brewer.pal(9, col.branch.colramp)
        }
        branch.ramp.start <- as.integer(floor(min(unlist(weights.branch))))
        branch.ramp.end <- as.integer(ceiling(max(unlist(weights.branch))))
        branch.ramp <- make.ramp(col.branch.colramp, branch.ramp.start, branch.ramp.end)
        col.branch <- lapply(names(phylog$paths), function(x) {
            val <- weights.branch[[x]]
            if (is.null(val)) {
                val <- branch.ramp.start
            }
            branch.ramp(val)
        })
        names(col.branch) <- names(phylog$paths)
    }

    # Draw branches
    for (i in 1:length(phylog$parts)) {
        w <- phylog$parts[[i]]
        but <- names(phylog$parts)[i]
        ang[but] <- mean(ang[w])

        if (length(col.branch) == 1) {
            col <- c(col.branch, col.branch)
        }
        else {
            col <- unlist(col.branch[w])
        }
        col[is.na(col)] <- col.branch.default

        if (length(lwd.branch) == 1) {
            lwd = c(lwd.branch, lwd.branch)
        }
        else {
            lwd = lwd.branch[w]
        }
        col[is.na(lwd)] <- lwd.branch.default

        # perpendicular to radius
        b <- range(ang[w])
        m <- ((b[2] - b[1]) / 2)
        a1.seq <- c(seq(b[1], b[2]-m, by = pi/180), b[2]-m)
        a2.seq <- c(seq(b[2]-m, b[2], by = pi/180), b[2])
        lines(dis[but] * cos(a1.seq), dis[but] * sin(a1.seq), lwd=lwd[1], col=col[1])
        lines(dis[but] * cos(a2.seq), dis[but] * sin(a2.seq), lwd=lwd[2], col=col[2])

        # parallel to radius
        x1 <- dis[w] * cos(ang[w])
        y1 <- dis[w] * sin(ang[w])
        x2 <- dis[but] * cos(ang[w])
        y2 <- dis[but] * sin(ang[w])
        segments(x1, y1, x2, y2, col=col, lwd=lwd)
    }

    if (is.null(labelptr.pad)) {
        labelptr.pad <- d.rayon
    }

    r.ptr <- ifelse(ragged.leaves, dist.leaves, circle) + labelptr.pad
    xptr <- r.ptr * cos(alpha)
    yptr <- r.ptr * sin(alpha)

    r.leaf <- r.ptr + pad.leaves
    xleaf <- r.leaf * cos(alpha)
    yleaf <- r.leaf * sin(alpha)

    # Draw leaves and leaf labels
    if (show.leaves) {
        retval$leaves <- list()
        for (i in 1:leaves.number) {
            drawn <- TRUE
            if (is.list(labels.leaves[i])) {
                if (is.integer(labels.leaves[[i]][[1]])) {
                    # This is wrong, but I'm not sure how to get the true point of a plotting symbol
                    pt.size <- max(convert.size(par("cin"))) * cex.leaf.points / 2.75
                    cols <- unlist(sapply(1:length(labels.leaves[[i]]), function(l)
                        rep(names(labels.leaves[[i]])[[l]], labels.leaves[[i]][[l]])))
                    npts <- length(cols)
                    pt.x <- NULL
                    pt.y <- NULL
                    row <- 1
                    while (npts > 0) {
                        row.r <- r.leaf + (pt.size * row)
                        if (length(row.r) > 1) {
                            row.r <- row.r[i]
                        }
                        max.pts <- floor((row.r * theta) / pt.size)
                        if (max.pts %% 2 == 0) {
                            max.pts <- max.pts - 1
                        }
                        t <- (theta * 0.9) / max.pts
                        n <- max(min(max.pts, npts), 1)
                        m <- .alternating(n)
                        pt.x <- c(pt.x, row.r * cos(alpha[i] + (t * m)))
                        pt.y <- c(pt.y, row.r * sin(alpha[i] + (t * m)))
                        npts <- npts - n
                        row <- row + 1
                    }
                    points(pt.x, pt.y, col=cols, pch=20, cex=cex.leaf.points)
                }
                else {
                    par(srt=alpha[i] * 360/2/pi)
                    for (col in names(labels.leaves[[i]])) {
                        lab <- labels.leaves[[i]][[col]]
                        text(xleaf[i], yleaf[i], lab, adj=0, col=col, cex=par("cex") * cex.leaves)
                    }
                }
            }
            else if (cex.leaves > 0 && nchar(labels.leaves[i]) > 0) {
                par(srt=alpha[i] * 360/2/pi)
                col.txt <- ifelse(length(fg.text.leaves) == 1, fg.text.leaves, fg.text.leaves[i])
                text(xleaf[i], yleaf[i], labels.leaves[i], adj=0, col=col.txt, cex=par("cex") * cex.leaves)
                retval$leaves[[i]] <- c(xleaf[i], yleaf[i])
            }
            else {
                drawn <- FALSE
            }
            if (drawn && labelptr.show) {
                segments(xptr[i], yptr[i], x[i], y[i], col=labelptr.col)
            }
            if (show.leaf.nodes) {
                pch <- ifelse(length(pch.leaves)==1, pch.leaves, pch.leaves[i])
                col <- ifelse(length(bg.leaves)==1, bg.leaves, bg.leaves[i])
                points(x[i], y[i], pch=pch, bg=col, cex=par("cex") * show.leaves)
            }
            if (number.leaves) {
                par(srt=alpha[i] * 360/2/pi)
                text(xptr[i], yptr[i], as.character(i), cex=cex.leaf.numbers)
            }
        }
    }

    # Draw center node
    points(0, 0, pch=21, cex=par("cex") * cex.center, bg=col.center)

    # Draw internal nodes and labels
    if (show.nodes) {
        delta <- strwidth(as.character(length(dist.nodes)), cex=par("cex") * ifelse(cex.nodes == 0, 0.3, cex.nodes))
        for (j in 1:(length(dist.nodes)-1)) {
            i <- names(dist.nodes)[j]
            #x1 <- (dis[i] - (adj.nodes * circle)) * cos(ang[i])
            #y1 <- (dis[i] - (adj.nodes * circle)) * sin(ang[i])
            x1 <- dis[i] * cos(ang[i])
            y1 <- dis[i] * sin(ang[i])
            bg <- ifelse(length(bg.nodes) == 1, bg.nodes, bg.nodes[j])
            pch <- ifelse(length(pch.nodes)==1, pch.nodes, pch.nodes[i])
            points(x1, y1, pch=pch, bg=bg, cex=par("cex") * show.nodes)
            if (cex.nodes > 0) {
                par(srt = (ang[i] * 360/2/pi + 90))
                fg <- ifelse(length(fg.nodes) == 1, fg.nodes, fg.nodes[j])
                text(x1, y1, labels.nodes[j], adj=0.5, cex=par("cex") * cex.nodes, col=fg)
            }
        }
    }

    # Draw a legend
    draw.node.leg <- !is.null(bg.colramp.nodes) && legend.nodes
    draw.branch.leg <- !is.null(col.branch.colramp) && legend.branch
    if (draw.node.leg && draw.branch.leg) {
        .radial.phylog.legend(bg.colramp.nodes, node.ramp.start, node.ramp.end, legend.title.nodes, 0.1, 0.4)
        .radial.phylog.legend(col.branch.colramp, branch.ramp.start, branch.ramp.end, legend.title.branch, 0.6, 0.9)
    }
    else if (draw.node.leg) {
        .radial.phylog.legend(bg.colramp.nodes, node.ramp.start, node.ramp.end, legend.title.nodes)
    }
    else if (draw.branch.leg) {
        .radial.phylog.legend(col.branch.colramp, branch.ramp.start, branch.ramp.end, legend.title.branch)
    }

    if (!is.null(pdf.file)) {
        dev.off()
    }

    return(invisible(retval))
}

.alternating <- function(n) {
    if (n %% 2 == 0) {
        x <- (n / 2) - 0.5
    }
    else {
        x <- floor(n / 2)
    }
    seq(-x, x)
}

.radial.phylog.legend <- function(grad, ramp.start, ramp.end, legend.title, x1=0.3, x2=0.7) {
    par(srt=0)
    ncol <- length(grad)
    y2 <- grconvertY(0, 'nfc', 'user')

    par(xpd=NA)
    x1 <- grconvertX(x1, 'ndc', 'user')
    x2 <- grconvertX(x2, 'ndc', 'user')
    x <- seq(x1, x2, (x2 - x1) / ncol)
    y1 <- grconvertY(0, 'ndc', 'user')
    y2 <- y1 + (0.6 * (y2 - y1))

    rect(x[1:(length(x)-1)], y1, x[2:length(x)], y2, col=grad, border=NA)
    ymid <- y1 + ((y2 - y1) / 2)
    text(x1, ymid, as.character(ramp.start), pos=2, cex=0.8)
    text(x2, ymid, as.character(ramp.end), pos=4, cex=0.8)
    text(x1 + ((x2 - x1) / 2), y2, legend.title, pos=3, cex=0.8)
}