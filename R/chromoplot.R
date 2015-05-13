split.gp <- function(gp, n) {
    gps <- list()
    for (i in 1:n) {
        l <- list()
        for (name in names(gp)) {
            if (length(gp[[name]]) > 1) {
                l[[name]] <- gp[[name]][i]
            }
            else {
                l[[name]] <- gp[[name]]
            }
        }
        gps[[i]] <- l
    }
    gps
}

clip.values <- function(v, lim) {
    if (is.vector(v)) {
        .do.clip(v, lim)
    }
    else {
        apply(v, 2, .do.clip, lim)
    }
}

.do.clip <- function(v, lim) {
    v[v < lim[1]] <- lim[1]
    v[v > lim[2]] <- lim[2]
    v
}

seq.units <- list(Mb=1000000, kb=100000, bp=1)
plot.functions <- list(
    bar=function(d, gp, ylim) grid.rect(d$x1, 0, d$x2 - d$x1, clip.values(d$y, ylim),
        default.units='native', just=c('left','bottom'), gp=gp),
    stackedbar=function(d, gp, ylim) {
        cols <- names(d)
        ycols <- sort(cols[grep('y\\d', cols, perl=TRUE)])
        n <- length(ycols)
        gps <- split.gp(gp, n)
        y1 <- rep(0, nrow(d))
        for (i in 1:n) {
            height <- clip.values(d[,ycols[i]], ylim)
            g <- do.call(gpar, gps[[i]])
            grid.rect(d$x1, y1, d$x2 - d$x1, height, default.units='native', just=c('left','bottom'), gp=g)
            y1 <- y1 + height
        }
    },
    lines=function(d, gp, ylim) {
        cols <- names(d)
        ycols <- sort(cols[grep('y\\d', cols, perl=TRUE)])
        n <- length(ycols)
        gps <- split.gp(gp, n)
        for (i in 1:n) {
            g <- do.call(gpar, gps[[i]])
            grid.lines(d$x, clip.values(d[,ycols[i]], ylim), default.units='native', gp=g)
        }
    },
    markers=function(d, gp, ylim) {
        if ("y" %in% names(d)) {
            y <- d$y
        }
        else {
            y <- 1
        }
        grid.segments(d$x, 0, d$x, y, default.units='native', gp=gp)
    },
    blocks=function(d, gp, ylim) {
        fills <- sort(unique(d$fill))
        n <- length(fills)
        if (length(gp) == 0) {
            g <- gpar()
        }
        else {
            gps <- split.gp(gp, n)
        }
        for (i in 1:n) {
            f <- fills[i]
            cur.d <- d[d$fill == f,,drop=FALSE]
            if (length(gp) > 0) {
                g <- do.call(gpar, gps[[i]])
            }
            g$fill <- f
            if (!("col" %in% names(g))) {
                g$col <- NA
            }
            grid.rect(cur.d$x1, 0, cur.d$x2 - cur.d$x1, clip.values(cur.d$y, ylim),
                default.units='native', just=c('left','bottom'), gp=g)
        }
    },
    point=function(d, gp, ylim) {
        pch <- 1
        size <- unit(1, 'char')
        if ('pch' %in% colnames(d)) {
            pch <- d$pch
        }
        if ('size' %in% colnames(d)) {
            size <- unit(d$size, 'char')
        }
        if ('col' %in% colnames(d)) {
            col <- d$col
            if (is.null(gp)) {
                gp <- gpar(col=col)
            }
            else {
                gp$col <- col
            }
        }
        grid.points(d$x, clip.values(d$y, ylim), default.units='native', pch=pch, size=size, gp=gp)
    },
    meanline=function(d, gp, ylim) {
        m <- mean(d$y)
        grid.segments(min(d$x), m, max(d$x), m, default.units='native', gp=gp)
    },
    curve=function(d, gp, ylim) {
        ss <- supsmu(d$x, d$y)
        grid.lines(ss$x, clip.values(ss$y, ylim), default.units='native', gp=gp)
    },
    segs=function(d, gp, ylim) {
        grid.segments(d$x1, clip.values(d$y, ylim), d$x2, clip.values(d$y, ylim), default.units='native', gp=gp)
    },
    steps=function(d, gp, ylim) {
        cols <- names(d)
        ycols <- sort(cols[grep('y\\d', cols, perl=TRUE)])
        n <- length(ycols)
        gps <- split.gp(gp, n)
        for (i in 1:n) {
            height <- clip.values(d[,ycols[i]], ylim)
            g <- do.call(gpar, gps[[i]])
            grid.segments(d$x1, height, d$x2, height, default.units='native', gp=g)
            grid.segments(d[1:(nrow(d)-1),'x2'], height[1:(nrow(d)-1)],
                          d[1:(nrow(d)-1),'x2'], height[2:nrow(d)], default.units='native', gp=g)
        }
    },
    jitter=function(d, gp, ylim) {
        pch <- 1
        size <- unit(1, 'char')
        if ('pch' %in% colnames(d)) {
            pch <- d$pch
        }
        if ('size' %in% colnames(d)) {
            size <- unit(d$size, 'char')
        }
        mid <- sum(ylim) / 2
        y <- jitter(rep(mid, nrow(d)), amount=mid - ylim[1])
        grid.points(d$x, clip.values(y, ylim), default.units='native', pch=pch, size=size, gp=gp)
    }
)

resolve.func <- function(name) {
    func <- name
    if (!is.null(name) && !is.function(name)) {
        func <- plot.functions[[name]]
    }
    if (is.null(func)) {
        stop(paste("Invalid plot type", name))
    }
    func
}

draw.genome <- function(chromosomes, karyotype, ylim, title=NULL, xint=25, xunit='Mb',
        ylab.right=NULL, yscale=FALSE, yscale.label=TRUE, ygrid=FALSE, chrm.layout=NULL) {
    n <- length(chromosomes)
    hsize <- c(3, 1, 2)
    vsize <- c(0, 1, 4)
    if (!is.null(title))
        vsize[1] <- 2
    if (!is.null(ylab.right))
        hsize[3] <- 1
    #if (yscale)
    #    hsize[3] <- hsize[3] + 2
    pushViewport(viewport(layout=grid.layout(3, 3,
        heights=unit(vsize, c('lines', 'null', 'lines')),
        widths=unit(hsize, c('lines', 'null', 'lines')))))
    if (!is.null(title))
        grid.text(title, y=unit(convertY(unit(1, 'npc'), 'lines', TRUE) - 1, 'lines'))
    grid.text(paste("Genomic Position (", xunit, ")"), y=unit(1, 'lines'))
    grid.text("Chromosome", x=unit(1, 'lines'), rot=90)
    if (!is.null(ylab.right))
        grid.text(ylab.right, x=unit(convertX(unit(1, 'npc'), 'lines', TRUE) - 1, 'lines'), rot=90)

    sizes <- karyotype$chrm[chromosomes,'size']
    xint <- xint * seq.units[[xunit]]
    xmax <- as.integer(ceiling(max(sizes) / xint) * xint)
    xticks <- seq(0, xmax, xint)
    ypad <- (ylim[2] - ylim[1]) * 0.1 # amount of padding between plot and chromosome borders
    ybox <- c(ylim[1] - ypad, ylim[2] + ypad)
    pushViewport(viewport(layout=grid.layout(n, 1), xscale=c(0, xmax), layout.pos.row=2, layout.pos.col=2, name='data'))
    if (ygrid) {
        grill <- unit(xticks[2:(length(xticks)-1)], 'native')
        grid.segments(grill, 0, grill, 1, gp=gpar(lty=2))
    }

    for (i in 1:n) {
        chrm <- chromosomes[i]
        size <- sizes[i]
        pushViewport(viewport(layout=grid.layout(3, 1, heights=c(0.1, 0.8, 0.1)),
            layout.pos.row=i, layout.pos.col=1))
        pushViewport(viewport(layout=chrm.layout, layout.pos.row=2, layout.pos.col=1,
            xscale=c(0, xmax), yscale=ybox, name=chrm))
        grid.rect(x=0, width=unit(size, 'native'), just='left', gp=gpar(col='white', fill='white', alpha=0.8))
        grid.segments(0, ybox, size, ybox, 'native')
        grid.curve(0, ybox[1], 0, ybox[2], 'native', -1, square=FALSE, ncp=100)
        grid.curve(size, ybox[1], size, ybox[2], 'native', 0.7, square=FALSE, ncp=100)
        grid.text(karyotype$chrm[chrm, 'num'], x=unit(-0.5, 'lines'), rot=90, just='top', gp=gpar(fontsize=9))
        if (i == 1 && (!is.logical(yscale) || yscale)) {
            yl <- ylim
            if (!is.logical(yscale)) {
                yl <- yscale
            }
            grid.yaxis(yl, label=yscale.label, main=FALSE, gp=gpar(fontsize=9))
        }
        upViewport()

        if (i == n) {
            pushViewport(viewport(layout.pos.row=3, layout.pos.col=1, xscale=c(0, xmax), yscale=ybox))
            grid.xaxis(xticks, as.character(round(xticks / seq.units[[xunit]])), gp=gpar(fontsize=9))
        }

        upViewport()
    }

    upViewport()
}

# Takes a data frame where rownames are labels and columns are col, pch and lty
draw.legend <- function (legend, gp=gpar(), vp=NULL) {
    labels <- rownames(legend)
    nkeys <- length(labels)

    legend.layout <- grid.layout(nkeys, 3,
        widths=unit.c(unit(1, "lines"),
            max(unit(rep(1, nkeys), "strwidth", as.list(labels))), unit(0.5, "lines")),
        heights=unit.pmax(unit(1, "lines"),
            unit(0.5, "lines") + unit(rep(1, nkeys), "strheight", as.list(labels))))

    fg <- frameGrob(layout=legend.layout, vp=vp, gp=gp)

    for (i in 1:nkeys) {
        col <- 'black'
        if ('col' %in% names(legend)) {
            col <- legend[i, 'col']
        }

        pch <- 'pch' %in% names(legend)
        lty <- 'lty' %in% names(legend)

        if (pch) {
            fg <- placeGrob(fg, pointsGrob(0.5, 0.5, pch=legend[i, 'pch'],
                gp=gpar(col=col, cex=0.75)), col=1, row=i)
        }

        if (lty) {
            fg <- placeGrob(fg, linesGrob(x=c(0.1, 0.9), y=0.5,
                gp=gpar(lty=legend[i,'lty'], col=col)), col=1, row=i)
        }

        if (!pch && !lty) {
            fg <- placeGrob(fg, rectGrob(0.5, 0.5, 0.8, 0.8, gp=gpar(fill=col)), col=1, row=i)
        }

        fg <- placeGrob(fg, textGrob(labels[i], x=0, y=0.5, just=c("left", "center"), gp=gpar(cex=0.75)), col= 2, row=i)
    }

    grid.draw(fg)
}

read.karyotype <- function(kfile='~/data/karyotype.mouse.mm9.txt') {
    tab <- read.table(kfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    chrms <- tab[,1] == 'chr'
    bands <- tab[,1] == 'band'
    return(list(
        chrm=data.frame(num=tab[chrms,4], size=tab[chrms,6], name=tab[chrms,7], row.names=tab[chrms,3]),
        band=lapply(by(tab[bands,], tab[bands,2], identity), function(x)
            data.frame(id=x[,4], start=x[,5], end=x[,6], name=x[,7], row.names=x[,3]))))
}

# Converts a two-column data frame (chromosome, position) or a list
# (names=chromosome, values=vector of positions) into input for plot.genome.data.
# e.g. plot.genome.data(as.bins(data, binsize), 'bar', ...)
positions.as.bins <- function(data, binsize=500000, rename.chrm=FALSE) {
    if (is.data.frame(data)) {
        bins <- by(floor(data[,2] / binsize) * binsize, data[,1],
            function(x) unlist(by(x, x, length, simplify=FALSE)))
    }
    else {
        bins <- lapply(data, function(x) table(floor(x / binsize) * binsize))
    }
    bins <- compact(bins)
    df <- do.call(rbind, lapply(names(bins), function(n) {
        x <- bins[[n]]
        start <- as.integer(names(x))
        data.frame(chrm=n, x1=start, x2=start + binsize - 1, y=as.vector(x))
    }))
    if (rename.chrm) {
        df$chrm <- paste("mm", df$chrm, sep="")
    }
    df
}

# Converts a three-column data frame (chromosome, position, value) into input for
# plot.genome.data, where values are aggregated using fun.
values.as.bins <- function(data, fun, binsize=500000, rename.chrm=FALSE) {
    bins <- by(data.frame(floor(data[,2] / binsize) * binsize, data[,3], stringsAsFactors=FALSE),
        data[,1], function(x) tapply(x[,2], x[,1], fun))
    bins <- compact(bins)
    df <- do.call(rbind, lapply(names(bins), function(n) {
        x <- bins[[n]]
        start <- as.integer(names(x))
        data.frame(chrm=n, x1=start, x2=start + binsize - 1, y=as.vector(x))
    }))
    if (rename.chrm) {
        df$chrm <- paste("mm", df$chrm, sep="")
    }
    df
}

#' Plot data on one or more chromosomes as a function of genomic position. Mapping must map columns
#' of data to the inputs required by plotting functions (x1, x2, y1, etc). type must one of the
#' plotting functions in .plot.functions or a function that accepts two arguments:
#' a data frame (subset of data for a single chromosome) and a gpar.
#' NOTE: It is strongly recommended to use a non-pdf device type (specified using the device
#' parameter) when data is a list due to the graphics antialiasing that most viewers do by default.
plot.genome.data <- function(data, mapping=NULL, types=c('point'), karyotype=NULL, chromosomes=NULL, ylim=NULL,
        gps=list(gpar()), legend=NULL, legend.vp=NULL, device="pdf", dev.file=NULL, dev.width=7, dev.height=7,
        yscale=FALSE, yscale.label=TRUE, sample.borders=FALSE, ...) {
    islist <- !is.data.frame(data) && is.list(data)

    if (is.null(karyotype)) {
        karyotype <- read.karyotype()
    }
    else if (is.character(karyotype)) {
        karyotype <- read.karyotype(karyotype)
    }

    if (!is.null(mapping)) {
        if (islist) {
            for (i in 1:length(data)) {
                data[[i]]$data <- eval.mapping(data[[i]]$data, mapping)
            }
        }
        else {
            data <- eval.mapping(data, mapping)
        }
    }

    if (is.null(chromosomes)) {
        if (islist) {
            ch <- unique(unlist(lapply(data, function(x) x$data$chrm)))
        }
        else {
            ch <- data$chrm
        }
        chromosomes <- intersect(rownames(karyotype$chrm), ch)
    }
    
    if (length(chromosomes) == 0) {
        stop("No chromosomes to plot. Make sure rownames in kayotype match chromosome names in data.")
    }

    funcs <- sapply(types, resolve.func)
    names(funcs) <- types
    chrm.layout <- NULL

    if (islist) {
        n <- length(data)
        heights <- unlist(lapply(data, function(x) x$height))
        w <- is.null(heights)
        s <- sum(heights)
        if (is.null(s)) s <- 0
        if (any(w)) {
            heights[w] <- (1.0 - s) / sum(w)
        }
        s <- sum(heights)
        # scale so that sum = 1.0
        if (s != 1.0) {
            heights <- heights * (1.0 / s)
        }
        chrm.layout <- grid.layout(n, 1, heights=heights)

        for (i in 1:length(data)) {
            # TODO: correctly calculate ylim for stackedbar
            if (!('ylim' %in% names(data[[i]]))) {
                if (!is.null(ylim)) {
                    data[[i]]$ylim <- ylim
                }
                else if ('y' %in% names(data[[i]]$data)) {
                    data[[i]]$ylim <- c(min(data[[i]]$data$y), max(data[[i]]$data$y))
                }
                else {
                    data[[i]]$ylim <- c(0,1)
                }
            }

            if ('yscale' %in% names(data[[i]]) && data[[i]]$yscale) {
                cs <- n * (1 - cumsum(heights))
                yscale <- c(cs[i], ifelse(i == 1, n, cs[i-1]))
                if (yscale.label) {
                    yscale.label <- data[[i]]$ylim
                }
            }

            if (!('gps' %in% names(data[[i]]))) {
                data[[i]]$gps <- gps
            }

            if ('funcs' %in% names(data[[i]])) {
                data[[i]]$funcs <- resolve.func(data[[i]]$func)
            }
            else {
                data[[i]]$funcs <- funcs
            }
        }

        ylim <- c(0, length(data))
    }
    else {
        if (is.null(funcs)) {
            stop(paste("Invalid plot type", type))
        }

        if (is.null(ylim)) {
            if ('y' %in% names(data)) {
                ylim <- c(min(data$y), max(data$y))
            }
            else {
                ylim <- c(0,1)
            }
        }
    }

    if (!is.null(dev.file)) {
        if (device == "pdf") {
            pdf(dev.file, dev.width, dev.height)
        }
        else {
            match.fun(device)(dev.file, dev.width, dev.height, units="in", res=300)
        }
        on.exit(dev.off())
    }

    grid.newpage()
    draw.genome(chromosomes, karyotype, ylim, chrm.layout=chrm.layout, yscale=yscale, yscale.label=yscale.label, ...)

    for (chrm in chromosomes) {
        downViewport(chrm)

        if (islist) {
            for (i in 1:length(data)) {
                w <- data[[i]]$data$chrm == chrm
                if (any(w)) {
                    parent <- current.viewport()
                    pushViewport(viewport(layout.pos.row=i, layout.pos.col=1,
                        xscale=parent$xscale, yscale=data[[i]]$ylim, name=paste(chrm,i)))
                    if (sample.borders) {
                        grid.rect(0, data[[i]]$ylim[1], karyotype$chrm[chrm,'size'],
                            data[[i]]$ylim[2] - data[[i]]$ylim[1], gp=gpar(col="lightgray", lwd=0.1),
                            default.units="native", just=c("left","bottom"))
                    }
                    for (j in 1:length(data[[i]]$funcs)) {
                        fname <- names(data[[i]]$funcs)[j]
                        f <- data[[i]]$funcs[[fname]]
                        gp <- data[[i]]$gp
                        if (fname %in% names(gp)) {
                            gp <- gp[[fname]]
                        }
                        f(data[[i]]$data[w,], gp, data[[i]]$ylim)
                    }
                    upViewport()
                }
            }
        }
        else {
            w <- data$chrm == chrm
            if (any(w)) {
                for (i in 1:length(funcs)) {
                    fname <- names(funcs)[[i]]
                    f <- funcs[[fname]]
                    gp <- gps
                    if (fname %in% names(gp)) {
                        gp <- gp[[fname]]
                    }
                    else {
                        gp <- gp[[1]]
                    }
                    f(data[w,], gp, ylim)
                }
            }
        }

        seekViewport('data')
    }

    if (!is.null(legend)) {
        if (is.null(legend.vp)) {
            legend.vp <- viewport(unit(0.8, 'npc'), unit(0.1, 'npc'))
        }
        draw.legend(legend, vp=legend.vp)
    }
}

# Plot a summary of a single variable across the genome for one or more samples. The data matrix
# has one row for each sample and one column for each possible value of the variable. Cell values
# are the fraction of the genome for which the given sample has the given value. Colors is a vector
# of one color per column.
plot.genome.summary <- function(data, colors, title=NULL, row.pad=0, pdf.file=NULL, pdf.width=7, pdf.height=7) {
    if (!is.null(pdf.file)) {
        pdf(pdf.file, pdf.width, pdf.height)
        on.exit(dev.off())
    }

    n <- nrow(data)
    hsize <- c(5, 1, 1)
    vsize <- c(0, 1, 3)
    if (!is.null(title))
        vsize[1] <- 2
    pushViewport(viewport(layout=grid.layout(3, 3,
        heights=unit(vsize, c('lines', 'null', 'lines')),
        widths=unit(hsize, c('lines', 'null', 'lines')))))
    if (!is.null(title))
        grid.text(title, y=unit(convertY(unit(1, 'npc'), 'lines', TRUE) - 1, 'lines'))
    grid.text(paste("Fraction of Genome"), y=unit(1, 'lines'))

    # draw sample labels
    pushViewport(viewport(layout=grid.layout(n, 1), xscale=c(0, 1), layout.pos.row=2, layout.pos.col=1, name='data'))
    for (i in 1:n) {
        pushViewport(viewport(layout=grid.layout(1, 1, heights=1),
            layout.pos.row=i, layout.pos.col=1, xscale=c(0, 1), yscale=c(0, 1)))
        name <- rownames(data)[i]
        grid.text(name, y=0, gp=gpar(cex=0.7))
        upViewport()
    }
    upViewport()

    # draw sample data
    pushViewport(viewport(layout=grid.layout(n, 1), xscale=c(0, 1), layout.pos.row=2, layout.pos.col=2, name='data'))
    for (i in 1:n) {
        pushViewport(viewport(layout=grid.layout(1, 1, heights=1),
            layout.pos.row=i, layout.pos.col=1, xscale=c(0, 1), yscale=c(0, 1)))
        x <- 0
        for (j in 1:ncol(data)) {
            grid.rect(
                x=unit(x, "native"), width=unit(data[i,j], "native"),
                y=unit(row.pad, "native"), height=unit(1 - (2 * row.pad), "native"),
                just='left', gp=gpar(col=NA, fill=colors[j]))
            x <- x + data[i,j]
        }
        upViewport()
    }
    upViewport()
}
