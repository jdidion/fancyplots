# Plot BAF and LRR across the genome (or the specified chromosomes) for the specified samples.
# If A and B alleles are not specified for snps, the BAF plot will only be grayscale. Extra
# arguments are passed directly to plot.sample.copy.number.metrics(). Plotting of LRR and
# Z-scores can be controlled by plot.lrr and plot.z params. If density.band.window is set to
# NULL, trend lines will be plotted instead of density bands.
plot.copy.number.metrics <- function(metrics, img.path, img.format="jpeg", sampleIDs=NULL, ...) {
    if (is.null(sampleIDs)) {
        sampleIDs <- names(metrics)
    }
    # If we're creating a PDF, create a single document with each plot a separate page
    if (img.format == "pdf") {
        pdf(img.path, 16, 12)
        for (sample in sampleIDs) {
            plot.sample.copy.number.metrics(metrics[[sample]], ...)
        }
        dev.off()
    }
    # Otherwise create one file per plot
    else {
        img.func <- get(img.format)
        for (sample in sampleIDs) {
            name <- gsub("[\\/%]", "_", metrics[[sample]]$name)
            img.func(file.path(img.path, paste(name, img.format, sep=".")), 16, 12, units="in", res=300)
            plot.sample.copy.number.metrics(metrics[[sample]], ...)
            dev.off()
        }
    }
}

# Plot BAF, LRR and Z-scores for a single sample.
plot.sample.copy.number.metrics <- function(metrics, ...) {
    # Only plot SNPs for which we have metrics
	snps <- metrics$snps
	mets <- metrics$metrics
    m <- match(rownames(snps), rownames(mets))
	snps <- snps[!is.na(m),]
    mets <- mets[m[!is.na(m)],]
    data <- data.frame(chrm=snps$chrm, pos=snps$pos, call=mets$call, 
        baf=mets$baf, lrr=mets$lrr, z=mets$z, row.names=rownames(snps), stringsAsFactors=FALSE)
    plot.sample.copy.number.data(metrics$name, data, ...)
}

plot.sample.copy.number.data <- function(name, data, chromosomes=c(1:19,"X","Y"),
        zoom=NULL, tick.interval=NULL, plot.baf=TRUE, plot.lrr=TRUE, lrr.ylim=c(-2,2),
        plot.z=FALSE, z.ylim=c(-1, 1), z.d=5, z.r=4, z.intervals=1000.0,
        plot.trend=TRUE, density.band.window=100, density.band.frac=0.75, 
        other.plots=NULL, other.data=NULL, ...) {

    # Only plot SNPs on specified chromosomes
    if (!is.null(chromosomes)) {
        data <- data[!is.na(match(data$chrm, chromosomes)),]
    }
    chromosomes <- unique(data$chrm)

    # Convert chromosome positions to absolute genomic positions
    chrms <- factor(data$chrm, levels=chromosomes)
    chrm.sizes <- as.numeric(tapply(data$pos, chrms, max))
    genome.size <- sum(chrm.sizes)
    names(chrm.sizes) <- chromosomes
    cs <- c(0, cumsum(as.numeric(chrm.sizes[-length(chrm.sizes)])))
    names(cs) <- chromosomes
    x <- cs[as.character(data$chrm)] + data$pos

    xlim <- c(0, genome.size)
    if (length(chromosomes) == 1 && !is.null(zoom)) {
        xlim <- zoom
    }

    nplots <- plot.baf + plot.lrr + plot.z
    if (!is.null(other.plots)) {
        nplots <- nplots + length(other.plots)
    }
    if (nplots > 1) {
        par(mfrow=c(nplots, 1), mar=c(0.25,5,0.25,1), oma=c(3, 0, 3, 0))
    }

    if (plot.baf) {
        .plot.baf(x, data$baf, data$call, chrms, xlim=xlim, ...)
    }

    if (plot.lrr) {
        .plot.scatter.with.trend.and.outliers(x, data$lrr, chrms, lrr.ylim, xlim=xlim,
            ylab="Log R Ratio", density.band.window, density.band.frac, plot.trend, ...)
    }

    if (plot.z) {
        z <- hyperlog(data$z, b=0, d=z.d, r=z.r, intervals=z.intervals)$y
        .plot.scatter.with.trend.and.outliers(x, z, chrms, z.ylim, xlim=xlim,
            xlab="Position", ylab="HL(Z-score)", density.band.window, density.band.frac, plot.trend, ...)
    }

    if (!is.null(other.plots)) {
        for (fn in other.plots) {
            fn(data, other.data, chromosomes=cs, xlim=xlim, ...)
        }
    }

    if (is.null(tick.interval)) {
        tick.pos <- cs[1:length(chromosomes)]
        tick.name <- chromosomes
    }
    else {
        tick.pos <- NULL
        tick.name <- NULL
        for (i in 1:length(cs)) {
            ticks <- seq(1, chrm.sizes[i], tick.interval)
            tick.pos <- c(tick.pos, cs[i] + ticks)
            tick.name <- c(tick.name, chromosomes[i], ticks[2:length(ticks)] / 1000000)
        }
    }

    if (!is.null(name)) {
        title(name, outer=TRUE)
    }

    axis(1, tick.pos, tick.name)
}

.plot.baf <- function(x, baf, calls, chrms, A.col=c("lightblue","darkblue"), B.col=c("pink", "red"),
        H.col=c("mediumpurple1", "purple4"), N.col=c("lightgray", "darkgray"), mid.col='black', ...) {
    plot(0:1, type="n", ylim=c(0, 1), xaxt="n", xlab="", ylab="B Allele Frequency", ...)
    segments(0, 0.5, max(x), 0.5, col=mid.col, lty=2)

    odd <- as.integer(chrms) %% 2 == 1
    even <- as.integer(chrms) %% 2 == 0

    col <- rep(N.col[1], length(x))

    A <- calls == "A"
    B <- calls == "B"
    H <- calls == "H"

    col[even] <- N.col[2]
    col[odd & A] <- A.col[1]
    col[even & A] <- A.col[2]
    col[odd & H] <- H.col[1]
    col[even & H] <- H.col[2]
    col[odd & B] <- B.col[1]
    col[even & B] <- B.col[2]

    segments(x[1], 0.5, x[length(x)], 0.5, lty=2)
    points(x, baf, col=col, pch=20)
}

.plot.scatter.with.trend.and.outliers <- function(x, y, chrms, ylim, density.band.window,
        density.band.frac, plot.trend=!is.null(density.band.window),
        point.col=c("lightgray", "darkgray"), zero.col="black",
        trend.col="red", density.fill.col=rgb(1,0,0,0.33), ...) {
    odd <- as.integer(chrms) %% 2 == 0
    even <- as.integer(chrms) %% 2 == 1

    outliers.low <- y < ylim[1]
    y[outliers.low] <- ylim[1]
    outliers.high <- y > ylim[2]
    y[outliers.high] <- ylim[2]

    col <- rep(point.col[1], length(x))
    col[even] <- point.col[2]
    col[outliers.high | outliers.low] <- trend.col

    plot(0:1, type="n", ylim=ylim, xaxt="n", ...)
    points(x, y, col=col, pch=20)

    segments(0, 0, max(x), 0, col=zero.col, lty=2)

    if (plot.trend) {
        # Plot "trend lines" for each chromosome. Currently using the supersmoother function,
        # but could also use a kernel smoother, regression, splines or loess. See:
        # http://statistics.berkeley.edu/classes/s133/Smooth-a.html
        # http://statistics.berkeley.edu/classes/s133/Smooth1a.html
        fits <- by(data.frame(x=x, y=y), chrms, function(m) supsmu(m$x, m$y))
        for (n in names(fits)) {
            lines(fits[[n]]$x, fits[[n]]$y, col=trend.col, lwd=2)
        }
    }
    if (!is.null(density.band.window)) {
        # Use a sliding window to determine the upper and lower bounds of y for density.band.frac SNPs
        bounds <- by(data.frame(x=x, y=y), chrms, function(m) {
            # Calculate the upper and lower quantiles for all windows
            rem <- (1 - density.band.frac) / 2
            if (nrow(m) <= density.band.window) {
                q <- quantile(m$y, c(rem, 1 - rem))
                m$low <- q[1]
                m$high <- q[2]
            }
            else {
                last.win <- nrow(m) - density.band.window + 1
                m$low <- NA
                m$high <- NA
                m[1:last.win, c("low", "high")] <- t(sapply(seq(1, last.win), function(start) {
                    end <- start + density.band.window - 1
                    quantile(m$y[start:end], c(rem, 1 - rem))
                }))
                # Copy the value of the last window to cover the remaining x values
                m[(last.win + 1):nrow(m), c("low", "high")] <- m[last.win, c("low", "high")]
            }
            m
        }, simplify=FALSE)

        # Draw filled band lines
        for (n in names(bounds)) {
            m <- bounds[[n]]
            lines(m$x, m$low, col=trend.col, lwd=0.5)
            lines(m$x, m$high, col=trend.col, lwd=0.5)
            polygon(c(m$x, rev(m$x)), c(m$low, rev(m$high)), border=NA, col=density.fill.col)
        }
    }
}
