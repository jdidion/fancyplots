# Plot H vs N. If fancy=TRUE, contour lines will be plotted (this can take a while to compute).
plot.HN <- function(HN, fancy=FALSE, k=11, H.lim=NULL, N.lim=NULL, col="red", legend.args=NULL,
        ylab="N", ...) {
    if (is.null(H.lim)) {
        H.lim <- c(0, max(HN[,1]))
    }
    if (is.null(N.lim)) {
        N.lim <- c(0, max(HN[,2]))
    }
    H.mean <- mean(HN[,1])
    N.mean <- mean(HN[,2])
    plot(0:1, type="n", xlab="H", ylab=ylab, xlim=H.lim, ylim=N.lim, ...)
    if (fancy) {
        my.cols <- rev(brewer.pal(k, "RdYlBu"))
        ## compute 2D kernel density
        z <- kde2d(HN[,1], HN[,2], n=50)
        # Make the base plot
        points(HN, pch=19, cex=0.4, col=col)
        # Draw the colored contour lines
        contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
        # Add lines for the mean of X and Y
        abline(h=N.mean, v=H.mean, col="gray", lwd=1.5)
        # Add the correlation coefficient to the top left corner
        legend("topleft", paste("R=", round(cor(HN)[1,2],3)), bty="n")
    }
    else {
        segments(H.mean, N.lim[1], H.mean, N.lim[2], lty=2)
        text(H.mean, N.lim[2] * 0.99, round(H.mean, 3), pos=4)
        segments(H.lim[1], N.mean, H.lim[2], N.mean, lty=2)
        text(H.lim[2] * 0.98, N.mean, round(N.mean, 3), pos=3)
        points(HN, col=rgb.alpha(col, 0.2), pch=20)
    }
    if (!is.null(legend.args)) {
        do.call(legend, legend.args)
    }
}

# If hist=TRUE, histograms will be drawn, otherwise boxplots.
plot.HN.with.dist <- function(HN, H.lim=c(0,max(HN[,1])), N.lim=c(0,max(HN[,2])), hist=FALSE, ...) {
    par(mar=c(4,4,1,1), oma=c(0,0,0,0))
    layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE), widths=c(2.5,1), heights=c(2.5,1))
    plot.HN(HN, H.lim=H.lim, N.lim=N.lim, ylab="", ...)
    result <- list()
    if (hist) {
        result$nhist <- hist(HN[,2], breaks=seq(0, N.lim[2] + 0.01, 0.01), plot=FALSE)
        result$hhist <- hist(HN[,1], breaks=seq(0, H.lim[2] + 0.01, 0.01), plot=FALSE)
        plot(0:1, type="n", xlab="", ylab="N", xlim=c(0, max(result$nhist$counts)), ylim=N.lim,
                yaxt="n")
        rect(0, result$nhist$breaks[1:(length(result$nhist$breaks) - 1)],
             result$nhist$counts, result$nhist$breaks[2:length(result$nhist$breaks)])
        my <- max(result$hhist$counts)
        ms <- seq(0, -my, length.out=7)
        plot(0:1, type="n", xlab="", ylab="", xlim=H.lim, ylim=c(-my, 0),
            xaxt="n", yaxt="n")
        axis(2, at=ms, labels=as.character(-ms))
        rect(result$hhist$breaks[1:(length(result$hhist$breaks) - 1)], 0,
             result$hhist$breaks[2:length(result$hhist$breaks)], -result$hhist$counts)
    }
    else {
        result$nbox <- boxplot(HN[,2], ylab="N", ylim=N.lim)
        result$hbox <- boxplot(HN[,1], horizontal=TRUE, ylim=H.lim)
    }
    invisible(result)
}

# sample.HN is a five-column matrix: sample name, <strain1>, <strain2>, H, N.
# control.HN is the same, but can only have two rows. The first is the left-side
# control, the second is the right-side control.
plot.conc.HN <- function(sample.HN, control.HN=NULL, bar.width=0.01) {
    sample.HN <- sample.HN[order(sample.HN[,3]),]
    rng <- range(sample.HN[,3])
    xlim <- c(rng[1] - (bar.width * 5), rng[2] + (bar.width * 5))
    ymax <- max(apply(sample.HN[,4:5],1,sum)) * 1.1

    layout(matrix(c(1,2), nrow=2), heights=c(1,3))

    par(mar=c(2,3,2,3), oma=c(6,0,0,0))
    plot.new()
    plot.window(xlim=xlim, ylim=c(0,1))
    polygon(c(rng[1], rng[2], rng[2], rng[1]), c(1, 1, rng[2], rng[1]), col="blue")
    polygon(c(rng[1], rng[2], rng[2], rng[1]), c(0, 0, rng[2], rng[1]), col="yellow")
    text(rng[1]-0.02, rng[1] + ((1-rng[1])/2), colnames(sample.HN)[2], srt=90, xpd=NA)
    text(rng[2]+0.02, rng[2] / 2, colnames(sample.HN)[3], srt=-90, xpd=NA)
    axis(3,c(rng[1],rng[2]),c(1-rng[1],1-rng[2]))
    axis(1,c(rng[1],rng[2]),rng)

    #par(mar=c(6,3,2,3), oma=c(0,0,0,0))
    plot(0:1, type="n", xlim=xlim, ylim=c(0, ymax), xaxt="n")
    if (is.null(control.HN)) {
        axis(1, sample.HN[,3], sample.HN[,1], las=3, cex.axis=0.5)
    }
    else {
        axis(1, c(xlim[1] + (bar.width * 1.5), sample.HN[,3], xlim[2] - (bar.width * 1.5)),
                c(control.HN[1,1], sample.HN[,1], control.HN[2,1]), las=3, cex.axis=0.5)
        rect(xlim[1] + bar.width, 0, xlim[1] + (2 * bar.width), control.HN[1,4], col="gray")
        rect(xlim[1] + bar.width, control.HN[1,4], xlim[1] + (2 * bar.width),
            control.HN[1,4] + control.HN[1,5], col="green")
        rect(xlim[2] - bar.width, 0, xlim[2] - (2 * bar.width), control.HN[2,4], col="gray")
        rect(xlim[2] - bar.width, control.HN[2,4], xlim[2] - (2 * bar.width),
            control.HN[2,4] + control.HN[2,5], col="green")
    }

    rect(sample.HN[,3] - (bar.width / 2), 0,
        sample.HN[,3] + (bar.width / 2), sample.HN[,4], col="gray")
    rect(sample.HN[,3] - (bar.width / 2), sample.HN[,4],
        sample.HN[,3] + (bar.width / 2), sample.HN[,4] + sample.HN[,5], col="green")

    legend("top", c("H","N"), fill=c("gray","green"), xpd=NA, xjust=0.5, horiz=TRUE)
}

plot.intensities <- function(intens, dim=c(3,3), ...) {
    n <- dim[1] * dim[2]
    ks <- sapply(1:ncol(intens), function(i) {
        if (i %% n == 1) {
            par(mfrow=dim, mar=c(2.5,2,1,1))
        }
        plot.intensity(colnames(intens)[i], intens[,i], ...)
    })
    names(ks) <- colnames(intens)
    ks
}

plot.intensity <- function(main, intens, ylim=NULL, control.mean=1.099, control.sd=0.443, legend=FALSE) {
    h <- hist(intens, breaks=seq(round(min(intens), 1) - 0.1, round(max(intens), 1) + 0.1, 0.1),
        plot=FALSE)
    n <- length(h$breaks)
    y <- h$density
    xlim <- c(h$breaks[1], h$breaks[n])
    if (is.null(ylim)) {
        ylim <- c(0, max(y))
    }
	ks <- ks.test(intens, "pnorm", mean=control.mean, sd=control.sd)
	plot(0:1, type="n", xlim=xlim, ylim=ylim, xlab="Intensity", ylab="Fraction of SNPs", main=main)
    rect(h$breaks[1:(n-1)], 0, h$breaks[2:n], y, col="lightgray")
    x <- seq(xlim[1], xlim[2], 0.1)
    lines(x, dnorm(x, mean=control.mean, sd=control.sd), col="blue", lwd=2)
    lines(h$mids, y, col="red", lwd=2)
    text(xlim[2] * 0.9, ylim[2], paste("KS =", round(ks$statistic,4)), pos=2)
    if (legend) {
        legend("topright", legend=c("Controls", "Sample"), fill=c("blue", "red"))
    }
    return(ks)
}

plot.N.markers <- function(snps, HN, fun=mean, main=NULL) {
    N <- values.as.bins(cbind(snps[,1:2], HN[,2]), fun, rename.chrm=TRUE)
    N[is.na(N)] <- 0
    plot.genome.data(N, type="bar", gp=gpar(col=NA, fill="black"), yscale=TRUE, title=main)
}
