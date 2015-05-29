# Functions to plot results from a GWAS studies

#' Plot results of BIMBAM run. Draws manhattan plots based on BF. Expects
#' at least 4 column data frame: pos, chr, bf, pv.
manhattan.bimbam <- function(bb, outdir, main='', outfile='bayes_manhattan.pdf', prior=NULL, 
        threshold=c(0.99, 0.95), as.prob=FALSE, ...) {
    if (is.character(bb)) {
        bb <- read.bimbam(bb)
    }
    dir.create(outdir, showWarnings=FALSE)
    if (is.null(prior)) {
        # use an unbiased prior by default
        prior <- 1 / nrow(bb)
    }
    bb$ppa <- sapply(bb$bf, ppa, prior)
    annotate <- bb[bb$ppa > max(threshold), 'rsnum']
    if (as.prob) {
        manhattan.bayes(bb, mapping(SNP=rsnum, CHR=chr, BP=pos, STAT=ppa), annotate=annotate, 
            threshold=threshold, main=main, pdf.file=paste(outdir, outfile, sep='/'), ...)
    }
    else {
        manhattan.bayes(bb, mapping(SNP=rsnum, CHR=chr, BP=pos, STAT=bf), annotate=annotate, 
            threshold=bf(threshold, prior), main=main, pdf.file=paste(outdir, outfile, sep='/'), ...)
    }
}

#' Manhattan plot for results from Bayesian analysis. Requires at least 3 column data frame: 
#' CHR, BP, and STAT, corresponding to chromosome number, genomic coordinate, and Bayes factor. 
#' Missing values (where regression failed to converge, for instance) should be recoded NA before 
#' reading into R. The plot alternates between two colors.
manhattan.bayes <- function(d, mapping=NULL, colors=c("steelblue", "gray50"), chrms=1:19, 
        sub="Bayesian Association", ...) {
    if (!is.null(mapping)) {
        d <- eval.mapping(d, mapping)
    }
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "STAT" %in% names(d))) {
        stop("Make sure your data frame contains columns CHR, BP, and STAT")
    }
    if (any(chrms)) {
        d <- d[d$CHR %in% chrms, ]
    }
    
    d <- subset(na.omit(d[order(d$CHR, d$BP), ]), (STAT > 0)) # remove na's, sort, and keep only BF > 0
    
    manhattan(d, colors, "log10(BF)", sub=sub, ...)
}
        
#' Manhattan plot for results from frequentist analysis. Requires at least 3 column data frame: 
#' CHR, BP, and P, corresponding to chromosome number, genomic coordinate, and p-value. Missing 
#' values (where regression failed to converge, for instance) should be recoded NA before reading 
#' into R. The plot alternates between two colors.
#'
#' @note Copied and modified from Stephen Turner (http://GettingGeneticsDone.blogspot.com/).
manhattan.freq <- function(d, mapping=NULL, colors=c("steelblue", "gray50"), chrms=1:19, 
        threshold=c(1e-5, 1e-8), sub="Frequentist Association", ...) {
    if (!is.null(mapping)) {
        d <- eval.mapping(d, mapping)
    }
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) {
        stop("Make sure your data frame contains columns CHR, BP, and P")
    }
    if (any(chrms)) {
        d <- d[d$CHR %in% chrms, ]
    }
    
    d <- subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$STAT <- -log10(d$P)
    
    if (any(threshold < 1)) {
        threshold <- -log10(threshold)
    }
    
    manhattan(d, colors, expression(-log[10](italic(p))), threshold=threshold, sub=sub, ...)
}

cycle <- function(v, i) {
    l <- length(v)
    x <- i %% l
    if (x == 0) {
        v[l]
    }
    else {
        v[x]
    }
}

alternate.colors <- function(d, colors) {
    chrs <- unique(d$CHR)
    colors <- rep(colors, length.out=length(chrs))
    col <- rep('black', nrow(d)) 
    for (i in 1:length(chrs)) {
        col[d$CHR==chrs[i]] <- colors[i]
    }
    col
}

#' Generic Manhattan plotting function. Requires a data frame with at least 3 columns: 
#' CHR, BP and STAT.
manhattan <- function(d, colors, ylab, chr.sizes=NULL, ymax=NULL, cex.x.axis=1, annotate=NULL, anno.col='green', 
        main=NULL, sub=NULL, threshold=NULL, lty.threshold=1, col.threshold='red', 
        dev.type=c("jpg","pdf"), dev.file=NULL, width=20, height=6, ...) {
    dev.type <- match.arg(dev.type)
    d$pos <- NA
    ticks <- NULL
    lastbase <- 0
    xlim <- NULL
    
    if (!is.null(chr.sizes)) {
        chr.sizes.cs <- cumsum(as.numeric(chr.sizes))
        xlim <- c(0, max(chr.sizes.cs))
    }
    
    if (is.data.frame(d)) {
        if (!("COL" %in% colnames(d))) {
            d$COL <- alternate.colors(d, colors)
        }
    }
    else {
        for (i in 1:length(d)) {
            d[[i]]$COL <- alternate.colors(d[[i]], colors[[i]])
        }
        d <- do.call(rbind, d)
    }
    
    if (is.null(ymax)) {
        ymax <- ceiling(max(c(d$STAT, threshold)))
    }
    
    chrs <- unique(d$CHR)
    chr.idxs <- 1:length(chrs)
    n <- length(chrs)
    if (n == 1) {
        d$pos <- d$BP
        ticks <- floor(length(pos)) / 2 + 1
    } 
    else {
        for (i in chr.idxs) {
            if (i == 1) {
                d[d$CHR==chrs[i], ]$pos <- d[d$CHR==chrs[i], ]$BP
            } 
            else {
                if (!is.null(chr.sizes)) {
                    lastbase <- chr.sizes.cs[i-1]
                }
                else {
                    lastbase <- lastbase + tail(subset(d,CHR==chrs[i-1])$BP, 1)
                }
                d[d$CHR==chrs[i], ]$pos <- d[d$CHR==chrs[i], ]$BP + lastbase
            }
            if (!is.null(chr.sizes)) {
                ticks <- c(ticks, ifelse(i==1,0,chr.sizes.cs[i-1]) + round(chr.sizes[i]/2))    
            }
            else {
                ticks <- c(ticks, d[d$CHR==chrs[i], ]$pos[floor(length(d[d$CHR==chrs[i], ]$pos)/2)+1])
            }
        }
    }
    
    if (is.null(xlim)) {
        xlim <- c(0, max(d$pos))
    }
    
    if (!is.null(dev.file)) {
        if (dev.type == "jpg") {
            jpeg(dev.file, res=300, units="in", width=width, height=height)
        }
        else {
            pdf(dev.file, width=width, height=height)
        }
    }
    
    if (n == 1) {
        with(d, plot(pos, STAT, xlim=xlim, ylim=c(0,ymax), main=main, sub=sub, ylab=ylab, 
            xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }   
    else {
        with(d, plot(pos, STAT, xlim=xlim, ylim=c(0,ymax), main=main, sub=sub, ylab=ylab, 
            xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=as.character(chrs), cex.axis=cex.x.axis)
        for (i in chr.idxs) {
            with(d[d$CHR==chrs[i], ], {
                u <- unique(COL)
                if (length(u) == 1) {
                    points(pos, STAT, pch=21, col=u[1], bg=u[1], ...)
                }
                else {
                    points(pos, STAT, pch=21, col=COL, bg=COL, ...)
                }
            })
        }
    }
    
    if (!is.null(annotate)) {
        d.annotate <- d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, STAT, pch=21, col=anno.col, bg=anno.col, ...)) 
    }
    
    if (!is.null(threshold)) {
        for (th in seq_along(threshold)) {
            abline(h=threshold[th], col=cycle(col.threshold, th), lty=cycle(lty.threshold, th))
        }
    }
    
    if (!is.null(dev.file)) {
        dev.off()
    }
}

#' Base graphics qq plot
qq <- function(pvector, pdf.file=NULL, ...) {
    if (!is.numeric(pvector)) {
        stop("D'oh! P value vector is not numeric.")
    }
    pvector <- pvector[!is.na(pvector) & pvector < 1 & pvector > 0]
    o <- -log10(sort(pvector, decreasing=F))
    e <- -log10(ppoints(length(pvector)))
    
    if (!is.null(pdf.file)) {
        pdf(pdf.file)
    }
    
    plot(e, o, pch=19, cex=1, xlab=expression(Expected~~-log[10](italic(p))), 
        ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
    abline(0, 1, col="red")
    
    if (!is.null(pdf.file)) {
        dev.off()
    }
}

#' Effect size plotting function. Requires a data frame with at least 3 columns: CHR, BP and BETA.
plot.effects <- function(d, chrm=NULL, rng=NULL, colors=c("steelblue", "gray50"), ylab="Effect Size", 
        ylim=NULL, cex.x.axis=1, annotate=NULL, main=NULL, sub=NULL, pdf.file=NULL, width=20, height=6, ...) {
    if (!is.null(chrm)) {
        d <- d[d$CHR==chrm,]
        if (!is.null(rng)) {
            d <- d[d$BP >= rng[1] & d$BP <= rng[2],]
        }
    }

    d$pos <- NA
    ticks <- NULL
    lastbase <- 0
    colors <- rep(colors, max(d$CHR))[1:max(d$CHR)]
    if (is.null(ylim)) {
        ylim <- range(d$BETA)
    }

    n <- length(unique(d$CHR))
    if (n == 1) {
        d$pos <- d$BP
        ticks <- floor(nrow(d)) / 2 + 1
    } 
    else {
        u <- unique(d$CHR)
        for (i in 1:length(u)) {
            ch <- u[i]
            if (i == 1) {
                d[d$CHR==ch, ]$pos <- d[d$CHR==ch, ]$BP
            } 
            else {
                lastbase <- lastbase + tail(subset(d,CHR==u[i-1])$BP, 1)
                d[d$CHR==ch, ]$pos <- d[d$CHR==ch, ]$BP + lastbase
            }
            ticks <- c(ticks, d[d$CHR==ch, ]$pos[floor(length(d[d$CHR==ch, ]$pos)/2)+1])
        }
    }

    if (!is.null(pdf.file)) {
        pdf(pdf.file, width=width, height=height)
    }

    if (n == 1) {
        with(d, plot(pos, BETA, type="l", ylim=ylim, main=main, sub=sub, ylab=ylab, 
            xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }   
    else {
        with(d, plot(pos, BETA, ylim=ylim, main=main, sub=sub, ylab=ylab, 
            xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), cex.axis=cex.x.axis)
        icol <- 1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ], lines(pos, BETA, col=colors[icol], ...))
            icol <- icol + 1
        }
    }

    if (!is.null(pdf.file)) {
        dev.off()
    }
}