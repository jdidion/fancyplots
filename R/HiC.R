# Load data from homer

read.tag.info <- function(f) {
    tab <- read_tsv(f, col_names=T)
    N <- length(tab[[1]]) - 1
    chrms <- data.frame(numTags=tab[[3]][10:N], row.names=tab[[1]][10:N], stringsAsFactors=F)
    metrics <- do.call(rbind, sapply(tab[[1]][2:9], str_split, '='))
    rownames(metrics) <- metrics[,1]
    metrics <- metrics[,-1]
    metrics <- c(numTags=tab[[3]][1], metrics)
    metrics <- data.frame(value=metrics, row.names=names(metrics), stringsAsFactors=F)
    list(chrms=chrms, metrics=metrics)
}

read.multiple.tag.info <- function(d) {
    libs <- list.dirs(d, recursive=FALSE)
    lib.names <- sapply(libs, basename)
    info.files <- sapply(libs, file.path, 'tagInfo.txt')
    dat <- lapply(info.files, read.tag.info)
    chrms <- lapply(dat, "[[", 'chrms')
    chrm.names <- unique(unlist(lapply(chrms, rownames)))
    chrms <- do.call(cbind, lapply(chrms, function(x) x[chrm.names,]))
    chrms[is.na(chrms)] <- 0
    rownames(chrms) <- chrm.names
    colnames(chrms) <- lib.names
    metrics <- do.call(cbind, lapply(dat, "[[", 'metrics'))
    colnames(metrics) <- lib.names
    list(chrms=chrms, metrics=metrics)
}

# Plot data from homer

plot.local.distribution <- function(f, main=NULL, span=0.01) {
    tab <- read.table(f, sep="\t", header=T)
    w <- tab[,1] >= -1000 & tab[,1] <= 1000
    ss <- supsmu(tab[w,1], tab[w,3], span=span)
    plot(ss, type='l', col="blue", lwd=2,
        xlab="Distance from 5' of 1st Read", ylab="Counds per bp",
        main=main)
    lines(supsmu(tab[w,1], tab[w,2], span=span), type='l', col='orange', lwd=2)
    y <- max(ss$y)
    x <- ss$x[ss$y == y]
    segments(x,0,x,y,col='red')
    text(x,round(y/50),paste0("Median frag len: ", round(x)),pos=4)
}

plot.freq.distribution <- function(f, ...) {
    title <- read.table(f, sep="\t", header=F, nrows = 1, stringsAsFactors=F)[1,2]
    tab <- read.table(f, sep="\t", header=T, nrows = 300000)
    w <- tab[,2]>0
    m <- min(tab[w, 2])
    tab[1, 1] <- 1
    tab[!w, 2] <- 10^(floor(log10(m)) - 1)
    plot(smooth.spline(log10(tab[,1]), log10(tab[,2]), ...), type='l', xaxt='n', frame.plot = F, 
        xlab="Log10(Distance between regions)", ylab="Log10(Fraction of reads) (1kb intervals)")
    axis(3)
    title(title, line=3)
}