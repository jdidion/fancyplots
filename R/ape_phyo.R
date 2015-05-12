#' Plot a newick tree using ADE.
plottree <- function(infile, ...) {
    plottree.phylo(read.tree(infile), ...)
}

#' Convenience function to plot an ADE phylo.
plottree.phylo <- function(p, type='phylogram', names.file=NA, pdf.file=NULL, pdf.width=20,
        pdf.height=20, draw.branch.labels=FALSE, convert.branch.labels=NULL, ...) {
    labs = p$tip.label
    if (!is.na(names.file)) {
        labs = read.table(names.file, sep='=', header=FALSE, colClasses='character', row.names=1)[p$tip.label, 1]
        p$tip.label <- labs
    }
    if (!is.null(pdf.file)) {
        pdf(pdf.file, pdf.width, pdf.height)
        on.exit(dev.off())
    }
    plot(p, type, ...)
    if (draw.branch.labels) {
        labs <- p$edge.length
        if (!is.null(convert.branch.labels)) {
            labs <- convert.branch.labels(labs)
        }
        edgelabels(labs, ...)
    }
}
