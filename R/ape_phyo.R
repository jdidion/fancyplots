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

phylo2df <- function(p, make.binary=FALSE) {
    if (is.character(p)) {
        p <- read.tree(p)
    }
    if (make.binary && !is.binary.tree(p)) {
        p <- multi2di(p, random = TRUE)
    }
    nodes <- data.frame(
        node=p$edge[,2],
        parent=p$edge[,1],
        dist=p$edge.length,
        label=paste0("node", p$edge[,2]),
        leaf=FALSE
    )
    ntips <- length(p$tip.label)
    m <- match(1:ntips, p$edge[,2])
    nodes[m, "label"] <- p$tip.label
    nodes[m, "leaf"] <- TRUE
    nodes
}

# This doesn't work the way I expect it to. Need
# to figure out the right metric for weighting nodes.
summarize <- function(p, K=10, err=0.1, out.dir=tempfile(), open.browser=interactive()) {
    if (!is.data.frame(p)) {
        p <- phylo2df(p)
    }
    # Weight leaf nodes as 0 and internal nodes as 1
    p$weight <- ifelse(p$leaf, 0, 1)
    # add a synthetic root node
    root <- setdiff(p$parent, p$node)
    p <- rbind(data.frame(node=root, parent=0, weight=1, label='root'), 
               p[,c("node","parent","weight","label")])
    # compute summary tree
    e <- optimal(
        node=p$node, 
        parent=p$parent, 
        weight=p$weight, 
        label=p$label, 
        K=K, 
        epsilon=err)
    # display summary tree
    json <- prepare.vis(
        tree.list = e$summary.trees,
        labels = e$data[, "label"],
        tree = e$tree,
        legend.width = 120,
        node.width = 150,
        node.height = 12,
        units = "# of descendants",
        print.weights = TRUE,
        legend.color = "lightsteelblue",
        color.level = 2)
    draw.vis(
        json.object=json, 
        out.dir=out.dir, 
        open.browser=open.browser)
    invisible(e)
}
