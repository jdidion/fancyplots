scale_colour_viridis <- function (..., type = "seq", palette = 1) {
    discrete_scale("colour", "viridis", viridis, ...)
}

scale_fill_viridis <- function (..., type = "seq", palette = 1) {
    discrete_scale("fill", "viridis", viridis, ...)
}

grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom", legend.text=element_text(size=6)))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    gridExtra::grid.arrange(
        do.call(gridExtra::arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
}