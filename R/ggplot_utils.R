GeomRugAlt <- proto(Geom, {
    draw <- function(., data, scales, coordinates, ...) {
        rugs <- list()
        data <- coordinates$transform(data, scales)
        if (!is.null(data$x)) {
            rugs$x <- with(data, segmentsGrob(x0 = unit(x, 
              "native"), x1 = unit(x, "native"), y0 = unit(0.955, 
              "npc"), y1 = unit(1, "npc"), gp = gpar(col = alpha(colour, 
              alpha), lty = linetype, lwd = size * .pt)))
        }
        if (!is.null(data$y)) {
            rugs$y <- with(data, segmentsGrob(y0 = unit(y, 
              "native"), y1 = unit(y, "native"), x0 = unit(0.955, 
              "npc"), x1 = unit(1, "npc"), gp = gpar(col = alpha(colour, 
              alpha), lty = linetype, lwd = size * .pt)))
        }
        gTree(children = do.call("gList", rugs))
    }
    objname <- "rug_alt"
    desc <- "Marginal rug plots"

    default_stat <- function(.) StatIdentity

    default_aes <- function(.) aes(colour = "black", size = 0.5, 
        linetype = 1, alpha = 1)

    guide_geom <- function(.) "path"

    examples <- function(.) {
        p <- ggplot(mtcars, aes(x = wt, y = mpg))
        p + geom_point()
        p + geom_point() + geom_rug_alt()
        p + geom_point() + geom_rug_alt(position = "jitter")
    }
})
geom_rug_alt <- GeomRugAlt$build_accessor()
