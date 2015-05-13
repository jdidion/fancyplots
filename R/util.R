#' Convert a 2D size from one unit to another.
convert.size <- function(dim, from="inches", to="user") {
    x <- grconvertX(c(0,1), from, to)
    y <- grconvertY(c(0,1), from, to)
    c(dim[1]  * (x[2] - x[1]), dim[2] * (y[2] - y[1]))
}

#' Generate a color ramp function for a given pallate and intensity range.
make.ramp <- function(pal, start=1, end=100) {
    g <- (end - start + 1) / length(pal)
    return(function(n) pal[ceiling((n - start + 1) / g)])
}

#' Create an arbitrary-sized palette from an RColorBrewer palette.
#'
#' @param name RColorBrewer palette name
#' @param cb.colors number of colors to use from the RColorBrewer palette
#' @param pal.colors number of colors to generate from the RColorBrewer palette
make.cb.pal <- function(name, cb.colors=8, pal.colors=255) {
    ramp <- colorRampPalette(brewer.pal(cb.colors, name))
    ramp(pal.colors)
}

#' Modifiy colors (names or rgb vectors) with an alpha value.
rgb.alpha <- function(r, alpha=1.0) {
    if (is.character(r)) {
        r <- col2rgb(r)
    }
    if (is.matrix(r)) {
        if (length(alpha) == 1) {
            apply(r, 2, rgb.alpha, alpha)
        }
        else {
            sapply(1:ncol(r), function(i) rgb.alpha(r[,i], alpha[i]))
        }
    }
    else {
        rgb(r[1], r[2], r[3], alpha * 255, maxColorValue=255)
    }
}

