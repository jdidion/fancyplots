library(grid)

my.as.char <- function(obj) {
	if (is.null(obj) || is.na(obj)) {
		''
	}
	else {
		as.character(obj)
	}
}

my.signif <- function(vec, digits=6) {
	vec = signif(vec, digits)
    om = floor(log10(abs(vec)))   # order of magnitude
    dp = digits-om-1              # no. decimal points
    dp[which(dp<0)] = 0           # use zero when dp is negative
    for(i in 1:length(vec)) {
       	if(is.na(vec[i])) { next }
        vec[i] = sprintf(paste("%.",dp[i],"f", sep=""), vec[i])
    }
    return(vec)
}

format.sci <- function(obj, ndig=6) {
	if (is.na(obj)) {
		''
	}
	else {
		format(my.signif(obj, ndig), scientific=TRUE)
	}
}

plot.functions <- list(
	c=function(obj, par=list(), ...) do.call(grid.text, c(my.as.char(obj), par)),
	s=function(obj, par=list(), ...) do.call(grid.text, c(sprintf(as.character(obj), ...), par)),
	n=function(obj, par=list(), ...) do.call(grid.text, c(ifelse(is.na(obj), '', prettyNum(obj, ...)), par)),
	e=function(obj, par=list(), ...) do.call(grid.text, c(format.sci(obj, ...), par))
)

#' Print a table as a figure using grid graphics.
#' 
#' @param data matrix or data frame
#' @param plotters list of plotters to use (will be recycled if length is shorter than
#' `ncol(data)`). Each element is either a name of a pre-defined plot function, or a
#' list with `fn` set to the plotter function name and additional elements being 
#' arguments that will be passed to the plotting function.
#' @param col.widths vector of column widths
#' @param table.border logical; whether to draw a border around the table.
#' @param table.border.gpar list of graphical parameters to use when drawing borders.
#' @param cell.mar size of margin around values in table cells.
#' @param cell.border logical; whether to draw borders around cells.
#' @param cell.border.gpar list of graphical parameters to use when drawing cell borders.
#' @param header.gpar list of graphical parameters to use when drawing column headers.
#' @param default.gpar list of default graphical parameters to use when drawing cell values.
#'
#' @note If you are using knitr, it is much better to use pander to format tables.
plot.table <- function(data, plotters=list('c'), col.widths=NULL, table.border=FALSE, 
        table.border.gpar=NULL,  cell.mar=0.05, cell.border=FALSE, cell.border.gpar=NULL, 
        header.gpar=gpar(fontface='bold'), default.gpar=gpar()) {
	plotters <- lapply(plotters, function(x) {
		if (!is.list(x) && is.character(x)) {
			x <- list(fn=plot.functions[[x]])
		}
		else if (is.list(x) && is.character(x$fn)) {
			x$fn <- plot.functions[[x$fn]]
		}
		x
	})

	if (is.null(col.widths)) {
		col.widths <- unit(rep(1, ncol(data)), "null")
	}
	lo <- grid.layout(nrow(data)+1, ncol(data), widths=col.widths)
	pushViewport(viewport(layout=lo))
	if (table.border) {
		if (is.null(table.border.gpar)) {
			table.border.gpar <- default.gpar
		}
		grid.rect(gp=table.border.gpar)
	}
	for (trow in 1:(nrow(data)+1)) {
		for (col in 1:ncol(data)) {
			if (trow == 1) {
				pushViewport(viewport(layout.pos.row=1, layout.pos.col=col))
				row = 1
			}
			else {
				row <- trow - 1

				ix <- ((col - 1) %% length(plotters)) + 1
				p <- plotters[[ix]]
			
				vp.args <- list(layout.pos.row=trow, layout.pos.col=col)
				if ('vp' %in% names(p)) {
					vp.args <- c(vp.args, p$vp)
				}
				pushViewport(do.call(viewport, vp.args))
			}

			if (cell.border) {
				left   <- ifelse(col == 1, cell.mar, cell.mar / 2)
				bottom <- ifelse(row == nrow(data), cell.mar, cell.mar / 2)
				right  <- ifelse(col == ncol(data), 1 - (left + cell.mar), 1 - (left + cell.mar / 2))
				top    <- ifelse(row == 1, 1 - (bottom + cell.mar), 1 - (bottom + cell.mar / 2))
				if (is.null(cell.border.gpar)) {
					cell.border.gpar <- default.gpar
				}
				grid.rect(left, bottom, right, top, just=c('left','bottom'), gp=cell.border.gpar)
			}
			
			if (trow == 1) {
				grid.text(colnames(data)[col], gp=header.gpar)
			}
			else {
				fn.args <- list(data[row,col])
				if (!('par' %in% names(p)) || is.null(p$par)) {
					p$par <- list(gp=default.gpar)
				}
				else if (!('gp' %in% names(p$par)) || is.null(p$par$gp)) {
					p$par$gp <- default.gpar
				}
				fn.args <- c(fn.args, p[-which(names(p) == 'fn')])
				do.call(p$fn, fn.args)
			}
			popViewport()
		}
	}
}