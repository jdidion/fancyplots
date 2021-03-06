% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/chromoplot.R
\name{plot.genome.data}
\alias{plot.genome.data}
\title{Plot data on one or more chromosomes as a function of genomic position. Mapping must map columns
of data to the inputs required by plotting functions (x1, x2, y1, etc). type must one of the
plotting functions in .plot.functions or a function that accepts two arguments:
a data frame (subset of data for a single chromosome) and a gpar.
NOTE: It is strongly recommended to use a non-pdf device type (specified using the device
parameter) when data is a list due to the graphics antialiasing that most viewers do by default.}
\usage{
plot.genome.data(data, mapping = NULL, types = c("point"),
  karyotype = NULL, chromosomes = NULL, ylim = NULL, gps = list(gpar()),
  legend = NULL, legend.vp = NULL, device = "pdf", dev.file = NULL,
  dev.width = 7, dev.height = 7, yscale = FALSE, yscale.label = TRUE,
  sample.borders = FALSE, ...)
}
\description{
Plot data on one or more chromosomes as a function of genomic position. Mapping must map columns
of data to the inputs required by plotting functions (x1, x2, y1, etc). type must one of the
plotting functions in .plot.functions or a function that accepts two arguments:
a data frame (subset of data for a single chromosome) and a gpar.
NOTE: It is strongly recommended to use a non-pdf device type (specified using the device
parameter) when data is a list due to the graphics antialiasing that most viewers do by default.
}

