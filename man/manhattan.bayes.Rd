% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/gwas.R
\name{manhattan.bayes}
\alias{manhattan.bayes}
\title{Manhattan plot for results from Bayesian analysis. Requires at least 3 column data frame:
CHR, BP, and STAT, corresponding to chromosome number, genomic coordinate, and Bayes factor.
Missing values (where regression failed to converge, for instance) should be recoded NA before
reading into R. The plot alternates between two colors.}
\usage{
manhattan.bayes(d, mapping = NULL, colors = c("steelblue", "gray50"),
  chrms = 1:19, sub = "Bayesian Association", ...)
}
\description{
Manhattan plot for results from Bayesian analysis. Requires at least 3 column data frame:
CHR, BP, and STAT, corresponding to chromosome number, genomic coordinate, and Bayes factor.
Missing values (where regression failed to converge, for instance) should be recoded NA before
reading into R. The plot alternates between two colors.
}

