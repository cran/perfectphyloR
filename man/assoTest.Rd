% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/association_tests_utility_functions.R
\name{assoTest}
\alias{assoTest}
\title{Compute the user-provided association statistics}
\usage{
assoTest(Dx, Dy, hapMat, nperm, method, xlab = "", ylab = "",
  main = "")
}
\arguments{
\item{Dx}{A numeric matrix of pairwise distances.}

\item{Dy}{A second numeric matrix of pairwise distances.}

\item{hapMat}{An object of \code{hapMat}.}

\item{nperm}{The number of times to permute the rows and columns of \code{Dy}.}

\item{method}{Association methods. Use "HHG" for HHG test, "dCor" for dCor test, "Mantel" for 
mantel test, "RV" for RV test.}

\item{xlab}{An optional character string for the label on the x-axis(none by default).}

\item{ylab}{An optional character string for the label on the y-axis(none by default).}

\item{main}{An optional character string for title(none by default).}
}
\value{
A list of permutation p-values, corresponding statistics and a plot of association profile 
        over SNV location.
}
\description{
This function computes the user-given association statistics and return the association
profile over a genomic region.
}
\keyword{internal}
