% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MFT.filterdata.R
\name{MFT.filterdata}
\alias{MFT.filterdata}
\title{MFT.filterdata}
\usage{
MFT.filterdata(x, filterwidth = NULL, filtersigma = NULL)
}
\arguments{
\item{x}{numeric vector, input sequence of random variables.}

\item{filterwidth}{postive interger, < length(x)/2, number of data points left and right of the current value that are taken into account for Gaussian smoothing.}

\item{filtersigma}{numeric, > 0, standard deviation of Gassian kernel.}
}
\value{
invisible
\item{xfiltered}{filtered data (for filtering the first and last (filterwidth many) data points of the original series cannot be evaluated and are omited)}
\item{xraw}{orignal data, but the first and last (filterwidth many) data point are omitted}
\item{xtrend}{trend that is removed by filtering. That is xfiltered = xraw - xtrend}
\item{x}{orignal data}
\item{filterwidth}{number of data points left and right of the current value that are taken into account for Gaussian smoothing}
\item{filtersigma}{standard deviation of the Gaussian kernel}
}
\description{
Naive routine to remove trend from the data.
}
\examples{
set.seed(0)
# Normally distributed sequence with negative trend
x <- rnorm(1000,mean=seq(5,0,length.out=1000))
MFT.filterdata(x)
MFT.filterdata(x,filterwidth=200,filtersigma=200)

}
\references{
Michael Messer, Hendrik Backhaus, Albrecht Stroh and Gaby Schneider (2019+). Peak detection in times series
}
\seealso{
\code{\link{MFT.peaks}, \link{plot.MFT}, \link{summary.MFT}, \link{MFT.rate}, \link{MFT.variance}, \link{MFT.mean}}
}
\author{
Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
}
