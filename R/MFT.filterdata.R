#' MFT.filterdata
#'
#' Naive routine to remove trend from the data.
#'
#' @param x numeric vector, input sequence of random variables.
#' @param filterwidth postive interger, < length(x)/2, number of data points left and right of the current value that are taken into account for Gaussian smoothing.
#' @param filtersigma numeric, > 0, standard deviation of Gassian kernel.
#
#' @return invisible
#' \item{xfiltered}{filtered data (for filtering the first and last (filterwidth many) data points of the original series cannot be evaluated and are omited)}
#' \item{xraw}{orignal data, but the first and last (filterwidth many) data point are omitted}
#' \item{xtrend}{trend that is removed by filtering. That is xfiltered = xraw - xtrend}
#' \item{x}{orignal data}
#' \item{filterwidth}{number of data points left and right of the current value that are taken into account for Gaussian smoothing}
#' \item{filtersigma}{standard deviation of the Gaussian kernel}
#' 
#' @seealso \code{\link{MFT.peaks}, \link{plot.MFT}, \link{summary.MFT}, \link{MFT.rate}, \link{MFT.variance}, \link{MFT.mean}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Hendrik Backhaus, Albrecht Stroh and Gaby Schneider (2019+). Peak detection in times series 

#' @examples
#' set.seed(0)
#' # Normally distributed sequence with negative trend
#' x <- rnorm(1000,mean=seq(5,0,length.out=1000))
#' MFT.filterdata(x)
#  # Set additional parameters
#' MFT.filterdata(x,filterwidth=200,filtersigma=200)
#' 
#' @rdname MFT.filterdata
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


MFT.filterdata <- function(x,filterwidth=NULL,filtersigma=NULL)
{
  if(is.null(filterwidth)){filterwidth <- length(x)/10}
  if(is.null(filtersigma)){filtersigma <- length(x)/20}
  if(filtersigma < 0){stop("filtersigma must be potive")}
  if(filterwidth >= length(x)/2){stop("Choose filterwidth < half of the length of input data")}
  if(filterwidth %% 1 != 0){stop("filterwidth must be postive integer")}
  
  gausskernel <- dnorm(-filterwidth:filterwidth,sd=filtersigma)
  probweights <- gausskernel / sum(gausskernel)     # Need probability weights, sum = 1
  xtrendNA    <- filter(x,filter=probweights)                 # smoothed data = trend (with NAs at beginning and end due to filtering, number of NAs equals filterwidth at beginning and end)
  xtrend      <- xtrendNA[-which(is.na(xtrendNA))]  # remove NAs of trend (at beginning and end)
  xraw        <- x[-which(is.na(xtrendNA))]         # remove first and last values of x (xraw same length as xtrend)
  xfiltered   <- xraw - xtrend # remove trend from raw data
  
  # Graphical representation 
  par(mfcol=c(2,1),cex=1.3,mar=c(2,5,2,1))
  plot(1:length(x),x,axes=FALSE,xlim=c(1,length(x)),type="l",main="original data",xlab="",ylab="")
  axis(1,at=c(1,length(x))); axis(2)
  plot((1+filterwidth):(length(xfiltered)+filterwidth),xfiltered,axes=FALSE,xlim=c(1,length(x)),type="l",main="filtered data",xlab="",ylab="")
  axis(1,at=c(filterwidth+1,length(x)-filterwidth)); axis(2)
  
  # Output
  list(xfiltered=xfiltered,xraw=xraw,xtrend=xtrend,x=x,filterwidth=filterwidth,filtersigma=filtersigma)
}



