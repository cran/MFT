#' MFT.m_est
#'
#' Naive routine for the estimation of the order of serial correlation (m-dependence) in point processes.
#'
#' @param Phi point process, vector of time stamps
#' @param n positive integer, number of life times used in segments for estimation of serial correlation
#' @param maxlag non-negative integer, maximal lag up to which serial correlations are calculated
#' @param alpha numeric, in (0,1), significance level
#' @param plot logical, if TRUE, estimation procedure is plotted
#' 
#'
#' @return
#' \item{m_est}{non-negative integer, estimated order of serial correlation (m-dependence)}
#' 
#' @seealso \code{\link{MFT.rate}, \link{plot.MFT}, \link{summary.MFT}, \link{MFT.variance}, \link{MFT.mean}, \link{MFT.peaks}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Kaue M. Costa, Jochen Roeper and Gaby Schneider (2017).
#' Multi-scale detection of rate changes in spike trains with weak dependencies. Journal of Computational Neuroscience, 42 (2), 187-201.
#' <doi:10.1007/s10827-016-0635-3>
#' 
#' @examples
#' # 1. Independent life times (m=0)
#' set.seed(117)
#' n <- 5000
#' Phi1 <- cumsum(rexp(n,3.5))
#' Phi2 <- cumsum(rexp(n,5))
#' Phi3 <- cumsum(rexp(n,2))
#' Phi  <- c(Phi1[Phi1<=200],Phi2[Phi2>200 & Phi2<400],Phi3[Phi3>400 & Phi3<700])
#' MFT.m_est(Phi)
#' 
#' # 2. Point process simulated according to model
#' # X_i = a_0 X_i + a_1 X_{i-1} + ... a_m X_{i-m}
#' # with life times X_i gamma-distributed, 2 change points and true m = 3.
#' set.seed(210)
#' Tt <- 3000
#' m <- 3
#' a <- c(1,0.5,0.25,0.125)
#' mu <- c(0.5,1,2)/(sum(a))
#' sigmaX <- sqrt(0.225/(sum(a^2)))
#' shape <- mu^2/sigmaX^2; rate <- mu/sigmaX^2
#' len <- 10000
#' # build auxiliary processes
#' X1 <- rgamma(len,rate=rate[1],shape=shape[1]); M1 <- embed(X1,m+1)
#' v1 <- cumsum(as.vector(M1 %*% a)); v1 <- v1[v1<Tt]
#' X2 <- rgamma(len,rate=rate[2],shape=shape[2]); M2 <- embed(X2,m+1)
#' v2 <- cumsum(as.vector(M2 %*% a)); v2 <- v2[v2<Tt]
#' X3 <- rgamma(len,rate=rate[3],shape=shape[3]); M3 <- embed(X3,m+1)
#' v3 <- cumsum(as.vector(M3 %*% a)); v3 <- v3[v3<Tt]
#' # build final point process with cps at 100 and 200
#' Phi <- c(v1[v1<Tt/3],v2[v2>Tt/3 & v2<(2/3)*Tt],v3[v3>(2/3)*Tt])
#' # estimate m
#' MFT.m_est(Phi)
#' 
#' @rdname MFT.m_est
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


MFT.m_est <-
function(Phi,n=200,maxlag=10, alpha=0.05,plot=TRUE)
{
	if (!all(diff(Phi)>=0)) stop("Invalid input process: Strictly increasing sequence of events required")	
	if (length(which(diff(Phi)==0))>0) cat("Warning: Some intervals have length zero")	
	if (!n%%1 == 0 | n <= 0) {
        stop("Invalid choice of number of life times per window: n must be a positive integer")
    }
    if (!maxlag%%1 == 0 | maxlag <= 0) {
        stop("Invalid choice of number of maximal lag: maxlag must be a positive integer")
    }
   if (alpha * (1 - alpha) <= 0) {
        stop("Invalid choice of significance level: alpha must be in (0,1)")
    }
  if (!is.logical(plot)) {
        stop("Invalid plot option: plot must be logical")
    }
 if (length(Phi)<n*2) stop("Invalid choice of n relative to process length: choose n smaller than half the number of events.")
 if (n<maxlag+5) stop("Invalid choice of n: n must be at least maxlag+5")
 
	sercor<-function(Phi,maxlag) 
	{	
		d<-diff(Phi)
		int<-array(dim=c(2,maxlag))
		ser<-cor(d[-1],d[1:(length(d)-1)]) 
		int[,1]<-cor.test(d[-1],d[1:(length(d)-1)])$conf[1:2]
		p<-cor.test(d[-1],d[1:(length(d)-1)])$p.value
		for (i in 2:maxlag)
		{	
			ser[i]<-cor(d[-c(1:i)],d[1:(length(d)-i)])
			int[,i]<-cor.test(d[-c(1:i)],d[1:(length(d)-i)])$conf[1:2]
			p[i]<-cor.test(d[-c(1:i)],d[1:(length(d)-i)])$p.value
		}
		return(list(ser,int,p))
	}

 

	n.wind<-floor(length(Phi)/n)
	sc<-array(NA,dim=c(n.wind, maxlag))
	for (i in 1:n.wind)
	{
		sptr<-Phi[((i-1)*n+1):(i*n)]
		sc[i,]<-sercor(sptr,maxlag)[[1]]
	}
	p<-c()
	for (i in 1:maxlag) p[i]<-wilcox.test(sc[,i])$p.value
	if (length(which(p>alpha))<1) {cat("Warning: all considered serial correlations were significantly different from zero, maxlag may be too small", sep="\n"); m_est<-maxlag}
	if (length(which(p>alpha))>=1) {m_est<-min(which(p>alpha))-1}
	if (plot==TRUE)
	{
	par(mar=c(4.1,4.3,0.5,0.5),cex.axis=1.1,cex.lab=1.3)
  plot(1:maxlag,sc[1,],ylim=c(-1,1),pch=19,cex=0.7,ylab="serial correlation",xlab="lag",col="grey",bty="n")
	for (i in 2:n.wind) points(jitter(1:maxlag),sc[i,],pch=19,cex=0.7,col="grey")
	abline(h=0)
	abline(v=m_est+.5,col="red",lwd=3)
	points(1:maxlag,apply(sc,2,median),pch=18,cex=1.5)
	text(maxlag/2,0.8,bquote(hat(m) == .(m_est)),cex=2)
	}
	return(m_est)
}
