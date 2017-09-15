#' MFT.m_est
#'
#' Naive routine for the estimation of the order of serial correlation (m-dependence) in point processes.
#'
#' @param Phi spike train, vector of time stamps
#' @param n positive integer, number of life times used in segments for estimation of serial correlation
#' @param maxlag non-negative integer, maximal lag up to which serial correlations are calculated
#' @param plot logical, if TRUE, estimation procedure is plotted
#'
#' @return
#' \item{m_est}{non-negative integer, estimated order of serial correlation (m-dependence)}
#' 
#' @seealso \code{\link{MFT.rate}, \link{MFT.variance}, \link{MFT.mean}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Kaue M. Costa, Jochen Roeper and Gaby Schneider (2017).
#' Multi-scale detection of rate changes in spike trains with weak dependencies. Journal of Computational Neuroscience, 42 (2), 187-201.
#' <doi:10.1007/s10827-016-0635-3>
#' 
#' @examples
#' # 1. Homogeneuous Poisson process (m = 0)
#' Phi <- cumsum(rexp(2000,rate=5)) 
#' MFT.m_est(Phi)
#' 
#' # 2. Point process simulated according to model
#' # X_i = a_0 X_i + a_1 X_{i-1} + ... a_m X_{i-m}
#' # with life times X_i gamma-distributed and true m = 3.
#' m   <- 3
#' Tt  <- 300
#' # generate coefficients a_i
#' c   <- 0.5
#' aa  <- c(1,c); for (i in 2:11) aa[i] <- c*aa[i-1]
#' a   <- aa[1:(m+1)]
#' # build auxiliary processes
#' muX1 <- 0.05/(sum(a)); muX2 <- 0.1/(sum(a)); muX3 <- 0.2/(sum(a))
#' sigmaX <- sqrt(0.0225/(sum(a^2)))
#' shape1 <- muX1^2/sigmaX^2;rate1 <- muX1/sigmaX^2
#' shape2 <- muX2^2/sigmaX^2;rate2 <- muX2/sigmaX^2
#' shape3 <- muX3^2/sigmaX^2;rate3 <- muX3/sigmaX^2
#' len  <- 10000
#' X1   <- rgamma(10000,rate=rate1,shape=shape1); M1 <- embed(X1,m+1)
#' Phi1 <- cumsum(as.vector(M1 %*% a)); Phi1 <- Phi1[Phi1<Tt]
#' X2   <- rgamma(10000,rate=rate2,shape=shape2); M2 <- embed(X2,m+1)
#' Phi2 <- cumsum(as.vector(M2 %*% a)); Phi2 <- Phi2[Phi2<Tt]
#' X3   <- rgamma(10000,rate=rate3,shape=shape3); M3 <- embed(X3,m+1)
#' Phi3 <- cumsum(as.vector(M3 %*% a)); Phi3 <- Phi3[Phi3<Tt]
#' # build final point process
#' Phi  <- sort(c(Phi1[Phi1<Tt/3],Phi2[Phi2>Tt/3 & Phi2<(2/3)*Tt],Phi3[Phi3>(2/3)*Tt]))
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
function(Phi,n=200,maxlag=10,plot=TRUE)
{
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
		sptr<-Phi[((i-1)*n+1):(i*n+1)]
		sc[i,]<-sercor(sptr,maxlag)[[1]]
	}
	p<-c()
	for (i in 1:maxlag) p[i]<-wilcox.test(sc[,i])$p.value
	m_est<-min(which(p>0.05))-1
	if (plot==TRUE)
	{
plot(1:maxlag,sc[1,],ylim=c(-1,1),pch=19,cex=0.7,ylab="serial correlation",xlab="lag",col="grey",bty="n")
	for (i in 2:n.wind) points(jitter(1:maxlag),sc[i,],pch=19,cex=0.7,col="grey")
	abline(h=0)
	abline(v=m_est+.5,col="red",lwd=3)
	points(1:maxlag,apply(sc,2,median),pch=18,cex=1.5)
	text(m_est,0,expression(hat(m)),cex=3)
	}
	return(m_est)
}
