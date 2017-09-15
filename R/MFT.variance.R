
#' MFT.variance
#'
#' The multiple filter test for variance change detection in point processes on the line. 
#'
#' @param Phi numeric vector of increasing events, input point process
#' @param rcp vector, rate CPs of Phi (if MFT for the rates is used: as CP[,1]), default: constant rate
#' @param S	numeric, start of time interval, default: Smallest multiple of d that lies beyond min(Phi)
#' @param E numeric, end of time interval, default: Smallest multiple of d that lies beyond max(Phi), needs E > S.
#' @param autoset.d_H	logical, automatic choice of window size H and step size d
#' @param d numeric, > 0, step size delta at which processes are evaluated. d is automatically set if autoset.d_H = TRUE
#' @param H vector, window set H, all elements must be increasing ordered multiples of d, the smallest element must be >= d and the largest =< (T/2). H is automatically set if autoset.d_H = TRUE
#' @param alpha numeric, in (0,1), significance level
#' @param sim integer, > 0, No of simulations of limit process (for approximation of Q), default = 10000
#' @param method either "asymptotic", or "fixed", defines how threshold Q is derived, default: "asymptotic". If "asymptotic": Q is derived by simulation of limit process L (Brownian motion); possible set number of simulations (sim). If "fixed": Q may be set automatically (Q)
#' @param Q	numeric, rejection threshold, default: Q is simulated according to sim and alpha.
#' @param perform.CPD logical, if TRUE change point detection algorithm is performed
#' @param plot.CPD logical, if TRUE CPD-scenario is plotted. Only active if perform.CPD == TRUE
#' @param plot.var logical, should the variance histogram be plotted? Only possible, if plot.CPD=TRUE
#' @param print.output logical, if TRUE results are printed to the console	
#' @param col "gray" or vector of colors of length(H). Colors for (R_ht) plot, default: NULL -> rainbow colors from blue to red. 
#' @param ylab1 character, ylab for 1. graphic
#' @param ylab2 character, ylab for 2. graphic 
#' @param cex.legend 	numeric, size of annotations in plot
#' @param cex.diamonds numeric, size of diamonds that indicate change points
#' @param wid integer,>0, width of bars in variance histogram	
#' @param main 	logical, indicates if title and subtitle are plotted
#' @param plot.Q 	logical, indicates if rejection threshold Q is plotted
#' @param plot.M	logical, indicates if test statistic M is plotted
#' @param plot.h 	logical, indicates if a legend for the window set H is plotted
#'
#' @return invisible
#' \item{M}{test statistic}
#' \item{varQ}{rejection threshold}
#' \item{sim}{number of simulations of the limit process (approximation of Q)}
#' \item{CP}{set of change points estmated by the multiple filter algorithm, increasingly ordered in time}
#' \item{var}{estimated variances between adjacent change points}
#' \item{S}{start of time interval}
#' \item{E}{end of time interval}
#' \item{H}{window set}
#' \item{d}{step size delta at which processes were evaluated}
#' \item{alpha}{significance level}
#' 
#' @examples 
#' # Rate and variance change detection in Gamma process 
#' # (rate CPs at t=30 and 37.5, variance CPs at t=37.5 and 52.5) 
#' set.seed(51)
#' mu <- 0.03; sigma <- 0.01
#' p1 <- mu^2/sigma^2; lambda1 <- mu/sigma^2
#' p2 <- (mu*0.5)^2/sigma^2; lambda2 <- (mu*0.5)/sigma^2
#' p3 <- mu^2/(sigma*1.5)^2; lambda3 <- mu/(sigma*1.5)^2
#' p4 <- mu^2/(sigma*0.5)^2; lambda4 <- mu/(sigma*0.5)^2
#' Phi<- cumsum(c(rgamma(1000,p1,lambda1),rgamma(500,p2,lambda2),
#' rgamma(500,p3,lambda3),rgamma(300,p4,lambda4)))
#' # rcp  <- MFT.rate(Phi)$CP[,1] # MFT for the rates
#' rcp <- c(30,37.5) # but here we assume known rate CPs
#' MFT.variance(Phi,rcp=rcp) # MFT for the variances
#' 
#' @seealso \code{\link{MFT.rate}, \link{MFT.mean}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Stefan Albert, Michael Messer, Julia Schiemann, Jochen Roeper and Gaby Schneider (2017) 
#' Multi-scale detection of variance changes in renewal processes in the presence of rate change points.
#' Journal of Time Series Analysis, <doi:10.1111/jtsa.12254>
#' 
#' @rdname MFT.variance
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end

MFT.variance <-
function(Phi,rcp=NULL,autoset.d_H=TRUE,S=NULL,E=NULL,d=NULL,H=NULL,alpha=0.05,sim=10000,method="asymptotic",Q=NA,
                       perform.CPD=TRUE,plot.CPD=TRUE,plot.var=TRUE,print.output=TRUE,col=NULL,ylab1=expression(abs(G[list(h,t)])),ylab2=expression(widehat(sigma)^2),cex.legend=1.2,cex.diamonds=1.4,wid=NULL,
                       main=TRUE,plot.Q=TRUE,plot.M=TRUE,plot.h=TRUE)
{
  ###
  ### Set Parameters
  ###
  
  if(is.null(S)  & is.null(E))	{S <- min(Phi); E <- max(Phi)}
  if(is.null(S)  & !is.null(E))	{S <- min(Phi)}
  if(!is.null(S) & is.null(E))	{E <- max(Phi)}	
  if(E-S <= 0){stop("Invalid choice of S and E: Need S < E.")}
  
  if( alpha*(1-alpha) <= 0){stop("Invalid choice of significance level: alpha must be in (0,1)")}
  if(!sim%%1==0 | sim <= 0 | is.logical(sim)) {stop("Invalid choice of number of simulations: sim must be a positive integer")}
  if ((sim<10000)&(method=="asymptotic")) {warning("Number of simulations for derivation of Q might be too low")}
  step <- d; result2 <- NA
  
  if (min(diff(Phi))<=0) {stop("Not all lifetimes are positive")}
  if(autoset.d_H){   # Automatical setting of step and window sizes
    if (length(Phi[Phi>S & Phi<=E])<=300) stop("Not enough events! Need more than 300")
    H   <- signif(150 / (length(Phi[Phi>S & Phi<=E])/(E-S)),digits=1)  	# rounds to first significant
    step   <- (H / 20)			 									# step size	
    S   <- floor(S/step)*step									# floor S to next d
    Tt  <- ceiling((E-S)/step)*step 					# ceil Tt to next d
    E   <- S + Tt														  # set end
    if(Tt/5 > H){H <- seq(H,Tt/5,H)}			 		# Window choice	
  }#end-if-set.H_and_d
  
  if(!autoset.d_H){
    if(is.null(d) | is.null(H)){stop("if autoset.d_H is FALSE, the step size d and vector of window sizes H must be set")}
    if(d<=0){stop("Invalid choice of step size d: Need d > 0")}	
    if ((2*max(H))>(E-S)) stop("Window too large!")
    S   <- floor(S/step)*step												# floor S to next d
    Tt  <- ceiling((E-S)/step)*step 				      	# ceil Tt to next d
    E   <- S + Tt		
    if(!all(c(diff(c(d,H,(Tt/2))))>0) | !all(H%%d == 0)){stop("Invalid choice of window set: H needs to be an increasing ordered vector of multiples of d, with min(H) > d and max(H) <= Tt/2")}
  }#end-if-!set.H_and_d 
  
  
  
  # Using only process between E and S, but save original process data
  Phi <- Phi[Phi>=S & Phi<=E]-S
  if (length(rcp)>0) {if ((min(rcp)<=S)|(max(rcp)>=E)) {warning("Rate change point(s) outside of the process")}}
  rcp <- rcp[rcp>=S & rcp<=E]
  S.old <- S; E.old <- E
  S   <- floor((S-S.old)/step)*step												# floor S to next d
  Tt  <- ceiling((E-S-S.old)/step)*step 				      	# ceil Tt to next d
  E   <- S + Tt	
  
  if (length(rcp)==0) {rcp <- 0; est.mean <- mean(diff(Phi)); cat("","\n"); print("MFT for the variances assumes a constant rate")} # mean rate is derived
  else {
    if (anyDuplicated(rcp)>0) {stop("Duplicate in rcp-vector")}
    rcp <- rcp-S.old
    est.mean <- 1/(hist(Phi,breaks=c(S,rcp,E),plot=FALSE)$counts / diff(c(S,rcp,E)))
    rcp <- c(0,rcp)
  }
  l <- length(rcp)
  
  
  if(! method %in% c("aymptotic","fixed")){"Invalid choice of method: method must be 'asymptotic' or 'fixed'"}
  if(method == "fixed" & !is.numeric(Q)){stop("Invalid choice of Q: Q must be positive real number")}
  if(method == "fixed" & !(is.numeric(Q) & Q>0)){cat("Warning: non-positive threshold might be inappropriate. Possibly choose Q > 0",sep="\n") }
  if(method != "fixed" & !is.na(Q)){cat("Warning: Q is derived by simulation. In order to set Q manually, choose: method = 'fixed'",sep="\n")}
  
  
  ###
  ### Simulation of Threshold Q, depending on T, H and alpha
  ###
  
  sim.Q <- function(sim=sim,Tt=Tt,H=H,alpha=alpha,d=d){
    
    maxh <- function(h=h,W=W,Tt=Tt,d=d){ # Given Brownian motion max(L_ht) is calculated for fixed h in H
      Wt_plus <- W[(2*h/d +1):(round(Tt/d,0)+1)] 
      Wt		<- W[(h/d +1):(round((Tt-h)/d,0) +1)] 
      Wt_minus<- W[1:(round((Tt-2*h)/d,0) +1)]	
      return(max((1/sqrt(2*h))*abs(Wt_plus - 2*Wt + Wt_minus)))
    }#end-maxh
    
    sim.maxH <- function(Tt=Tt,H=H,d=d){ # Simulates Brownian motion and calculates max(L_ht) for all h in H
      W  <- c(0,cumsum(rnorm(round(Tt/d,0),sd=sqrt(d))))
      return(vapply(as.matrix(H),FUN=maxh,W=W,Tt=Tt,d=d,numeric(length(1))))
    }#end-maxH	
    
      re 		<- matrix(replicate(sim,sim.maxH(Tt=Tt,H=H,d=d)),nrow=sim,byrow=TRUE)	
      maxmax 	<- apply(t(re), MARGIN=2,max)
      Q 		<- quantile(maxmax,1-alpha)
      return(Q)
   
    
  }#end-sim.parameter	
  
  #if (!is.na(Q)) {warning("Input Q is used. Q is not simulated.")}
  if (method=="asymptotic") Q <- sim.Q(Tt=Tt,H=H,alpha=alpha,d=step,sim=sim) # Simulate threshold
  
  if (l>0)  tc <- c(1,rep(0,l-1)) # Saving rate-CPs
  else tc <- c()
  
  ltc <- length(tc)
  if (ltc>1) 
  {
    for (i in 2:(ltc))
    {
      tc[i] <- which(Phi>rcp[i])[1]
    } #end-for
  } #end-if
  muz <- est.mean[1]
  # Creating mean-valued vector
  if ( ltc>1) 
  {
    for (i in 1:(ltc-1))
    {
      muz <- c(muz,rep(est.mean[i],tc[i+1]-tc[i]-1),NA)  
    } #end-for
  } #end-if
  
  muz <- c(muz,rep(est.mean[ltc],length(Phi)-tc[ltc])); muz<-muz[-1]
  hilfsm <- matrix(NA,length(H),30) #Auxiliary matrix for CPs
  if(is.null(col)){colors <- rainbow(length(H),start=2/3,end=0)} else{if(length(col) == 1 & col[1] == "gray"){colors <- gray.colors(length(H),0.8,0)}
    else{colors <- col; if(length(col)!=length(H)){stop("Argument col must be either NULL, ''gray'', or a vector of colors of length of H")}}
  } # Vector of colors for different windows  
  GhtAll <- matrix(NA,length(H),1.1*ceiling((E-S)/(step))); G<-GhtAll
  
  bar <- (mean(diff(Phi)))*100
  if (is.null(wid)) wid <- floor(max(Phi))/(floor(floor(max(Phi))/bar)) 
  if (wid>((E-S)/2)) {stop("width of variance bars too large")}
  
  ###
  ### Function to derive Ght for fixed h
  ###
  
  deriveGht <- function(t,Phi,h,tc,est.mean) {
    nl <- length(Phi[t-h<Phi & Phi<t])     # Counting number of events
    nr <- length(Phi[t<Phi & Phi<t+h])
    phi.t <- which(Phi>=t) [1]
    el <- which(Phi>t-h) [1]
    if (is.na(el)==TRUE) el <- length(Phi)
    er <- max(which(Phi<t+h))
    quaddevl <- 0; quaddevr <- 0;  lifet <- diff(Phi)     # Calculating lifetimes.
    quaddevl <- ((lifet[el:(nl+el-1)]-est.mean[el:(nl+el-1)])^2); quaddevl <- na.omit(quaddevl)
    quaddevr <- ((lifet[(er-nr):(er-1)]-est.mean[(er-nr):(er-1)])^2); quaddevr <- na.omit(quaddevr)
    nl <- nl-length(na.omit(unique(est.mean[el:(nl+el-1)]))); nr <- nr-length(na.omit(unique(est.mean[(er-nr):(er-1)])))
    # Calculating estimated variances
    if (nl<=1) sigmal <- 0
    else     sigmal <- sum(quaddevl)/(nl-1)
    if (nr<=1) sigmar <- 0
    else     sigmar <- sum(quaddevr)/(nr-1)  
    # Calculating norming factor.
    if (nr<=1) {corr <- 0}   else {corr <- sum((quaddevr-sigmar)^2)/(nr-1)}
    if (nl<=1) {corl <- 0}   else  {corl <- sum((quaddevl-sigmal)^2)/(nl-1)}
    if (nl<=1) {hmul <- 0} else {hmul <- h/mean(diff(Phi[(el):(el+nl-1)]))}
    if (nr<=1) {hmur <- 0} else {hmur <- h/mean(diff(Phi[(er-nr):(er-1)]))}
    corb <- sqrt(sum(corl)/(hmul-1)+sum(corr)/(hmur-1))
    # Calculate test statistic  
    if (min(c(hmul,hmur))<=1) dht <- 0
    else dht <- (sigmal-sigmar)/corb
    
    return(dht)
  }#end-deriveGht
  
  ###
  ### Function to derive Ght for fixed H
  ###
  
  ght.mat <- function(h,Phi,S,E,d,tc,est.mean)
  {
    ght <- vapply(seq(S+h,E-h,by=d),deriveGht,Phi=Phi,h=h,tc=tc,est.mean=est.mean,numeric(1)) 
  }#end-ght.mat
  
  ###
  ### Function to estimate variance in a window of size h
  ###
  
  var.window <- function (t,Phi,h,tc,est.mean)
  {
    # Calculation of lifetimes
    lifet <- diff(Phi); h <- h/2; n <- length(Phi[t-h<Phi & Phi<t+h])
    el <- which(Phi>t-h) [1]
    if (is.na(el)==TRUE) el <- length(Phi)-1  
    er <- max(which(Phi<t+h))
    if (er<2) {er <- length(Phi)}
    if (is.na(er)==TRUE) {er <- length(Phi)}
    quaddev <- (lifet[el:(er-2)]-est.mean[el:(er-2)])^2   #Estimating variances
    quaddev <- na.omit(quaddev)
    n <- n-length(na.omit(unique(est.mean[el:(er-2)])))
    resu <- sum(quaddev)/(n-2)
    return(resu)
  }#end-var.window
  
  ###
  ### Plotting G_ht-processes for every h and t
  ###
  
  plot.ght <- function(ght,H,step,colors,Q,cex.legend,plot.Q=TRUE,plot.M=TRUE,plot.h=TRUE,ylab1) {
    max.ght <- max(abs(ght),na.rm=TRUE)
    xm <- c(); x_max <- c()
    par(mar=c(0.5,4,4,0.5),cex=1)
    
    for (i in 1:length(H))
    {
      Dz <- ght[1,]; Dz <- Dz[is.na(Dz)==FALSE]
      border <- 1.2*max(Q,max.ght) ;
      # Save max values for each i to determine M    
      g <- na.omit(abs(ght[i,])); xm[i] <- which(g==max(g))[1]; x_max[i] <- abs(g[xm[i]])
      if (i==1) #For i=1 create new plot, otherwise add lines
      { 
        plot(abs(ght[i,]), main="", xlab="", ylab="", type="l", ylim=c(-0.6,border),xaxt="n",
             xlim=c(0,(length(Dz)+2*H[1]/step)), cex.axis=cex.legend, font.main=1, 
             cex.lab=cex.legend, cex.main=1, bty="n", col=colors[1],yaxt="n")
        lang <- length(ght[i,])
        mtext(ylab1,side=2,padj=-2.6,cex=cex.legend)
        if (plot.h) {legend(x="topright",inset=c(0,0.02),legend=as.character(round(H,0)),pch=19,col=colors,cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),bty="n",title=expression(paste(h[i], " =")),pt.cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),title.adj=0.2,xpd=TRUE)}
        if (plot.Q) {
          if (main) mtext(bquote(Q %~~% .(round(Q,2)) ),adj=1,padj=-0.5,cex=cex.legend)
          lines(c(0,length(Dz)+2*H[1]/step),c(Q,Q),col="black", lty=2, lwd=1); text(S+15,Q,"Q",pos=3,cex=1.0*cex.legend)
        } #end-if
      } #end-if
      if (i>1) lines(abs(ght[i,]),col=colors[i])
    } #end-for
    if (plot.M) {
      if (main) mtext(bquote(M %~~% .(round(max(abs(ght),na.rm=TRUE),2)) ),adj=0,padj=-0.5, cex= cex.legend); 
      posmax <- which(x_max==max(x_max)) #Print M  
      points(xm[posmax]+H[posmax]/step,x_max[posmax],pch=4,cex=cex.legend,lwd=1)
      text(xm[posmax]+(1.3*H[posmax])/step,x_max[posmax],paste("M"),cex=1.0*cex.legend)
    }
    axis(2,cex.axis=cex.legend) 
    
    
  }#end-plot.ght
  
  ###
  ### Plotting of variances
  ###
  
  plot.variances <- function(Phi,S,E,S.old,wid,tc,muz,l1,lg,varianz,true.var,xl,h,step,colors,hvalue,cex.legend,ylab2)
  {    
    zeitpl <- seq(S+(wid/2),E-(wid/2),by=wid)
    varf <- sapply(zeitpl,var.window,Phi=Phi,h=wid,tc=tc,est.mean=muz)
    lzpl <- length(zeitpl)
    par(mar=c(2,4,0.5,0.5),cex=1)
    
    yl <- 1.2*max(varf[is.na(varf)==FALSE & is.finite(varf)==TRUE])
    plot(1,main="",xlab="",ylim=c(-0, yl),
         ylab="",type="n", xaxt="n", cex.axis=cex.legend, cex.lab=1.0, cex.main=1.0,bty="n",
         xlim=c(0,0+(lzpl)*xl),yaxt="n")
    myTicks <- axTicks(2)
    axis(2, at = myTicks, labels = formatC(myTicks, format = 'e',digits=1),cex.axis=cex.legend)
    
    for (m in 1:lzpl)
    {
      lines(c(m-1,m-1)*xl,c(0,varf[m])); lines(c(m-1,m)*xl,c(varf[m],varf[m])) 
      lines(c(m,m)*xl,c(varf[m],0)) 
    } #end-for
    lines(c(0,lzpl)*xl,c(0,0))
    mtext(ylab2,side=2,padj=-2.6,cex=cex.legend)
    const.var <- mean((diff(Phi)-muz)^2,na.rm=TRUE)
    if (l1==0) lines(c(0,length(varf[is.na(varf)==FALSE])*xl),rep(const.var,2),lwd=2, col="red") 
    else 
      for (k in 1:l1) 
      {
        lg[k] <- ((varianz[3*k-2])*xl*lzpl)/E
        if (varianz[3*k-1]<=varianz[3*k]) 
        { 
          arrows(lg[k],varianz[3*k-1],lg[k],varianz[3*k],col="red", lwd=2, length=0)
          arrows(lg[k],1.8*max(varf),lg[k],1.1*varianz[3*k],col=colors[which(h==hvalue[k])],xpd=TRUE,length=0.1,lwd=2)
        }
        else 
        {
          arrows(lg[k],varianz[3*k-1],lg[k],varianz[3*k],col="red", lwd=2, length=0)
          arrows(lg[k],1.8*max(varf),lg[k],1.1*varianz[3*k-1],col=colors[which(h==hvalue[k])],xpd=TRUE,length=0.1,lwd=2)
        }
        
        if (k==1) 
        {
          lines(c(0.0*xl,lg[k]),c(rep(varianz[3*k-1],2)),lwd=2,col="red")
        }
        else
        {
          lines(c(lg[k-1],lg[k]),c(rep(varianz[3*k-1],2)),lwd=2,col="red")
        }
      } #end-for
    if (l1>0) lines(c((varianz[3*l1-2])*xl*lzpl/E,xl*lzpl),
                    c(rep(varianz[3*l1],2)),lwd=2,col="red") 
    axis(1,S+c(0.0*xl,lg,length(varf[is.na(varf)==FALSE])*xl),c(round(S+S.old,2),round(true.var+S.old,2),round(E+S.old,2)),cex.axis=cex.legend)
  }#end-plot.variances
  
  ###
  ### Change-Point-Detection-Function
  ###
  
  cpa <- function(ght,Q,H,step)
  {
    # Application of SWD for every h 
    for (i in 1:length(H))
    {    
      D <- ght[i,]; D[is.na(D)] <- 0    # Comparison with threshold Q
      cp <- which(abs(D)>Q); k <- 1; varianz <- c(rep(0,50)); D1 <- D; mac <- 0
      while (length(cp)>0) {
        ma <- which(abs(D1)==max(abs(D1))) [1]
        if (mac==ma) break
        mac <- ma
        varianz[k] <- ma
        k <- k+1
        D1[max(0,(ma-H[i]/step)):(min(length(D1),(ma+H[i]/step))+1)] <- 0         # Setting neighbourhood points to zero
        cp <- which(abs(D1)>Q)
      } #end-while
      
      hilfsv <- varianz[varianz!=0]; hilfsv <- sort(hilfsv)             # Saving variance cps on auxiliary vector
      if (length(hilfsv)>0) hilfsm[i,1:length(hilfsv)] <- hilfsv      # If there is CP it is marked in the matrix
    } #end-for
    
    # Decision over final CPs is met
    cpcand <- 0 ; cpf <- rep(NA,50); hvalue <- rep(NA,50); hilfcpf <- sort(hilfsm[1,]); acp <- length(hilfcpf)
    # CPs of smallest windows are saved.
    if (acp>0) {cpf[1:acp] <- hilfcpf 
    hvalue[1:acp] <- H[1]} #end-if
    
    # CPs of large windows are checked and saved if necessary.
    if (length(H)>1)
      for (i in 2:length(H))
      {
        hilfcpf <- sort(hilfsm[i,])
        for (j in 1:length(hilfcpf))
        {
          cpcand <- hilfcpf[j]
          cpfwna <- c(0,cpf[is.na(cpf)==FALSE])
          cplength <- length(cpfwna[(cpfwna>cpcand-H[i]/step)&(cpfwna<cpcand+H[i]/step)])
          # Checking, if there CPs without other CPs in their neighbourhood.
          if ((cplength ==0) & (length(cpcand)>0))
          {
            cpf[acp+1] <- cpcand; hvalue[acp+1] <- H[i]; acp <- acp+1
          } #end-if
        } #end-for
      } #end-for
    
    # Only CPs on hilfsv
    hvalue <- na.omit(hvalue); hilfsv <- sort(cpf)
    # Corresponding window sizes are saved.
    rhvalue <- rep(0,length(hilfsv))
    if (length(hilfsv)>0)   
      for (i in 1:length(hilfsv))
      {
        rhvalue[i] <- hvalue[cpf==hilfsv[i]][1]
      } #end-for
    invisible(list(hilfsv=hilfsv,rhvalue=rhvalue))
  }#end-cpa
  
  ###
  ### Calculate variances between cp
  ###
  
  var.betCP <- function(Phi,hilfsv,step,muz, S.old)
  {
    lifet <- diff(Phi); TL <- length(lifet)
    for (k in 1:length(hilfsv))
    {
      ma <- hilfsv[k]*step
      maalt <- hilfsv[max(1,k-1)]*step
      manew <- hilfsv[min(k+1,length(hilfsv))]*step
      rma1alt <- which (Phi > maalt)[1] ; rma1new <- which (Phi > manew) [1];       realma<-which (Phi > ma) [1]
      if (is.na(realma)==TRUE) realma <- length(Phi)
      if (is.na(rma1new)==TRUE) rma1new <- length(Phi)
      if (is.na(rma1alt)==TRUE) rma1alt <- length(Phi)
      
      if (k==1) quad1 <- sum(na.omit((lifet[1:max(1,(realma-2))]-muz[1:max(1,(realma-2))])^2))/(realma-2)
      else quad1 <- sum(na.omit((lifet[rma1alt:realma-2]-muz[rma1alt:realma-2])^2))/(realma-rma1alt-2) #
      if (k==length(hilfsv)) quad2 <- sum(na.omit((lifet[(realma):(TL-1)]-muz[(realma):(TL-1)])^2))/(TL-realma-2) #
      else quad2 <- sum(na.omit((lifet[realma:rma1new-2]-muz[realma:rma1new-2])^2))/(rma1new-realma-2)#
      
      varianz[3*k-2] <- ma; varianz[3*k-1] <- quad1; varianz[3*k] <- quad2
      
    } #end-for
    return(varianz)
  }#end-var.betCP
  
  ###
  ###begin main part
  ###
  
  Ght <- sapply(as.matrix(H),FUN=ght.mat,Phi=Phi,S=S,E=E,d=step,
              tc=tc,est.mean=muz) # Derive Ght for all h and t
  if (length(H)==1) {Ght <- list(Ght[,1])}
  for (i in 1:length(H)) #Fill with NA
  {
    ml <- length(rep(NA,1*(H[i]/step)))
    GhtAll[i,1:ml] <- rep(NA,ml); GhtAll[i,((ml+1):(ml+length(Ght[[i]])))] <- Ght[[i]]
  } #end-for
  
  ###
  ###Plotting of Ght-processes
  ###
  
  if (plot.CPD==TRUE) 
  {
    layout(c(1,1,2))
    plot.ght(ght=GhtAll,H=H,step=step,colors=colors,Q=Q,cex.legend=cex.legend,plot.Q=plot.Q,plot.M=plot.M,plot.h=plot.h,ylab1=ylab1)  
  } #end-if
  
  result <- 1; G <- abs(GhtAll)
  if (max(G,na.rm=TRUE)<Q)   
  {
    result <- 0
    invisible(list(cp=NA))
    if (plot.CPD==TRUE) {if (main==TRUE) {title(main=list("0 variance change points detected",cex=1.0*cex.legend))}
      if (plot.var==TRUE)  
      { 
        Dz <- GhtAll[1,] ; Dz <- Dz[is.na(Dz)==FALSE]
        xl <- length(Dz)+2*H[1]/step
        plot.variances(Phi,S,E,S.old=S.old,wid=wid,tc,muz=muz,l1=0,lg=NA,varianz=NA,
                       true.var=NA,xl,H,step,colors,hvalue=NA,cex.legend=cex.legend,ylab2=ylab2)
        result3 <- "No variance CP found!"; result2 <- NA
        const.var <- mean((diff(Phi)-muz)^2)
      } #end-if
    }#end-if
  } #end-if
  else { 
    l1 <- 0; lg <- NA; true.var <- NA
    if (perform.CPD==TRUE)
    {  
      # Change-Point-Detection
      cp <- cpa(GhtAll,Q,H,step) 
      rhvalue <- cp$rhvalue; hilfsv <- cp$hilfsv; varianz <- 0
      # Estimate variances between detected CPs
      if (length(hilfsv)==0) realma <- 0 
      else  varianz<-var.betCP(Phi,hilfsv,step,muz,S.old)
      
      l1 <- length(varianz[varianz!=0])/3 ;lg <- rep(NA,l1); pos.h <- rep(NA,l1); true.var <- rep(NA,l1)
      
      # Draw the CPs
      if (l1==0) l1 <- 0 else
        for (k in 1:l1) 
        {
          pos.h[k] <- which(H==rhvalue[k]); true.var[k]<-varianz[3*k-2]
          lg[k] <- varianz[3*k-2]/step
          if (plot.CPD==TRUE) 
          {
            points(lg[k],abs(GhtAll[pos.h[k],round(lg[k])]),col=1,cex=cex.diamonds,lwd=1)
            points(lg[k],-0.3,pch=18,col=colors[pos.h[k]],xpd=TRUE,cex=cex.diamonds+0.4)  
            points(lg[k],-0.3,pch=5,col=1,xpd=TRUE,cex=cex.diamonds) 
            if (main==TRUE)
            {
              if (l1==1) {title(main=list("1 variance change point detected",cex=1.0*cex.legend))}
              else {title(main=list(paste(l1," variance change points detected"),cex=1.0*cex.legend))}
            } #end-if
          } #end-if
          
        } #end-for
      if (plot.var==FALSE) axis(1,S+c(0,lg,E/step),round(c(round(S,2),round(lg,2),round(E/step,2))*step,2),cex.axis=cex.legend)
      
      # Results are written in a result matrix.
      result2 <- matrix(0,l1+1,4)  
      if (3*l1>1) 
        for (i in 1: (ceiling(l1))) 
        {
          result2[i,] <- c(varianz[(3*i-2):(3*i)],rhvalue[i])
          result2[i,1] <- result2[i,1]+S.old
        } #end-if
      else
        result2[1,3] <- var.window((E-S)/2+S,Phi,(E-S)/2-wid/2,rcp,est.mean)
      
      # Writing the corresponding window sizes in the matrix
      result2 <- result2[-dim(result2)[1],]
      result2 <- signif(result2,3)
      first <- c("variance CP","var bef.","var after","h")
      if (is.vector(result2)) {names(result2) <- first}
      else {colnames(result2) <- first}
    }
    if (plot.CPD==TRUE&plot.var==TRUE)   # If necessary: Estimated variances in wid-intervals are calculated
    {  
      Dz <- GhtAll[1,]; Dz <- Dz[is.na(Dz)==FALSE]; xl <- length(Dz)+2*H[1]/step
      plot.variances(Phi,S,E,S.old=S.old,wid,tc,muz=muz,l1,lg,varianz,true.var,xl,H,step,colors,rhvalue,cex.legend,ylab2=ylab2)
    } #end-if
  }
  M <- max(abs(GhtAll),na.rm=TRUE); const.var <- mean((diff(Phi)-muz)^2,na.rm=TRUE)
  Tt <- signif(E-S,2)
  
  if (print.output) { # Print output in a short report
    if ((is.na(result2[1]))|(is.vector(result2))) result2 <- t(as.matrix(result2)); 
    cat("","\n")
    cat("MFT variances table",sep="\n"); cat("","\n")
    cat(paste("Tt = ",Tt, ", d = ",round(step,2)," and H = {",paste(round(H,2),collapse=", "),"}",sep="")); cat("","\n")
    if(M>Q){cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(Q,2),sep=""),sep="\n")}	else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(Q,2),sep=""),sep="\n")}; cat("","\n")
    if(perform.CPD){
      cat("CPD was performed: ")
      if(M<Q){cat("No change points detected")}
      if((dim(result2)[1]==1)&(M>Q)){cat(paste(dim(result2)[1],"change point detected at "))} 
      if((dim(result2)[1]>=2)&(M>Q)){cat(paste(dim(result2)[1],"change points detected at "))}
      if((dim(result2)[1]>0)&(M>Q)){cat(paste(round(result2[,1])),sep=", ")}; cat("","\n")
      if(M<Q){cat(paste("The estimated variance is",signif(const.var,3) ));  cat("","\n")} else{cat("The estimated variances are ")
        if(perform.CPD) {cat(paste(signif(c(result2[,2],result2[dim(result2)[1],3]),3)),sep=", ")}
        cat("","\n")
      }
    }#end-if-perform.CPD	
    if(!perform.CPD){cat("CPD was not performed")
      cat("","\n"); cat(paste("The estimated variance is",signif(const.var,3) ))}
  }#end-if-print.output  
  
  M <- max(abs(GhtAll),na.rm=TRUE); names(M) <- "test statistic M"; va <- "not estimated"
  if ((perform.CPD)&(M>Q)) {va <- signif(c(unname(result2[,2]),unname(result2[dim(result2)[1],3])),3)}
  if (M>Q) {invisible(list(CP=result2,var=va,varQ=Q,M=M,H=signif(H,2),alpha=alpha,d=step,sim=sim,S=signif(S.old,2),E=signif(E.old,2)))}
  else {invisible(list(CP=result3,var=signif(const.var,3),varQ=Q,M=M,H=signif(H,2),alpha=alpha,d=step,sim=sim,S=signif(S.old,2),E=signif(E.old,2)))}
  
}
