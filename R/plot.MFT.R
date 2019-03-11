#' plot.MFT
#'
#' Plot method for class 'mft'.
#'
#' @param x object of class MFT
#' @param col 	"gray" or vector of colors of length(H). Colors for (G_ht) plot, default: NULL -> rainbow colors from blue to red 
#' @param ylab1 character, ylab for 1. graphic
#' @param ylab2 character, ylab for 2. graphic 
#' @param cex.legend 	numeric, size of annotations in plot
#' @param cex.diamonds numeric, size of diamonds that indicate change points
#' @param main 	logical, indicates if title and subtitle are plotted
#' @param plot.Q 	logical, indicates if rejection threshold Q is plotted
#' @param plot.M	logical, indicates if test statistic M is plotted
#' @param plot.h 	logical, indicates if a legend for the window set H is plotted
#' @param breaks integer, >0, number of breaks in rate histogram	
#' @param wid integer, >0, width of bars in variance histogram	
#' @param ... additional parameters	
#'
#'
#' @examples 
#' # Rate change detection in Poisson process 
#' # with three change points (at t = 250, 600 and 680)
#' set.seed(0)
#' Phi1 <- runif(rpois(1,lambda=390),0,250)
#' Phi2 <- runif(rpois(1,lambda=380),250,600)
#' Phi3 <- runif(rpois(1,lambda=200),600,680)
#' Phi4 <- runif(rpois(1,lambda=400),680,1000)
#' Phi  <- sort(c(Phi1,Phi2,Phi3,Phi4)) 
#' mft  <- MFT.rate(Phi)
#' plot(mft)
#' 
#' 
#' @seealso \code{\link{MFT.rate}, \link{MFT.variance}, \link{MFT.mean}, \link{MFT.peaks}, \link{summary.MFT}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Marietta Kirchner, Julia Schiemann, Jochen Roeper, Ralph Neininger and Gaby Schneider (2014).
#' A multiple filter test for the detection of rate changes in renewal processes with varying variance. The Annals of Applied Statistics 8(4): 2027-67
#' <doi:10.1214/14-AOAS782>
#'
#' 
#' @rdname plot.MFT
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


plot.MFT<-function(x,col=NULL,ylab1=NULL,ylab2=NULL,cex.legend=1.2,cex.diamonds=1.4,
                   main=TRUE,plot.Q=TRUE,plot.M=TRUE,plot.h=TRUE,breaks=NULL,wid=NULL,...)
{
  MFT<-x
  S<-MFT$S; E<-MFT$E; CP<-MFT$CP; H<-MFT$H; alpha<-MFT$alpha; sim<-MFT$sim; M<-MFT$M; rescale<-MFT$rescale
  
  if(is.null(col)){col <- rainbow(length(H),start=2/3,end=0)} else{if(length(col) == 1 & col[1] == "gray"){col <- gray.colors(length(H),0.8,0)}
    else{col <- col; if(length(col)!=length(H)){stop("Argument col must be either NULL, ''gray'', or a vector of colors of length of H")}}
  } # Vector of colors for different windows  
  
  if (MFT$type=="mean")
  {
    X<-MFT$tech.var$X; 
    if (rescale) R_ht<-MFT$tech.var$R_ht
    else R_ht<-MFT$tech.var$G_ht
  
    Q<-MFT$Q; Tt<-MFT$Tt; means<-MFT$expectation

    plot.cpd <- function(X,R_ht,S,E,H,Q,M,CP,means,col,ylab1,ylab2,plot.Q,cex.legend,cex.diamonds,main,plot.h,rescale){
        
      if (MFT$perform.CPD) col.cp <- apply(as.matrix(CP[,2]),MARGIN=1, function(x){col[which(x == H)]}) # Vector that codes CP <-> col.
        
        # Graphic 1: R_ht	
        if(rescale){mi <- min (  c(min( sapply(R_ht, min ) , -7 ) )) - 1 }else{mi <- 0}
        ma <- max (  c(max( sapply(R_ht, max ) ,  12 ) )) + 1 
        
        layout(c(1,1,2))
        par(cex=1,mar=c(0.5,4,4,0.5))
        plot(1,type="n",ylim=c(mi,ma),xlim=c(S,E+(1/10)*(E-S)),xlab="",ylab="",axes=FALSE,cex.lab=cex.legend) 
        axis(2,cex.axis=cex.legend); mtext(ylab1,side=2,cex=cex.legend,padj=-2.5)
        
        for(i in 1:length(H)){ # R_ht-processes
          lines(seq(S+H[i],(E-H[i]),1),R_ht[[i]],col=col[i])
        }
        # End-for-i	
        
        if(plot.M){
          h_max <- which.max(sapply(R_ht,max)); x_max <- which.max(R_ht[[h_max]]); y_max <- R_ht[[h_max]][x_max]
          points((x_max-1)+H[h_max]+S,y_max,pch=4,cex=cex.legend); text((x_max-1)+H[h_max]+S,y_max,"M",pos=2,cex=cex.legend)
        }
        # End-if-plot.M	
        
        if(plot.Q){ # Rejection threshold Q
          lines(c(S,E),rep(Q,2),lty="dashed"); text(S+15,Q,"Q",pos=3,cex=cex.legend)
        }#end-if-plot.Q
        
        if (MFT$perform.CPD) 
        {
        if(dim(CP)[1] > 0){ # CPs 
          for(i in 1:(dim(CP)[1]) ){
            points(CP[i,1]+S,R_ht[[which(CP[i,2] == H)]][((CP[i,1]-(CP[i,2]-1)))],col= 1,cex=cex.diamonds) 
            points(CP[i,1]+S,mi,col=col[which(CP[i,2] == H)],pch=18,cex=cex.diamonds+0.4) 
            points(CP[i,1]+S,mi,col=1,pch=5,cex=cex.diamonds) 
          }#end-for
        }#end-if
        }#end-if
        
        if(is.logical(main)){ # Title
          if(main){ 
            if (MFT$perform.CPD) {if(dim(CP)[1] == 1)	{title(paste(dim(CP)[1]," change point detected"))} else{title(paste(dim(CP)[1]," change points detected"))}} 
            mtext(bquote(M %~~% .(round(M,2)) ),adj=0,padj=-0.5, cex= cex.legend )
            mtext(bquote(Q %~~% .(round(Q,2)) ),adj=1,padj=-0.5,cex=cex.legend )
          }#end-if 
        }#end-if
        
        if(plot.h != FALSE){ # Legend for h_i
          if(plot.h){plot.h <- "topright"}
          legend(x=plot.h,inset=c(0,0.02),legend=as.character(H),pch=19,col=col,cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),bty="n",title=expression(paste(h[i], " =")),pt.cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),title.adj=0.2)
        }#end-if
        
       
        # Graphic 2: Process X.
          par(mar=c(2,4,0.5,0.5))
          plot(X,axes=FALSE,xlim=c(0,length(X)+(1/10)*(Tt)),ylab = "",xlab = "Time",main="",cex.lab=cex.legend, cex = 0.5,pch = 16) 
          mtext(ylab2,side=2,cex=cex.legend,padj=-3.8)
          axis(2,cex.axis=cex.legend)
          axis(1,cex.axis=cex.legend,seq(0,length(X),length=5),round(seq(S,E,length=5))) 
          if (MFT$perform.CPD) {
          if(length(CP[,1])>0){
            for(i in 1:length(CP[,1])){
              arrows((CP[i,1]),max(X)*1.3,(CP[i,1]),max(X)*0.7,col=col.cp[i],length=0.1,lwd=2)	
            }#end-for
          }#end-if
          meanTimes <- diff(c(0,as.vector(CP[,1]),length(X))) 
          doubleMean <- rep(means,times = meanTimes)               # Mean vector for y values.
          lines(seq(from=1,to =Tt,by=1),doubleMean,type="s",col=2,lwd=3)     # Plots lines according to mean between change points.
          }#end-if

      }
      
      # End-plot.cpd
      if(is.null(ylab1) & rescale){ylab1  <- expression(paste(R[list(h,t)]))}
      if(is.null(ylab1) & !rescale){ylab1  <- expression(paste("|",G[list(h,t)],"|",sep =""))}
      if(is.null(ylab2)){ylab2  <- "X"}
      
      # Perform plot:
      plot.cpd(X=X,R_ht=R_ht,S=S,E=E,H=H,M=M,Q=Q,CP=CP,means=means,col=col,ylab1=ylab1,ylab2=ylab2,plot.Q=plot.Q,cex.legend=cex.legend,cex.diamonds=cex.diamonds,main=main,plot.h=plot.h,rescale=rescale)
      
   
  }
  if (MFT$type=="rate")
  {
    d<-MFT$d; Phi<-MFT$tech.var$Phi
    if (rescale) R_ht<-MFT$tech.var$R_ht
    else R_ht<-MFT$tech.var$G_ht
    Tt<-MFT$Tt; rescale<-MFT$rescale; Q<-MFT$Q; rate<-MFT$rate
    if (is.null(breaks)){breaks <- seq(S,E,length=round(Tt/H[1]))}

    plot.cpd.rate <- function(Phi,R_ht,S,E,d,H,Q,M,CP,rate,col,ylab1,ylab2,plot.Q,cex.legend,cex.diamonds,main,plot.h,breaks,rescale){
        
      if (MFT$perform.CPD)  {col.cp <- apply(as.matrix(CP[,2]),MARGIN=1, function(x){col[which(x == H)]})} # Vector that codes CP <-> col
        
        # Graphic 1: R_ht	
        
        if(rescale){mi <- min (  c(min( sapply(R_ht, min ) , -7 ) )) - 1 }else{mi <- 0}
        ma <- max (  c(max( sapply(R_ht, max ) ,  Q ) )) + 1  
        
        layout(c(1,1,2))
        par(cex=1,mar=c(0.5,4,4,0.5))
        plot(1,type="n",ylim=c(mi,ma),xlim=c(S,E+(1/10)*(E-S)),xlab="",ylab="",axes=FALSE,cex.lab=cex.legend) 
        axis(2,cex.axis=cex.legend); mtext(ylab1,side=2,cex=cex.legend,padj=-2.5)
        
        
        for(i in 1:length(H)){ # R_ht-processes
          lines(seq(S+H[i],(E-H[i]),d),R_ht[[i]],col=col[i])
        }#end-for-i	
        
        if(plot.M){
          h_max <- which.max(sapply(R_ht,max)); x_max <- which.max(R_ht[[h_max]]); y_max <- R_ht[[h_max]][x_max]
          points((x_max-1)*d+H[h_max]+S,y_max,pch=4,cex=cex.legend); text((x_max-1)*d+H[h_max]+S,y_max,"M",pos=2,cex=cex.legend)
        }#end-if-plot.M	
        
        if(plot.Q){ # Rejection threshold Q
          lines(c(S,E),rep(Q,2),lty="dashed"); text(S+15,Q,"Q",pos=3,cex=cex.legend)
        }#end-if-plot.Q
        
        if (MFT$perform.CPD) {
        if(dim(CP)[1] > 0){ # CPs 
          for(i in 1:(dim(CP)[1]) ){
            points(CP[i,1],R_ht[[which(CP[i,2] == H)]][((CP[i,1]-S-(CP[i,2]-d))/d)],col= 1,cex=cex.diamonds)
            points(CP[i,1],mi,col=col[which(CP[i,2] == H)],pch=18,cex=cex.diamonds+0.4)
            points(CP[i,1],mi,col=1,pch=5,cex=cex.diamonds)
          }#end-for
        }#end-if
        }#end-if
        
        if(is.logical(main)){ # Title
          if(main){ 
            if (MFT$perform.CPD) {if(dim(CP)[1] == 1)	{title(paste(dim(CP)[1]," rate change point detected"))} else{title(paste(dim(CP)[1]," rate change points detected"))}} 
            mtext(bquote(M %~~% .(round(M,2)) ),adj=0,padj=-0.5, cex= cex.legend )
            mtext(bquote(Q %~~% .(round(Q,2)) ),adj=1,padj=-0.5,cex=cex.legend )
          }#end-if 
        }#end-if
        
        if(plot.h != FALSE){ # Legend for h_i
          if(plot.h){plot.h <- "topright"}
          legend(x=plot.h,inset=c(0,0.02),legend=as.character(H),pch=19,col=col,cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),bty="n",title=expression(paste(h[i], " =")),pt.cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),title.adj=0.2)
        }#end-if
        
        
        # Graphic 2: Rate histogram
        par(mar=c(2,4,0.5,0.5))
        if (is.null(breaks)){breaks <- seq(S,E,length=round(Tt/H[1]))}
        histo 	  		<-  hist(Phi,breaks=breaks,plot=FALSE)
        bin       		<-  diff(histo$breaks[1:2])
        histo$counts  	<-  histo$counts / bin
        plot(histo,xlab="",ylab=ylab2,main="",axes=FALSE,xlim=c(S,E+(1/10)*(E-S)),ylim=c(0,1.3*max(histo$counts)),cex.lab=cex.legend) 
        if (MFT$perform.CPD) {axis(1,at=c(S,CP[,1],E),cex.axis=cex.legend); axis(2,cex.axis=cex.legend)}
        else {axis(1,at=c(S,E),cex.axis=cex.legend); axis(2,cex.axis=cex.legend) }
        if (MFT$perform.CPD) {
        lines(c(S,CP[,1],E),c(rate,rate[length(rate)]),type="s",col=2,lwd=2)
        if(length(CP[,1] > 0)){
          for(i in 1:length(CP[,1])){
            arrows(CP[i,1],2*max(histo$counts),CP[i,1],max(rate[i:(i+1)])*1.15,col=col.cp[i],length=0.1,lwd=2)	
          }#end-for
        }#end-if
        }#end-if
        
      }#end-plot.cpd.rate
      
      if(is.null(ylab1) & rescale){ylab1  <- expression(paste(R[list(h,t)]))}
      if(is.null(ylab1) & !rescale){ylab1  <- expression(paste("|",G[list(h,t)],"|",sep =""))}
      if(is.null(ylab2)){ylab2  <- "rate"}
      
      # Perform plot:
      plot.cpd.rate(Phi=Phi,R_ht=R_ht,S=S,E=E,d=d,H=H,M=M,Q=Q,CP=CP,rate=rate,col=col,ylab1=ylab1,ylab2=ylab2,plot.Q=plot.Q,cex.legend=cex.legend,cex.diamonds=cex.diamonds,main=main,plot.h=plot.h,breaks=breaks,rescale=rescale)
      
 }
  if (MFT$type=="variance")
  {
    tech.var<-MFT$tech.var; varQ<-MFT$varQ;   step<-MFT$d; d<-MFT$d; G_ht<-MFT$tech.var$G_ht; Tt<-MFT$Tt; Phi<-MFT$tech.var$Phi
   
    bar <- (mean(diff(Phi)))*100
    if (is.null(wid)) wid <- floor(max(Phi))/(floor(floor(max(Phi))/bar)) 
    if (wid>((E-S)/2)) {stop("width of variance bars too large")}
    
    if(is.null(ylab1)) {ylab1=expression(abs(G[list(h,t)]))}
    if(is.null(ylab2)) {ylab2=expression(widehat(sigma)^2)}
    
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
    ### Function to plot G_ht-processes for every h and t
    ###
    
    plot.ght <- function(ght,H,E,S,step,col,Q,cex.legend,plot.Q=TRUE,plot.M=TRUE,plot.h=TRUE,ylab1) {
      max.ght <- max(abs(ght),na.rm=TRUE)
      xm <- c(); x_max <- c();         lDz<-(E-S)/step
      par(mar=c(0.5,4,4,0.5),cex=1)
      for (i in 1:length(H))
      {
        border <- 1.2*max(Q,max.ght) ;
        # Save max values for each i to determine M    
        g <- na.omit(abs(ght[i,])); xm[i] <- which(g==max(g))[1]; x_max[i] <- abs(g[xm[i]])
        if (i==1) #For i=1 create new plot, otherwise add lines
        { 
          plot(abs(ght[i,]), main="", xlab="", ylab="", type="l", ylim=c(-0.6,border),xaxt="n",
               xlim=c(0,1.1*lDz), cex.axis=cex.legend, 
               cex.lab=cex.legend, bty="n", col=col[1],yaxt="n")
          lang <- length(ght[i,])
          mtext(ylab1,side=2,padj=-2.6,cex=cex.legend)
          if (plot.h) {legend(x="topright",inset=c(0,0.02),legend=as.character(round(H,0)),pch=19,col=col,cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),bty="n",title=expression(paste(h[i], " =")),pt.cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),title.adj=0.2,xpd=TRUE)}
          if (plot.Q) {
            if (main) mtext(bquote(Q %~~% .(round(Q,2)) ),adj=1,padj=-0.5,cex=cex.legend)
            lines(c(0,1.0*lDz),c(Q,Q),col="black", lty=2, lwd=1); text(S+15,Q,"Q",pos=3,cex=1.0*cex.legend)
          } #end-if
        } #end-if
        if (i>1) lines(abs(ght[i,]),col=col[i])
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
    ### Function to plot variances in small windows
    ###
    
    plot.variances <- function(Phi,S,E,S.old,wid,tc,muz,ncp,lg,varianz,true.var,xl,H,step,col,hvalue,cex.legend,ylab2,perform.CPD)
    {    
      zeitpl <- seq(S+(wid/2),E-(wid/2),by=wid)
      lzpl <- length(zeitpl)
      zeitpl <- seq(S+(wid/2),E-(wid/2),by=wid)
      varf <- sapply(zeitpl,var.window,Phi=Phi,h=wid,tc=tc,est.mean=muz)
      par(mar=c(2,4,0.5,0.5),cex=1)
      yl <- 1.2*max(varf[is.na(varf)==FALSE & is.finite(varf)==TRUE])
      plot(1,main="",xlab="",ylim=c(-0, yl),
           ylab="",type="n", xaxt="n", cex.axis=cex.legend, cex.lab=1.0, bty="n",
           xlim=c(0,0+1.1*(lzpl)*xl),yaxt="n")
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
      if (perform.CPD) {
      if (ncp==0) lines(c(0,length(varf[is.na(varf)==FALSE])*xl),rep(const.var,2),lwd=2, col="red") 
      else 
        for (k in 1:ncp) 
        {
          lg[k] <- ((varianz[3*k-2])*xl*lzpl)/E
          if (varianz[3*k-1]<=varianz[3*k]) 
          { 
            arrows(lg[k],varianz[3*k-1],lg[k],varianz[3*k],col="red", lwd=2, length=0)
            arrows(lg[k],1.8*max(varf),lg[k],1.1*varianz[3*k],col=col[which(H==hvalue[k])],xpd=TRUE,length=0.1,lwd=2)
          }
          else 
          {
            arrows(lg[k],varianz[3*k-1],lg[k],varianz[3*k],col="red", lwd=2, length=0)
            arrows(lg[k],1.8*max(varf),lg[k],1.1*varianz[3*k-1],col=col[which(H==hvalue[k])],xpd=TRUE,length=0.1,lwd=2)
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
      if (ncp>0) lines(c((varianz[3*ncp-2])*xl*lzpl/E,xl*lzpl),
                      c(rep(varianz[3*ncp],2)),lwd=2,col="red") 
      }#end-if
      axis(1,S+c(0.0*xl,lg,length(varf[is.na(varf)==FALSE])*xl),c(round(S+S.old,2),round(true.var+S.old,2),round(E+S.old,2)),cex.axis=cex.legend)
      
    }#end-plot.variances
   xl<- length(G_ht[1,][is.na(G_ht[1,]==FALSE)])+2*H[1]/step
   layout(c(1,1,2))
   plot.ght(ght=G_ht,H=H,E=E,S=S,step=step,col=col,Q=varQ,cex.legend=cex.legend,plot.Q=plot.Q,plot.M=plot.M,plot.h=plot.h,ylab1=ylab1)  
   
  
  #Plot variances in small windows    
   if (M<varQ)   
    {
      if (main==TRUE) {title("0 variance change points detected")
        plot.variances(Phi,S,E,S.old=tech.var$S.old,wid=wid,tc=tech.var$tc,muz=tech.var$muz,ncp=0,lg=NA,varianz=NA,
                         true.var=NA,xl=xl,H=H,step=step,col=col,hvalue=NA,cex.legend=cex.legend,ylab2=ylab2,perform.CPD=MFT$perform.CPD)
                     }#end-if
    } #end-if
    else { 
          ncp <- ifelse(is.null(dim(MFT$CP)),0,dim(MFT$CP)[1])[1]
          lg <- c(); pos.h <- c(); true.var <- c()
          # Draw the CPs
          if (MFT$perform.CPD) {
          if (ncp>0)  {
            varianz<-as.vector(t(MFT$CP[,1:3]))
            varianz[seq(1,length(varianz)-2,by=3)]<-varianz[seq(1,length(varianz)-2,by=3)]-tech.var$S.old  
            for (k in 1:ncp) 
            {
              pos.h[k] <- which(H==MFT$CP[k,4]); true.var[k]<-varianz[3*k-2]
              lg[k] <- varianz[3*k-2]/step
              points(lg[k],abs(G_ht[pos.h[k],round(lg[k])]),col=1,cex=cex.diamonds,lwd=1)
              points(lg[k],-0.3,pch=18,col=col[pos.h[k]],xpd=TRUE,cex=cex.diamonds+0.4)  
              points(lg[k],-0.3,pch=5,col=1,xpd=TRUE,cex=cex.diamonds) 
              if (main==TRUE)
              {
                if (ncp==1) {title("1 variance change point detected")}
                else {title(paste(ncp," variance change points detected"))}
              } #end-if
            } #end-for
        } #end-if
          }#end-if
        plot.variances(Phi,S,E,S.old=tech.var$S.old,wid=wid,tc=tech.var$tc,muz=tech.var$muz,ncp=ncp,lg=lg,varianz=varianz,col=col,H=H,xl=xl,hvalue=MFT$CP[,4],true.var=true.var,cex.legend=cex.legend,ylab2=ylab2,perform.CPD=MFT$perform.CPD)
    
      }
  }
  
  if (MFT$type=="peaks")
  {
    R_ht<-MFT$tech.var$D_ht; x<-MFT$tech.var$x
    Tt<-MFT$Tt;  Q<-MFT$Q; 
    if (MFT$perform.CPD) CP[,1]<-CP[,1]-S
    
    
      plot.cpd.peaks <- function(x, R_ht, S, E, H, Q, M, CP, 
                           col, ylab1, ylab2, plot.Q, cex.legend, cex.diamonds, 
                           main, plot.h, Tt) 
      {
        if (is.null(col)) {
          col <- rainbow(length(H), start = 2/3, end = 0)
        }
        else {
          if (length(col) == 1 & col[1] == "gray") {
            col <- gray.colors(length(H), 0.8, 0)
          }
          else {
            col <- col
            if (length(col) != length(H)) {
              stop("Argument col must be either NULL, ''gray'', or a vector of colors of length of H")
            }
          }
        }
        if (MFT$perform.CPD) {col.cp <- apply(as.matrix(CP[, 2]), MARGIN = 1, function(x) {
          col[which(x == H)]
        })}
        mi <- min(c(min(sapply(R_ht, min), -7))) - 1
        ma <- max(c(max(sapply(R_ht, max), 12))) + 1
        layout(c(1, 1, 2))
        par(cex = 1, mar = c(0.5, 4, 4, 0.5))
        plot(1, type = "n", ylim = c(mi, ma), xlim = c(S, 
                                                       E + (1/10) * (E - S)), xlab = "", ylab = "", 
             axes = FALSE, cex.lab = cex.legend)
        axis(2, cex.axis = cex.legend)
        mtext(ylab1, side = 2, cex = cex.legend, padj = -2.5)
        for (i in 1:length(H)) {
          lines(seq(S + H[i]-1, (E - H[i]), 1), R_ht[[i]], 
                col = col[i])
        }
        if (plot.M) {
          h_max <- which.max(sapply(R_ht, max))
          x_max <- which.max(R_ht[[h_max]])
          y_max <- R_ht[[h_max]][x_max]
          points((x_max - 1) + H[h_max] + S, y_max, pch = 4, 
                 cex = cex.legend)
          text((x_max - 1) + H[h_max] + S, y_max, "M", 
               pos = 2, cex = cex.legend)
        }
        if (plot.Q) {
          lines(c(S, E), rep(Q, 2), lty = "dashed")
          text(S + 15, Q, "Q", pos = 3, cex = cex.legend)
          if(MFT$two.sided){lines(c(S,E),rep(-Q,2),lty = "dashed")}
        }
        
        
          if (is.logical(main)) {
            if (main) {
                        mtext(bquote(M %~~% .(round(M, 2))), adj = 0, 
                    padj = -0.5, cex = cex.legend)
              mtext(bquote(Q %~~% .(round(Q, 2))), adj = 1, 
                    padj = -0.5, cex = cex.legend)
            }
          }
          
        if (plot.h != FALSE) {
            if (plot.h) {plot.h <- "topright"}
            legend(x = plot.h, inset = c(0, 0.02), legend = as.character(H), 
                   pch = 19, col = col, cex = min(cex.legend, 
                                                  cex.legend + 0.5 - 0.05 * length(H)), bty = "n", 
                   title = expression(paste(h[i], " =")), pt.cex = min(cex.legend, 
                                                                       cex.legend + 0.5 - 0.05 * length(H)), title.adj = 0.2)
          }
        
        if(any(is.na(CP))){CP <- matrix(numeric(0),nrow=0,ncol=2)}
        if (dim(CP)[1] > 0) {
          for (i in 1:(dim(CP)[1])) {
            points(CP[i, 1] + S, R_ht[[which(CP[i, 2] == 
                                               H)]][((CP[i, 1] - (CP[i, 2] - 1)))], col = 1, 
                   cex = cex.diamonds)
            points(CP[i, 1] + S, mi, col = col[which(CP[i, 
                                                        2] == H)], pch = 18, cex = cex.diamonds + 
                     0.4)
            points(CP[i, 1] + S, mi, col = 1, pch = 5, 
                   cex = cex.diamonds)
          }
        }
        if (is.logical(main)) {
          if (main) {
            if (dim(CP)[1] == 1) {
              title(paste(dim(CP)[1], " peak detected"))
            }
            else {
              title(paste(dim(CP)[1], " peaks detected"))
            }
            mtext(bquote(M %~~% .(round(M, 2))), adj = 0, 
                  padj = -0.5, cex = cex.legend)
            mtext(bquote(Q %~~% .(round(Q, 2))), adj = 1, 
                  padj = -0.5, cex = cex.legend)
          }
        }
        
        par(mar = c(2, 4, 0.5, 0.5))
        plot(x,type="l", axes = FALSE, xlim = c(0, length(x) + 
                                                    (1/10) * (Tt)), ylab = ylab2, xlab = "Time", 
               main = "", cex.lab = cex.legend, cex = 0.5, 
               pch = 16)
          axis(2, cex.axis = cex.legend)
          axis(1, cex.axis = cex.legend, at = floor(seq(1, length(x),length.out = 5)), floor(seq(S, E, length.out = 5)))
          if (length(CP[, 1]) > 0) {
            for (i in 1:length(CP[, 1])) {
              arrows((CP[i, 1]), max(x) * 1.3, (CP[i, 1]), 
                     max(x) * 0.9, col = col.cp[i], length = 0.1, 
                     lwd = 2)
            }
          }
      }#end-function-plot.cpd.peaks   
      
      if (is.null(ylab1)) {
        ylab1 <- expression(paste(D[list(h, t)],sep = ""))
      }
      if (is.null(ylab2)) {
        ylab2 <- "X"
      }
      
      plot.cpd.peaks(x = x, R_ht = R_ht, S = S, E = E, H = H, M = M, 
               Q = Q, CP = CP, col = col, 
               ylab1 = ylab1, ylab2 = ylab2, plot.Q = plot.Q, plot.h=plot.h,cex.legend = cex.legend, 
               cex.diamonds = cex.diamonds, main = main, Tt=Tt)
    } 

}