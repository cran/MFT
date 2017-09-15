
#' MFT.mean
#'
#' The multiple filter test for mean change detection in time series or sequences of random variables.
#'
#' @param Y numeric vector, input sequence of random variables.
#' @param rescale logical, if TRUE statistic G is rescaled to statistic R 
#' @param S	numeric, start of time interval, default: Smallest multiple of d that lies beyond min(Phi)
#' @param E numeric, end of time interval, default: Smallest multiple of d that lies beyond max(Phi), needs E > S.
#' @param autoset.H	logical, automatic choice of window size H
#' @param H vector, window set H, all elements must be increasing ordered multiples of d, the largest element must be =< (T/2). H is automatically set if autoset.H = TRUE
#' @param alpha numeric, in (0,1), significance level
#' @param sim integer, > 0, No of simulations of limit process (for approximation of Q), default = 10000
#' @param method either "asymptotic" or "fixed", defines how threshold Q is derived, default: "asymptotic", If "asymptotic": Q is derived by simulation of limit process L (Brownian motion); possible set number of simulations (sim), If "fixed": Q may be set automatically (Q)
#' @param Q	numeric, rejection threshold, default: Q is simulated according to sim and alpha.
#' @param perform.CPD logical, if TRUE change point detection algorithm is performed
#' @param plot.CPD logical, if TRUE CPD-scenario is plotted. Only active if perform.CPD == TRUE
#' @param print.output logical, if TRUE results are printed to the console	
#' @param col 	"gray" or vector of colors of length(H). Colors for (R_ht) plot, default: NULL -> rainbow colors from blue to red. 
#' @param ylab1 character, ylab for 1. graphic
#' @param ylab2 character, ylab for 2. graphic 
#' @param cex.legend 	numeric, size of annotations in plot
#' @param cex.diamonds numeric, size of diamonds that indicate change points
#' @param main 	logical, indicates if title and subtitle are plotted
#' @param plot.Q 	logical, indicates if rejection threshold Q is plotted
#' @param plot.M	logical, indicates if test statistic M is plotted
#' @param plot.h 	logical, indicates if a legend for the window set H is plotted
#' @param plot.mean logical, indicates if a legend of estimated rates is plotted
#' @param plot.cp logical, indicates if a legend of detected CPs is plotted
#' @param plot.process logical, indicates if there should be a plot of Y as second graphic.
#'
#' @return invisible
#' \item{M}{test statistic}
#' \item{Q}{rejection threshold}
#' \item{sim}{number of simulations of the limit process (approximation of Q)}
#' \item{CP}{set of change points estmated by the multiple filter algorithm, increasingly ordered in time}
#' \item{rate}{estimated mean values between adjacent change points}
#' \item{SWD}{sets of change points estimated from preprocessing single window detections}
#' \item{S}{start of time interval}
#' \item{E}{end of time interval}
#' \item{H}{window set}
#' \item{alpha}{significance level}
#' 
#' @seealso \code{\link{MFT.rate}, \link{MFT.variance}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' 
#' @references 
#' Michael Messer, Marietta Kirchner, Julia Schiemann, Jochen Roeper, Ralph Neininger and Gaby Schneider (2014).
#' A multiple filter test for the detection of rate changes in renewal processes with varying variance. The Annals of Applied Statistics 8(4): 2027-67
#' <doi:10.1214/14-AOAS782>
#' 
#' @examples 
#' # Normal distributed sequence with 3 change points of the mean (at n=100, 130, 350)
#' Y1 <- rnorm(400,0,1); Y2 <- rnorm(400,3,1); Y3 <- rnorm(400,5,1); Y4 <- rnorm(600,4.6,1)
#' Y  <- c(Y1[1:100],Y2[101:130],Y3[131:350],Y4[351:600])
#' MFT.mean(Y)
#' # Set additional parameters (window set)
#' MFT.mean(Y,autoset.H=FALSE,H=c(40,80,160))
#' 
#' @rdname MFT.mean
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


MFT.mean <-
function(Y,rescale=TRUE,autoset.H=TRUE,S=NULL,E=NULL,H=NULL,alpha=0.05,sim=10000,method="asymptotic",Q=NA,perform.CPD=TRUE,print.output=TRUE,
                plot.CPD=TRUE,col=NULL,ylab1=NULL,ylab2=NULL,cex.legend=1.2,cex.diamonds=1.4,
                main=TRUE,plot.Q=TRUE,plot.M=TRUE,plot.h=TRUE,plot.mean=FALSE,plot.cp=FALSE,
                plot.process = TRUE){
###
### Set Parameters
###
  
  if(is.null(S)  & is.null(E))	{S <- 1; E <- length(Y)}
  if(is.null(S)  & !is.null(E))	{S <- 1}
  if(!is.null(S) & is.null(E))	{E <- length(Y)}	
  if(E-S <= 0){stop("Invalid choice of S and E: Need S < E.")}  
  
  
  if(! method %in% c("aymptotic","fixed")){"Invalid choice of method: method must be 'asymptotic' or 'fixed'"}
  if(method == "fixed" & !is.numeric(Q)){stop("Invalid choice of Q: Q must be positive real number")}
  if(method == "fixed" & !(is.numeric(Q) & Q>0)){cat("Warning: non-positive threshold might be inappropriate. Possibly choose Q > 0",sep="\n") }
  if(method != "fixed" & !is.na(Q)){cat("Warning: Q is derived by simulation. In order to set Q manually, choose: method = 'fixed'",sep="\n")}
  
  


# Begin-if-autoset.H
  if(autoset.H){
    Tt  <- length(Y[S:E]) 											                        # Length of interval.
    frac <- c(1/16,1/8,3/16,1/4,1.25/4,1.5/4)
    H   <-  round(frac*Tt)                                              # Set window sizes.
  }
# End-if-set.H 

# Begin-if-!autoset.H
  if(!autoset.H){
    if(is.null(H)){stop("If autoset.H is FALSE, the vector of window sizes H must be set")}
    Tt  <- length(Y[S:E])	
    if(!all(c(diff(c(1,H,(Tt/2))))>0) | !all(H%%1 == 0)){stop("Invalid choice of window set: H needs to be an increasing ordered vector of integers, with min(H) > 1 and max(H) <= Tt/2")}
  }
# End-if-!autoset.H 

  Y <- Y[S:E]         # Shorten sequence to region of analysis.

  if( !sim%%1==0 | sim <= 0) {stop("Invalid choice of number of simulations: sim must be a positive integer")}
  if(sim < 10000){cat("Warning: Number of simulations might be too low",sep="\n")}
  if( alpha*(1-alpha) <= 0){stop("Invalid choice of significance level: alpha must be in (0,1)")}
  
###
### sim.parameter 
###
  
  sim.parameter <- function(sim=sim,Tt=Tt,H=H,alpha=alpha,rescale=rescale){
# Begin-maxSternh: Given Brownian motion W, max(L_ht) is calculated for fixed h in H.
    maxSternh <- function(h=h,W=W,Tt=Tt){ 
      Wt_plus <- W[(2*h+1):(Tt)]                              # W_(t+h)
      Wt		<- W[(1+h):(Tt-h)]                                # W_(t)
      Wt_minus<- W[1:(Tt-2*h)]	                              # W_(t-h)
      return(max((1/sqrt(2*h))*abs(Wt_plus - 2*Wt + Wt_minus))) # max_(t)|L_(h,t)|
    }
# End-maxSternh

# Begin-sim.maxH: Simulates Brownian motion and calculates max(L_ht) for all h in H.
    sim.maxH <- function(Tt=Tt,H=H){ 
      W  <- c(0,cumsum(rnorm((Tt-1),mean = 0, sd=1)))         # Standard Brownian Motion of length Tt.
      return(vapply(as.matrix(H),FUN=maxSternh,W=W,Tt=Tt,numeric(length(1))))   # Apply the calculation of the maximum of Lht.
    }
# End-sim.maxH

# Begin-sim.maxmax: Simulates Brownian motions and calculates mean(max|L_ht|), 
# sd(max|L_ht|) and Q.   
    sim.maxmax <- function(sim=sim,Tt=Tt,H=H,alpha=alpha,rescale=rescale){ 
      re 		<- matrix(replicate(sim,sim.maxH(Tt=Tt,H=H)),nrow=sim,byrow=TRUE)	  # Simulates the Browian motion sim times and calculates maxLht for every h in H.
      meanH 	<- apply(re,MARGIN=2,mean)       # Calculates the empirical mean of the set of max Lhts for all h in H.
      sdH   	<- apply(re,MARGIN=2,sd)         # Calculates the empirical standard deviation of the Lhts for all h in H.
      if (rescale) maxmax 	<- apply((t(re)-meanH)/sdH, MARGIN=2,max)  # Calculates M* sim times.
      if(!rescale){maxmax 	<- apply((t(re)), MARGIN=2,max)}  
      Q 		<- quantile(maxmax,1-alpha)     # Calculates the empirical alpha-Quantile M*.
      list(Q=Q,meanH=meanH,sdH=sdH)         # Writes the threshold Q, and the vectors of means and sd in a list.
    }
# End-sim.maxmax
    
    sim.maxmax(sim=sim,Tt=Tt,H=H,alpha=alpha,rescale=rescale)     # Calls function sim.maxmax.
    
  }
# End-sim.parameter	
  
  parameter <- sim.parameter(sim=sim,Tt=Tt,H=H,alpha=alpha,rescale=rescale) # Simulates parameter values.
  
  if(method=="fixed"){parameter$Q <- Q}   # If manually set Q not Null, then it is set as rejection threshold.
  
  
  ### 
  ### Calculate R_ht
  ###
  
  
  # Calculates G_ht for fixed t and h.
  ght <- function(t,Y,h){ 
    windowRight <- Y[(t+1):(t+h)]   # Values in right window.
    windowLeft <- Y[(t-h+1):t]
    meanRight <- mean(windowRight)  # Mean in right window.
    meanLeft <- mean(windowLeft)
    varianceRight <- var(windowRight) # Variance in right window.
    varianceLeft <- var(windowLeft)
    if (is.na(varianceLeft) | is.na(varianceRight) | varianceLeft <= 0 | varianceRight <= 0){g <- 0}
    else {  g <- (meanRight - meanLeft) / sqrt((varianceLeft + varianceRight)/h)} # Calculates G_ht
    return(g)
  }
  # End-Ght	
  
  # Calculate R_ht for all h and t.
  Rh <- function(h,Y,H,rescale,meanH,sdH){
    seq2 <- seq(1+h,length(Y)-h,by=1)    # All possible t values for window size h.
    G <- vapply(seq2,FUN = ght,Y=Y,h=h,numeric(1))	# Calculates all G_ht for window size h.
    if(rescale){return( (abs(G) - meanH[which(H==h)]) / sdH[which(H==h)] )}
    if(!rescale){return(abs(G))}
  }
  # End-Rht
 

  R_ht <- lapply(as.matrix(H),FUN=Rh,Y=Y,H=H,rescale=rescale,meanH=parameter$meanH,sdH=parameter$sdH) # R_ht values for every window size.
  Mh 	 <- vapply(R_ht,max,numeric(1))	      # Maximum for every window size.
  M	 <- max(Mh); names(M) <- "Test statistic M"   # Global Maximum M, Test statistic.

  
  ###
  ### Perform Change Point Detection (CPD)
  ###	
  
  if(perform.CPD){
    
    # Perform Single window detection (SWD) for a fixed window size h.
    SWD <- function(R,Tt,Q){ 
      h  <- (Tt - (length(R))) / 2      # Calculates window size h from R_ht.
      CP <- rep(NA,Tt); j <- 1 		      # Vector to safe change poinst. Set index to 1.
      while ( any(R > Q) ) {            # As long as the process R_ht > Q, search for change points.
        c <- which( R == max(R) )[1]    # Search for the index of the maximum value.
        CP[j] <- (c-1) + h              # Calculate inx in actual timeline.
        left 		  <- max (1, c-h+1)     # Boundary left h-neighbourhood.
        right  		  <- min ((Tt-2*h)+1, c+h-1)  # Boundary right h-neighbourhood.
        R[left:right] <- rep(Q-1,length(left:right)) # Set value in h-neighbourhood to insignificant value.
        j <- j + 1                      # Increase counter.
      } 
      # End-while
      if( all(is.na(CP)) ){return(numeric(0))}    # If all CP do not have a value, return 0.
      else{return(as.vector(na.omit(CP)))}        # Else return calculated CPs.
    }
    # End-SWD	
    
    SWD <- lapply(R_ht,FUN=SWD,Tt=Tt,Q=parameter$Q) # SWD for all h in H.
    
    # Combine SWD results.
    
    delete <- function(element,h2,c1){ # If c1 falls within the h2-neighbourhood of element, then element is deleted.
      hit <- any( ifelse( element - h2 < c1  &  c1 < element + h2, TRUE, FALSE) )  
      return(hit)
    }
    # End-delete
    
    outside <- function(c1,c2,h2){ # Outside is logical vector of length c2. TRUE -> Component deleted. FALSE -> Component remains.
      tot_hit <- apply(as.matrix(c2),MARGIN=1,FUN=delete,h2=h2,c1=c1) # For a fixed CP c1 and a fixed window size h2, test all c2s for overlaps.
      return(tot_hit)
    }
    # End-function
    
    # With "outside" run through all "h-areas" and delete "overlap-CPs".
    c		<- SWD[[1]] 					# First, accept CPs detected with smallest h.
    h_val 	<- rep(H[1],length(SWD[[1]])) 	# Corresponding h value.
    
    if(length(H)>1){              # If there are more than 1 window sizes in H, work through all. 
      for (k in 2:(length(H))){   # For windowsize with index 2 and up ...
        if(!is.na(SWD[[k]][1]) ){ # Test if there are CP detected with window size hk.
          put_out <- outside(c,SWD[[k]],H[k]) # Tests if CP in c are in a hk-neighbourhood of the CP in SWD[[k]]. Delete if TRUE.
          c 	 	<- c(c,SWD[[k]][put_out == FALSE]) # All CP which are not deleted, are added to the Set c.
          h_val	<- c(h_val,rep(H[k],length(c)-length(h_val))) # New h-value to test by.
        } # End-if
      } # End-for
    } # End-if-length(H)
    
    CP <- cbind(c,h_val)                        # Safes CPs together with window size of detection.
    if(dim(CP)[1]>1){CP <- CP[order(CP[,1]),]}  # Sorts CPs by their detection h value.
    colnames(CP) <- c("changepoint","h_value")  # Assigns cloumn names.
    if(dim(CP)[1] > 0) {rownames(CP) <- as.character(1:dim(CP)[1])} # Assigns row names as characters.
    
    # Calculats mean between change points.
    CPs <- c(1,as.vector(CP[,1]),length(Y))           # CP plus start and end
    meanCP <- rep(NA,(length(CPs)-1))         # Vector for means.
    for(i in 1:(length(CPs)-1)){      
      meanCP[i] <- mean(Y[CPs[i]:CPs[i+1]])   # Mean for every section. 
    }
    # End-for  
    
    names(meanCP) <- apply(as.matrix(1:length(meanCP)),1,function(index){paste("mu_",index,sep="")})
    
  }
  # End-perform.CPD

  
  ###
  #### plot of scenario
  ###
  
  if(plot.CPD & !perform.CPD){cat("Can not plot scenario. Need perform.CPD == TRUE",sep="\n")}
  
  if(plot.CPD & perform.CPD){
    
    plot.cpd <- function(Y,R_ht,S,E,H,Q,M,CP,meanCP,col,ylab1,ylab2,plot.Q,cex.legend,cex.diamonds,main,plot.h,plot.cp){
      
      if(is.null(col)){col <- rainbow(length(H),start=2/3,end=0)} else{if(length(col) == 1 & col[1] == "gray"){col <- gray.colors(length(H),0.8,0)}
        else{col <- col; if(length(col)!=length(H)){stop("Argument col must be either NULL, ''gray'', or a vector of colors of length of H")}}
      } # Vector of colors for different windows.
      col.cp <- apply(as.matrix(CP[,2]),MARGIN=1, function(x){col[which(x == H)]}) # Vector that codes CP <-> col.
      
      # Graphic 1: R_ht	
      mi <- min (  c(min( sapply(R_ht, min ) , -7 ) )) - 1 
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
      
      if(dim(CP)[1] > 0){ # CPs 
        for(i in 1:(dim(CP)[1]) ){
          points(CP[i,1]+S,R_ht[[which(CP[i,2] == H)]][((CP[i,1]-(CP[i,2]-1)))],col= 1,cex=cex.diamonds) 
          points(CP[i,1]+S,mi,col=col[which(CP[i,2] == H)],pch=18,cex=cex.diamonds+0.4) 
          points(CP[i,1]+S,mi,col=1,pch=5,cex=cex.diamonds) 
        }#end-for
      }#end-if
      
      if(is.logical(main)){ # Title
        if(main){ 
          if(dim(CP)[1] == 1)	{title(paste(dim(CP)[1]," change point detected"))} else{title(paste(dim(CP)[1]," change points detected"))} 
          mtext(bquote(M %~~% .(round(M,2)) ),adj=0,padj=-0.5, cex= cex.legend )
          mtext(bquote(Q %~~% .(round(Q,2)) ),adj=1,padj=-0.5,cex=cex.legend )
        }#end-if 
      }#end-if

      if(plot.h != FALSE){ # Legend for h_i
        if(plot.h){plot.h <- "topright"}
        legend(x=plot.h,inset=c(0,0.02),legend=as.character(H),pch=19,col=col,cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),bty="n",title=expression(paste(h[i], " =")),pt.cex=min(cex.legend,cex.legend+0.5-0.05*length(H)),title.adj=0.2)
      }#end-if
      
      if(plot.mean != FALSE){ # Legend for means
        legend(x=plot.mean,inset=c(0,0.02),legend=as.character(round(meanCP,2)),pch=19,cex=cex.legend,bty="n",title=expression(paste(hat(mean)[i] %~~% "")),title.adj=0.2)
      }#end-if
      
      if(dim(CP)[1] > 0){ # Legend for CPs	
        if(plot.cp != FALSE){
          legend(x=plot.cp,inset=c(0,0.02),legend=as.character(CP[,1]),col=col.cp,pch=18,cex=cex.legend,bty="n",text.col=col.cp,title=expression(paste(hat(c)[i]," = ")),title.col=1,pt.cex=cex.legend+0.4 )
          legend(x=plot.cp,inset=c(0,0.02),legend=as.character(CP[,1]),col=1,pch=5,cex=cex.legend,bty="n",text.col=0,title=expression(paste(hat(c)[i]," = ")),title.col=0,pt.cex=cex.legend ) 
        }#end-if
      }#end-if

     
      # Graphic 2: Process Y.
      if(plot.process==TRUE){
         par(mar=c(2,4,0.5,0.5))
         plot(Y,axes=FALSE,xlim=c(0,length(Y)+(1/10)*(Tt)),ylab = ylab2,xlab = "Time",main="",cex.lab=cex.legend, cex = 0.5,pch = 16) 
         axis(2,cex.axis=cex.legend)
         axis(1,cex.axis=cex.legend,seq(0,length(Y),length=5),round(seq(S,E,length=5))) 
         if(length(CP[,1])>0){
           for(i in 1:length(CP[,1])){
             arrows((CP[i,1]),max(Y)*1.3,(CP[i,1]),max(Y)*0.7,col=col.cp[i],length=0.1,lwd=2)	
           }#end-for
         }#end-if
         meanTimes <- diff(c(0,as.vector(CP[,1]),length(Y))) 
         doubleMean <- rep(meanCP,times = meanTimes)               # Mean vector for y values.
         lines(seq(from=1,to =Tt,by=1),doubleMean,type="s",col=2,lwd=3)     # Plots lines according to mean between change points.
        }
      
     
    }
  
    # End-plot.cpd
    if(is.null(ylab1) & rescale){ylab1  <- expression(paste(R[list(h,t)]))}
    if(is.null(ylab1) & !rescale){ylab1  <- expression(paste("|",G[list(h,t)],"|",sep =""))}
    if(is.null(ylab2)){ylab2  <- "meanCP"}
    
    # Perform plot:
    plot.cpd(Y=Y,R_ht=R_ht,S=S,E=E,H=H,M=M,Q=parameter$Q,CP=CP,meanCP=meanCP,col=col,ylab1=ylab1,ylab2=ylab2,plot.Q=plot.Q,cex.legend=cex.legend,cex.diamonds=cex.diamonds,main=main,plot.h=plot.h,plot.cp=plot.cp)
    
  }
  # End if-plot.CPD
  
  ###
  #### Output
  ###
  
  CP[,1]<-CP[,1]+S 
  if(print.output){
    cat("","\n")
    cat("MFT table",sep="\n"); cat("","\n")	
    cat(paste("Tt = ",Tt," and H = {",paste(H,collapse=", "),"}",sep="")); cat("","\n")	
    
    if(M>parameter$Q)
      {cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(parameter$Q,2),sep=""),sep="\n")}	
    else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(parameter$Q,2),sep=""),sep="\n")}; cat("","\n")
    
    if(perform.CPD){
      cat("CPD was performed: ")
      if(dim(CP)[1]==0){cat("No change points detected")}
      if(dim(CP)[1]==1){cat(paste(dim(CP)[1],"change point detected at "))} 
      if(dim(CP)[1]>1){cat(paste(dim(CP)[1],"change points detected at "))}
      if(dim(CP)[1]>0){cat(paste(sort(CP[,1])),sep=", ")}; cat("","\n"); cat("                   ")
      if(length(meanCP)==1){cat(paste("The estimated mean is",signif(meanCP,2)) )} else{cat("The estimated means are ")
        cat(paste(signif(meanCP,2)),sep=", ")}
    }
    # End-if-perform.CPD	
    if(!perform.CPD){cat("CPD was not performed")}
    
  }
  # End-if-print.output
  
  SWD<-lapply(SWD,function(x) x<-x+S) 
  if(perform.CPD==FALSE){CP<-NA; SWD <- NA; rate <- NA; meanCP<-NA}
  invisible(list(M=M,Q=parameter$Q,CP=CP,meanCP=meanCP,SWD=SWD,S=S,E=E,H=H,sim=sim,alpha=alpha))
  
  }
