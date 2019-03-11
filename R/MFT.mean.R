#' MFT.mean
#'
#' The multiple filter test for mean change detection in time series or sequences of random variables.
#'
#' @param X numeric vector, input sequence of random variables
#' @param S	numeric, start of time interval, default: NULL, if NULL then 1 is chosen
#' @param E numeric, end of time interval, default: NULL, if NULL then length(X) is chosen, needs E > S.
#' @param autoset.H	logical, automatic choice of window size H
#' @param H vector, window set H, all elements must be increasing, the largest element must be =< (T/2). H is automatically set if autoset.H = TRUE
#' @param alpha numeric, in (0,1), significance level
#' @param method either "asymptotic" or "fixed", defines how threshold Q is derived, default: "asymptotic", If "asymptotic": Q is derived by simulation of limit process L (Brownian motion); possible set number of simulations (sim), If "fixed": Q may be set manually (Q)
#' @param sim integer, > 0, No of simulations of limit process (for approximation of Q), default = 10000
#' @param rescale logical, if TRUE statistic G is rescaled to statistic R, default = FALSE 
#' @param Q	numeric, rejection threshold, default: Q is simulated according to sim and alpha.
#' @param perform.CPD logical, if TRUE change point detection algorithm is performed
#' @param print.output logical, if TRUE results are printed to the console	

#' @return invisible
#' \item{M}{test statistic}
#' \item{Q}{rejection threshold}
#' \item{method}{how threshold Q was derived, see 'Arguments' for detailed description}
#' \item{sim}{number of simulations of the limit process (approximation of Q)}
#' \item{rescale}{states whether statistic G is rescaled to R}
#' \item{CP}{set of change points estmated by the multiple filter algorithm, increasingly ordered in time}
#' \item{means}{estimated mean values between adjacent change points}
#' \item{S}{start of time interval}
#' \item{E}{end of time interval}
#' \item{Tt}{length of time interval}
#' \item{H}{window set}
#' \item{alpha}{significance level}
#' \item{perform.CPD}{logical, if TRUE change point detection algorithm was performed}
#' \item{tech.var}{list of technical variables with processes X and G_ht or R_ht}
#' \item{type}{type of MFT which was performed: "mean"}
#' 
#' @seealso \code{\link{plot.MFT}, \link{summary.MFT}, \link{MFT.rate}, \link{MFT.variance}, \link{MFT.peaks}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' 
#' @references 
#' Michael Messer, Stefan Albert and Gaby Schneider (2018). The multiple filter test for change point detection in time
#' series. Metrika <doi:10.1007/s00184-018-0672-1> 
#' 
#' @examples 
#' # Normal distributed sequence with 3 change points of the mean (at n=100, 155, 350)
#' set.seed(50)
#' X1   <- rnorm(400,0,1); X2 <- rnorm(400,3,1); X3 <- rnorm(400,5,1); X4 <- rnorm(600,4.6,1)
#' X    <- c(X1[1:100],X2[101:155],X3[156:350],X4[351:600])
#' mft  <- MFT.mean(X)
#' plot(mft)
#' # Set additional parameters (window set)
#' mft2 <- MFT.mean(X,autoset.H=FALSE,H=c(80,160,240))
#' plot(mft2)
#' 
#' 
#' @rdname MFT.mean
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


MFT.mean <-
function(X,autoset.H=TRUE,S=NULL,E=NULL,H=NULL,alpha=0.05,method="asymptotic",sim=10000,rescale=FALSE,Q=NA,perform.CPD=TRUE,print.output=TRUE){
###
### Set Parameters
###
  
  if(is.null(S)  & is.null(E))	{S <- 1; E <- length(X)}
  if(is.null(S)  & !is.null(E))	{S <- 1}
  if(!is.null(S) & is.null(E))	{E <- length(X)}	
  if(E-S <= 0){stop("Invalid choice of S and E: Need S < E.")}  
  
  
  if(! method %in% c("aymptotic","fixed")){"Invalid choice of method: method must be 'asymptotic' or 'fixed'"}
  if(method == "fixed" & !is.numeric(Q)){stop("Invalid choice of Q: Q must be positive real number")}
  if(method == "fixed" & !(is.numeric(Q) & Q>0)){cat("Warning: non-positive threshold might be inappropriate. Possibly choose Q > 0",sep="\n") }
  if(method != "fixed" & !is.na(Q)){cat("Warning: Q is derived by simulation. In order to set Q manually, choose: method = 'fixed'",sep="\n")}
  
  


# Begin-if-autoset.H
  if(autoset.H){
    Tt  <- length(X[S:E]) 											                        # Length of interval.
    if (Tt<=100) {stop("Automatic choice of windows failed. Time series is too short (need at least 101 random variables). Set H manually.")}
    H <- seq(50,min((Tt-1)/2,200),25)
  }
# End-if-set.H 

# Begin-if-!autoset.H
  if(!autoset.H){
    if(is.null(H)){stop("If autoset.H is FALSE, the vector of window sizes H must be set")}
    Tt  <- length(X[S:E])	
    if(!all(c(diff(c(1,H,(Tt/2))))>0) | !all(H%%1 == 0)){stop("Invalid choice of window set: H needs to be an increasing ordered vector of integers, with min(H) > 1 and max(H) <= Tt/2")}
  }
# End-if-!autoset.H 

  X <- X[S:E]         # Shorten sequence to region of analysis.

  if( !sim%%1==0 | sim <= 0) {stop("Invalid choice of number of simulations: sim must be a positive integer")}
  if(( method=="asymptotic" ) & ( sim < 10000 )){cat("Warning: Number of simulations might be too low",sep="\n")}
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
      if (rescale) {maxmax 	<- apply((t(re)-meanH)/sdH, MARGIN=2,max)}  # Calculates M* sim times.
      if (!rescale){maxmax 	<- apply((t(re)), MARGIN=2,max)}  
      Q 		<- quantile(maxmax,1-alpha)     # Calculates the empirical alpha-Quantile M*.
      list(Q=Q,meanH=meanH,sdH=sdH)         # Writes the threshold Q, and the vectors of means and sd in a list.
    }
# End-sim.maxmax
    
    sim.maxmax(sim=sim,Tt=Tt,H=H,alpha=alpha,rescale=rescale)     # Calls function sim.maxmax.
    
  }
# End-sim.parameter	
  
  parameter <- list()
  if(! method %in% c("aymptotic","fixed")){"Invalid choice of method: method must be 'asymptotic' or 'fixed'"}
  if(rescale | method == "asymptotic") parameter <- sim.parameter(sim=sim,Tt=Tt,H=H,alpha=alpha,rescale=rescale) # Simulates parameter values.
  if(method == "fixed" & !is.numeric(Q)){stop("Invalid choice of Q: Q must be positive real number")}
  if(method == "fixed") {parameter$Q <- Q}
  if(method == "fixed" & !rescale){sim <- NA}
  if( alpha*(1-alpha) <= 0){stop("Invalid choice of significance level: alpha must be in (0,1)")}
  

  ### 
  ### Calculate R_ht
  ###
  
  
  # Calculates G_ht for fixed t and h.
  ght <- function(t,Xp,h){ 
    X <- Xp
    windowRight <- X[(t+1):(t+h)]   # Values in right window.
    windowLeft <- X[(t-h+1):t]
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
  Rh <- function(h,Xp,H,rescale,meanH=NULL,sdH=NULL){
    seq2 <- seq(1+h,length(Xp)-h,by=1)    # All possible t values for window size h.
    G <- vapply(seq2,FUN = ght,Xp=Xp,h=h,numeric(1))	# Calculates all G_ht for window size h.
    if(rescale){return( (abs(G) - meanH[which(H==h)]) / sdH[which(H==h)] )}
    if(!rescale){return(abs(G))}
  }
  # End-Rht

  if (rescale) {R_ht <- lapply(as.matrix(H),FUN=Rh,Xp=X,H=H,rescale=rescale,meanH=parameter$meanH,sdH=parameter$sdH)} # R_ht values for every window size.
  if (!rescale){R_ht <- lapply(as.matrix(H),FUN=Rh,Xp=X,H=H,rescale=rescale)} # R_ht values for every window size.
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
    CPs <- c(1,as.vector(CP[,1]),length(X))           # CP plus start and end
    means <- rep(NA,(length(CPs)-1))         # Vector for means.
    for(i in 1:(length(CPs)-1)){      
      means[i] <- mean(X[CPs[i]:CPs[i+1]])   # Mean for every section. 
    }
    # End-for  
    
    names(means) <- apply(as.matrix(1:length(means)),1,function(index){paste("mu_",index,sep="")})
    
  }
  # End-perform.CPD

  
  ###
  #### Output
  ###
  
  if (perform.CPD) CP[,1]<-CP[,1]+S 
  if(print.output){
    cat("","\n")
    cat("MFT table",sep="\n"); cat("","\n")	
    cat(paste("Tt = ",Tt," and H = {",paste(H,collapse=", "),"}",sep="")); cat("","\n")	
    
    if(M>parameter$Q)
      {cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(parameter$Q,2),sep=""),sep="\n")}	
    else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(parameter$Q,2),sep=""),sep="\n")}; cat("","\n")
    
    if(method=="asymptotic"){cat("(Threshold Q derived by simulation) \n\n")}
    if(method=="fixed"){cat("(Threshold Q set by user) \n\n")}
    
    if(perform.CPD){
      cat("CPD was performed: ")
      if(dim(CP)[1]==0){cat("No change points detected")}
      if(dim(CP)[1]==1){cat(paste(dim(CP)[1],"change point detected at "))} 
      if(dim(CP)[1]>1){cat(paste(dim(CP)[1],"change points detected at "))}
      if(dim(CP)[1]>0){cat(paste(sort(CP[,1])),sep=", ")}; cat("","\n"); cat("                   ")
      if(length(means)==1){cat(paste("The estimated mean is",signif(means,2)) )} else{cat("The estimated means are ")
        cat(paste(signif(means,2)),sep=", ")}
    }
    # End-if-perform.CPD	
    if(!perform.CPD){cat("CPD was not performed")}
    
  }
  # End-if-print.output
  
  #if(perform.CPD) SWD<-lapply(SWD,function(x) x<-x+S) 
  if(perform.CPD==FALSE){CP<-NA; means<-NA}
  if (rescale) tech.var<-list(X=X,R_ht=R_ht,G_ht=NA)
  else tech.var<-list(X=X,R_ht=NA,G_ht=R_ht)
  
  mft<-list(M=M,Q=parameter$Q,method=method,rescale=rescale,CP=CP,expectation=means,S=S,E=E,Tt=Tt,H=H,sim=sim,alpha=alpha,perform.CPD=perform.CPD,tech.var=tech.var,type="mean")
  class(mft) <- "MFT"
  invisible(mft)
}
