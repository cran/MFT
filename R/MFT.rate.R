
#' MFT.rate
#'
#' The multiple filter test for rate change detection in point processes on the line. 
#'
#' @param Phi numeric vector of increasing events, input point process
#' @param m  non-negative integer, dependence parameter: serial corellation up to order m is respected for the estimation of rho
#' @param rescale logical, if TRUE statistic G is rescaled to statistic R 
#' @param cutout  logical, if TRUE for every point, for which the estimated rho becomes negative, the h-neighborhood of G (resp. R) is set to zero. This might only occur, if m > 0 
#' @param S	numeric, start of time interval, default: Smallest multiple of d that lies beyond min(Phi)
#' @param E numeric, end of time interval, default: Smallest multiple of d that lies beyond max(Phi), needs E > S.
#' @param autoset.d_H	logical, automatic choice of window size H and step size d
#' @param d numeric, > 0, step size delta at which processes are evaluated. d is automatically set if autoset.d_H = TRUE
#' @param H vector, window set H, all elements must be increasing ordered multiples of d, the smallest element must be >= d and the largest =< (T/2). H is automatically set if autoset.d_H = TRUE
#' @param alpha numeric, in (0,1), significance level
#' @param sim integer, > 0, No of simulations of limit process (for approximation of Q), default = 10000
#' @param method either "asymptotic", "bootstrap" or "fixed", defines how threshold Q is derived, default: "asymptotic", If "asymptotic": Q is derived by simulation of limit process L (Brownian motion); possible set number of simulations (sim), If "bootstrap": Q is deried by (Block)-Bootstrapping; possibly set number of simulations (sim) and blocksize (blocksize), If "fixed": Q may be set automatically (Q)
#' @param Q	numeric, rejection threshold, default: Q is simulated according to sim and alpha.
#' @param blocksize  NA or integer >= 1, if method == 'bootstrap', blocksize determines the size of blocks (number of life times) for bootstrapping
#                If NA: blocksize is autoset to ceiling((length(n)^(1/4))), while n denotes the number of life times of input
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
#' @param plot.cp logical, indicates if a legend of detected CPs is plotted
#' @param plot.rate logical, indicates if a legend of estimated rates is plotted
#' @param breaks integer, > 0, number of breaks in rate histogram	
#'
#' @return invisible
#' \item{M}{test statistic}
#' \item{Q}{rejection threshold}
#' \item{sim}{number of simulations of the limit process (approximation of Q)}
#' \item{blocksize}{size of blocks (number of life times) for bootstrapping (approximation of Q)}
#' \item{rescale}{states whether statistic G is rescaled to R}
#' \item{m}{order of respected serial correlation (m-dependence)}
#' \item{CP}{set of change points estmated by the multiple filter algorithm, increasingly ordered in time}
#' \item{rate}{estimated mean rates between adjacent change points}
#' \item{SWD}{sets of change points estimated from preprocessing single window detections}
#' \item{S}{start of time interval}
#' \item{E}{end of time interval}
#' \item{H}{window set}
#' \item{d}{step size delta at which processes were evaluated}
#' \item{alpha}{significance level}
#' \item{cutout}{states whether cutout was used (see arguments)}
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
#' MFT.rate(Phi,breaks=30)
#' 
#' 
#' @seealso \code{\link{MFT.m_est}, \link{MFT.variance}, \link{MFT.mean}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Marietta Kirchner, Julia Schiemann, Jochen Roeper, Ralph Neininger and Gaby Schneider (2014).
#' A multiple filter test for the detection of rate changes in renewal processes with varying variance. The Annals of Applied Statistics 8(4): 2027-67
#' <doi:10.1214/14-AOAS782>
#'
#' Michael Messer, Kaue M. Costa, Jochen Roeper and Gaby Schneider (2017).
#' Multi-scale detection of rate changes in spike trains with weak dependencies. Journal of Computational Neuroscience, 42 (2), 187-201.
#' <doi:10.1007/s10827-016-0635-3>
#' 
#' @rdname MFT.rate
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


MFT.rate <-
function(Phi,m=0,rescale=TRUE,cutout=TRUE,autoset.d_H=TRUE,S=NULL,E=NULL,H=NULL,d=NULL,alpha=0.05,sim=10000,method="asymptotic",Q=NA,blocksize=NA,perform.CPD=TRUE,plot.CPD=TRUE,print.output=TRUE,col=NULL,ylab1=NULL,ylab2=NULL,cex.legend=1.2,cex.diamonds=1.4,main=TRUE,plot.Q=TRUE,plot.M=TRUE,plot.h=TRUE,plot.rate=FALSE,plot.cp=FALSE,breaks=NULL){
  
  ###
  ### Set Parameters
  ###
  
  if(is.null(S)  & is.null(E))	{S <- min(Phi); E <- max(Phi)}
  if(is.null(S)  & !is.null(E))	{S <- min(Phi)}
  if(!is.null(S) & is.null(E))	{E <- max(Phi)}	
  if(E-S <= 0){stop("Invalid choice of S and E: Need S < E.")}
  
  if(autoset.d_H){
    H   <- signif(100 / (length(Phi[Phi>S & Phi<=E])/(E-S)),digits=1)  	# rounds to first significant
    d   <- (H / 20)			 											          # step size	
    S   <- floor(S/d)*d													        # floor S to next d
    Tt  <- ceiling((E-S)/d)*d 											    # ceil Tt to next d
    E   <- S + Tt														            # set E
    if(2*H > Tt){stop("Can not set parameter: intensity is too low / Phi too short")}
    if(Tt/5 > H){H <- seq(H,Tt/5,H)}			 							# Window choice	
  }#end-if-set.H_and_d
  
  if(!autoset.d_H){
    if(is.null(d) | is.null(H)){stop("if autoset.d_H is FALSE, the step size d and vector of window sizes H must be set")}
    if(d<=0){stop("Invalid choice of step size d: Need d > 0")}	
    S   <- floor(S/d)*d
    Tt  <- ceiling((E-S)/d)*d	
    E	  <- S + Tt
    if(!all(c(diff(c(d,H,(Tt/2))))>0) | !all(H%%d == 0)){stop("Invalid choice of window set: H needs to be an increasing ordered vector of multiples of d, with min(H) > d and max(H) <= Tt/2")}
  }#end-if-!set.H_and_d 
  
  #print(paste("S =",S,"min(Phi)/d)*d =",min(Phi)/d)*d,"ceiling(max(Phi)/d)*d =",ceiling(max(Phi)/d)*d,"E =",E)
  
  if(S < floor(min(Phi)/d)*d | E > ceiling(max(Phi)/d)*d){cat("Warning: Interval [S,E] does not suit time horizon of Phi",sep="\n") }
  
  Phi <- Phi[Phi > S & Phi <= E]
  
  if(!is.logical(rescale) | is.na(rescale)){stop("Invalid choice of rescaling parameter: rescale must be logical")}
  
  if(! method %in% c("aymptotic","bootstrap","fixed")){"Invalid choice of method: method must be 'asymptotic', 'bootstrap' or 'fixed'"}
  if(method == "fixed" & !is.numeric(Q)){stop("Invalid choice of Q: Q must be positive real number")}
  if(method == "fixed" & !(is.numeric(Q) & Q>0)){cat("Warning: non-positive threshold might be inappropriate. Possibly choose Q > 0",sep="\n") }
  if(method != "fixed" & !is.na(Q)){cat("Warning: Q is derived by simulation / bootstrapping. In order to set Q manually, choose: method = 'fixed'",sep="\n")}
  
  if(!is.logical(cutout)){stop("Invalid cutout option: cutout must be logical")}
  
  if(method == "bootstrap" & rescale == TRUE){stop("rescaling not available for bootstrapping: set rescale to FALSE")}
  if(method == "bootstrap" & !(is.na(blocksize) | is.numeric(blocksize)) ){stop("Invalid choice of blocksize for Bootstapping: blocksize must be positive integer")}
  if(method == "bootstrap" & !(is.na(blocksize) | (is.numeric(blocksize) & blocksize%%1 == 0 & blocksize>=1))){stop("Invalid choice of blocksize for bootstapping: blocksize must be positive integer")}
  if(method == "bootstrap" & sim > 1000){cat("Warning: High comupational time. Possibly reduce the number of bootstrap-simulations: sim",sep="\n")}
  
  if(!is.numeric(m)){stop("Invalid choice of dependence parameter: m must be non-negative integer")}
  if(m%%1!=0 | m<0){stop("Invalid choice of dependence parameter: m must be non-negative integer")}  
  
  if(!is.numeric(sim)){stop("Invalid choice of number of simulations: sim must be a positive integer")}
  if(!sim%%1==0 | sim <= 0) {stop("Invalid choice of number of simulations: sim must be a positive integer")}
  
  if((rescale | method == "asymptotic") & sim < 10000){cat("Warning: Number of simulations might be too low",sep="\n")}
  if( alpha*(1-alpha) <= 0){stop("Invalid choice of significance level: alpha must be in (0,1)")}
  
  ###
  ### sim.parameter 
  ###
  
  parameter <- list(Q=Q) # eventually updated depending on method
  
  if(rescale | method == "asymptotic"){# simulation of rescaling constants (if rescale == TRUE) or threshold Q (if method == "asymptotic")
    
    sim.parameter <- function(sim=sim,Tt=Tt,H=H,alpha=alpha,d=d,rescale=rescale){
      
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
      
      sim.maxmax <- function(sim=sim,Tt=Tt,H=H,alpha=alpha,d=d,rescale=rescale){ # Simulates Brownian motions and calculates mean(max|L_ht|), sd(max|L_ht|) and Q
        re 		<- matrix(replicate(sim,sim.maxH(Tt=Tt,H=H,d=d)),nrow=sim,byrow=TRUE)	
        meanH 	<- apply(re,MARGIN=2,mean)
        sdH   	<- apply(re,MARGIN=2,sd) 
        if(rescale){maxmax 	<- apply((t(re)-meanH)/sdH, MARGIN=2,max)}  
        if(!rescale){maxmax 	<- apply((t(re)), MARGIN=2,max)}  
        if(method != "fixed"){Q 		<- quantile(maxmax,1-alpha)}
        list(Q=Q,meanH=meanH,sdH=sdH)
      }#end-sim.maxmax
      
      sim.maxmax(sim=sim,Tt=Tt,H=H,alpha=alpha,d=d,rescale=rescale)
      
    }#end-sim.parameter	
    
    parameter <- sim.parameter(sim=sim,Tt=Tt,H=H,alpha=alpha,d=d,rescale=rescale) # Simulates parameter values
    
  }# end-if(rescale) 
  
  if(method == "fixed" & !rescale){sim <- NA}
  
  ### 
  ### Calculate G_ht (fixed h and t)
  ###
  
  numerator <- function(t,Phi,h){ # Calculates Numerator of G_ht for fixed t and h
    v 		<-  hist(Phi[Phi > t - h & Phi <= t + h],breaks=c(t-h,t,t+h),plot=FALSE)$counts 
    return(v[2] - v[1])
  }#end-ght
  
  
  rho_sq <- function(lt,m){ # Calculates hat_rho^2
    rho_squared <- var(lt) # NA, if length(lt) = 0
    if(length(lt)<=m){
      rho_squared <- NA  
    }
    if(m>0 & length(lt)>m){
      rho_l <- apply(as.matrix(c(1:m)),MARGIN=1,
                     function(lag,lt){(1/(length(lt)-lag)) * sum(lt[1:(length(lt)-lag)] * lt[(1+lag):length(lt)])}
                     ,lt) - mean(lt)^2
      rho_squared <- rho_squared + 2*sum(rho_l)
    }#end-if-m    
    rho_squared
  }#end-rho_sq
  
  
  denominator <- function(t,Phi,h){ # Calculates Denominator of G_ht for fixed t and h
    lt1   	<-  diff(Phi[Phi > t - h & Phi <= t])
    lt2   	<-  diff(Phi[Phi > t & Phi <= t + h])
    va1    	<- rho_sq(lt1,m)*h/mean(lt1)^3
    va2	   	<- rho_sq(lt2,m)*h/mean(lt2)^3
    if (is.na(va1) | is.na(va2) | va1 <= 0 | va2 <= 0) {erg <- 0}
    else {erg <- sqrt(va1 + va2)}
    return(erg)
  }#end-ght
  
  ####
  #### Calculate R_ht 
  ####
  
  Rh <- function(h,Phi,m,S,E,H,d,cutout,rescale,meanH,sdH){ # Calculate R_ht for all h and t
    numerator_vector    <- vapply(seq((S+h),(E-h),d),FUN=numerator,Phi=Phi,h=h,numeric(1))	
    denominator_vector  <- vapply(seq((S+h),(E-h),d),FUN=denominator,Phi=Phi,h=h,numeric(1))	
    G                   <- numerator_vector / denominator_vector
    G                   <- ifelse(is.nan(G) | G==Inf |G==-Inf,0,G)
    if(cutout){# Set all values of G to zero, that lie in an h-neighborhood of a negative rho 
      tau     <- seq((S+h),(E-h),d)
      sub_tau <- tau[denominator_vector<=0] # time points for which denominator is negative
      if(length(sub_tau)>0){
        all_timepoints  <- as.vector(apply(cbind(sub_tau-h,sub_tau+h),MARGIN=1,function(mat_row){seq(mat_row[1],mat_row[2],d)})) 
        settozero       <- intersect(union(all_timepoints,all_timepoints),tau)  # all time points to set to zero
        G               <- ifelse(tau%in%settozero,0,G) # Set corresponding G values to zero
      }#end-if-length-sub_tau>0 
    }#-if-cutout
    if(rescale){return( (abs(G) - meanH[which(H==h)]) / sdH[which(H==h)] )}
    if(!rescale){return(abs(G))}
  }#end-R
  
  R_ht <- lapply(as.matrix(H),FUN=Rh,Phi=Phi,m=m,S=S,E=E,H=H,d=d,rescale=rescale,cutout=cutout,meanH=parameter$meanH,sdH=parameter$sdH)
  M 	 <- max(vapply(R_ht,max,numeric(1))); names(M) <- "test statistic M"		
  
  ### 
  ### Eventually derive Q via Bootstrapping
  ###
  
  if(method == "bootstrap"){
    
    fun_M_boot <- function(Phi=Phi,blocksize=blocksize,S=S,E=E,H=H,d=d,rescale=rescale,cutout=cutout){# Evaluates one Bootstrap Version of M
      lt                <- diff(Phi)
      startblock        <- sample(1:(length(lt)-blocksize+1),ceiling(length(lt) / blocksize),replace=TRUE) 
      index_lt_boot     <- as.vector(apply(as.matrix(startblock),MARGIN=1,function(startvalue,blocksize){startvalue:(startvalue+blocksize-1)},blocksize))[1:length(lt)]
      Phi_boot          <- cumsum(lt[index_lt_boot]) # Bootstrapped process
      R_ht_boot         <- lapply(as.matrix(H),FUN=Rh,Phi=Phi_boot,S=S,E=E,H=H,d=d,rescale=rescale,cutout=cutout)
      max(vapply(R_ht_boot,max,numeric(1)))	# Bootstrapped statistic M
    }#end-fun_M_boot
    
    if(is.na(blocksize)){blocksize <- ceiling(length(Phi)^(1/4))} # autoset blocksize
    
    M_boot      <- replicate(sim,fun_M_boot(Phi=Phi,blocksize=blocksize,S=S,E=E,H=H,d=d,rescale=rescale,cutout=cutout))
    parameter$Q <- quantile(M_boot,1-alpha)
    
  }#end-if-method=="bootstrap"
  
  if(method != "bootstrap"){blocksize <- NA}
  
  ###
  ### Perform CPD
  ###	
  
  if(perform.CPD){
    
    # Perform SWD for fixed h
    SWD <- function(R,Tt,d,Q,S){ 
      h  <- (Tt - (length(R)-1)*d) / 2
      CP <- rep(NA,round(Tt/d,0)); j <- 1 		  
      while ( any(R > Q) ) {
        c <- which( R == max(R) )[1]  
        CP[j] <- (c-1)*d + h + S
        left 		  <- max (1, c-round(h/d,0)+1)  
        right  		  <- min (round((Tt-2*h)/d,0)+1, c+round(h/d,0)-1)
        R[left:right] <- rep(Q-1,length(left:right))
        j <- j + 1
      } #end-while
      if( all(is.na(CP)) ){return(numeric(0))} 
      else{return(as.vector(na.omit(CP)))}
    }#end-SWD	
    
    SWD <- lapply(R_ht,FUN=SWD,Tt=Tt,d=d,Q=parameter$Q,S=S) # SWD for all h in H
    
    # Combine SWD results
    delete <- function(element,h2,c1){ # Is element c2 deleted? (yes, if ball is hit)
      hit    <- any( ifelse( element - h2 < c1  &  c1 < element + h2, TRUE, FALSE) )  
      return(hit)
    }#end-delete
    
    outside <- function(c1,c2,h2){ # outside is logical vector of length c2. TRUE -> Component deleted. FALSE -> Component remains.
      tot_hit <- apply(as.matrix(c2),MARGIN=1,FUN=delete,h2=h2,c1=c1) 
      return(tot_hit)
    }#end-function
    
    # With "outside" run through all "h-areas" and delete "overlap-CPs"
    c		  <- SWD[[1]] 					# First, accept CPs detected with smallest h
    h_val <- rep(H[1],length(SWD[[1]])) 	# corresponding h value
    
    if(length(H)>1){
      for (k in 2:(length(H))){ 
        if(!is.na(SWD[[k]][1]) ){
          put_out <- outside(c,SWD[[k]],H[k]) # deleted if TRUE
          c 	 	<- c(c,SWD[[k]][put_out == FALSE])
          h_val	<- c(h_val,rep(H[k],length(c)-length(h_val))) 
        }#end-if
      }#end-for
    }#end-if-length(H)
    
    CP <- cbind(c,h_val); if(dim(CP)[1]>1){CP <- CP[order(CP[,1]),]}
    colnames(CP) <- c("changepoint","h_value"); if(dim(CP)[1] > 0) {rownames(CP) <- as.character(1:dim(CP)[1])}
    
    rate    	<- hist(Phi,breaks=c(S,CP[,1],E),plot=FALSE)$counts / diff(c(S,CP[,1],E)) # Calculate mean rate
    names(rate) <- apply(as.matrix(1:length(rate)),1,function(index){paste("lambda_",index,sep="")})
    
  }#end-perform.CPD
  
  ###
  #### plot of scenario
  ###
  
  if(plot.CPD & !perform.CPD){cat("Can not plot scenario. Need perform.CPD == TRUE",sep="\n")}
  
  if(plot.CPD & perform.CPD){
    
    plot.cpd <- function(Phi,R_ht,S,E,d,H,Q,M,CP,rate,col,ylab1,ylab2,plot.Q,cex.legend,cex.diamonds,main,plot.h,plot.rate,plot.cp,breaks){
      
      if(is.null(col)){col <- rainbow(length(H),start=2/3,end=0)} else{if(length(col) == 1 & col[1] == "gray"){col <- gray.colors(length(H),0.8,0)}
        else{col <- col; if(length(col)!=length(H)){stop("Argument col must be either NULL, ''gray'', or a vector of colors of length of H")}}
      } # Vector of colors for different windows
      col.cp <- apply(as.matrix(CP[,2]),MARGIN=1, function(x){col[which(x == H)]}) # Vector that codes CP <-> col
      
      # Graphic 1: R_ht	
      
      if(rescale){mi <- min (  c(min( sapply(R_ht, min ) , -7 ) )) - 1 }else{mi <- 0}
      ma <- max (  c(max( sapply(R_ht, max ) ,  Q ) )) + 1  # Eventuelle exchange Q with 12
      
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
      
      if(dim(CP)[1] > 0){ # CPs 
        for(i in 1:(dim(CP)[1]) ){
          points(CP[i,1],R_ht[[which(CP[i,2] == H)]][((CP[i,1]-S-(CP[i,2]-d))/d)],col= 1,cex=cex.diamonds)
          points(CP[i,1],mi,col=col[which(CP[i,2] == H)],pch=18,cex=cex.diamonds+0.4)
          points(CP[i,1],mi,col=1,pch=5,cex=cex.diamonds)
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
      
      if(plot.rate != FALSE){ # Legend for rate
        legend(x=plot.rate,inset=c(0,0.02),legend=as.character(round(rate,2)),pch=19,cex=cex.legend,bty="n",title=expression(paste(hat(lambda)[i] %~~% "")),title.adj=0.2)
      }#end-if
      
      if(dim(CP)[1] > 0){ # Legend for CPs	
        if(plot.cp != FALSE){
          legend(x=plot.cp,inset=c(0,0.02),legend=as.character(CP[,1]),col=col.cp,pch=18,cex=cex.legend,bty="n",text.col=col.cp,title=expression(paste(hat(c)[i]," = ")),title.col=1,pt.cex=cex.legend+0.4 )
          legend(x=plot.cp,inset=c(0,0.02),legend=as.character(CP[,1]),col=1,pch=5,cex=cex.legend,bty="n",text.col=0,title=expression(paste(hat(c)[i]," = ")),title.col=0,pt.cex=cex.legend ) 
        }#end-if
      }#end-if
      
      # Graphic 2: Rate histogram
      par(mar=c(2,4,0.5,0.5))
      if (is.null(breaks)){breaks <- seq(S,E,length=round(Tt/H[1]))}
      histo 	  		<-  hist(Phi,breaks=breaks,plot=FALSE)
      bin       		<-  diff(histo$breaks[1:2])
      histo$counts  	<-  histo$counts / bin
      plot(histo,xlab="",ylab=ylab2,main="",axes=FALSE,xlim=c(S,E+(1/10)*(E-S)),ylim=c(0,1.3*max(histo$counts)),cex.lab=cex.legend) 
      axis(1,at=c(S,CP[,1],E),cex.axis=cex.legend); axis(2,cex.axis=cex.legend)
      lines(c(S,CP[,1],E),c(rate,rate[length(rate)]),type="s",col=2,lwd=2)
      if(length(CP[,1] > 0)){
        for(i in 1:length(CP[,1])){
          arrows(CP[i,1],2*max(histo$counts),CP[i,1],max(rate[i:(i+1)])*1.15,col=col.cp[i],length=0.1,lwd=2)	
        }#end-for
      }#end-if
      
    }#end-plot.cpd
    
    if(is.null(ylab1) & rescale){ylab1  <- expression(paste(R[list(h,t)]))}
    if(is.null(ylab1) & !rescale){ylab1  <- expression(paste("|",G[list(h,t)],"|",sep =""))}
    if(is.null(ylab2)){ylab2  <- "rate"}
    
    # Perform plot:
    plot.cpd(Phi=Phi,R_ht=R_ht,S=S,E=E,d=d,H=H,M=M,Q=parameter$Q,CP=CP,rate=rate,col=col,ylab1=ylab1,ylab2=ylab2,plot.Q=plot.Q,cex.legend=cex.legend,cex.diamonds=cex.diamonds,main=main,plot.h=plot.h,plot.rate=plot.rate,plot.cp=plot.cp,breaks=breaks)
    
  }#end if-plot.CPD
  
  ###
  #### Output
  ###
  
  if(print.output){
    cat("","\n")
    cat("MFT table",sep="\n"); cat("","\n")	
    cat(paste("Tt = ",Tt,", d = ",d,", m = ",m," and H = {",paste(H,collapse=", "),"}",sep="")); cat("","\n")	
    
    if(M>parameter$Q){cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(parameter$Q,2),sep=""),sep="\n")}	else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(parameter$Q,2),sep=""),sep="\n")}#; cat("","\n")
    
    if(method=="bootstrap"){cat("(Threshold Q derived by bootstrapping) \n\n")}
    if(method=="asymptotic"){cat("(Threshold Q derived by simulation) \n\n")}
    if(method=="fixed"){cat("(Threshold Q set by user) \n\n")}
    
    
    if(perform.CPD){
      cat("CPD was performed: ")
      if(dim(CP)[1]==0){cat("No change points detected")}
      if(dim(CP)[1]==1){cat(paste(dim(CP)[1],"change point detected at "))} 
      if(dim(CP)[1]>1){cat(paste(dim(CP)[1],"change points detected at "))}
      if(dim(CP)[1]>0){cat(paste(sort(CP[,1])),sep=", ")}; cat("","\n")
      if(length(rate)==1){cat(paste("The estimated rate is",signif(rate,2) ))} else{cat("The estimated rates are ")
        cat(paste(signif(rate,2)),sep=", ")}
    }#end-if-perform.CPD	
    if(!perform.CPD){cat("CPD was not performed")}
    
  }#-end-if-print.output
  
  if(perform.CPD==FALSE){CP<-NA; SWD <- NA; rate <- NA}
  invisible(list(M=M,Q=parameter$Q,method=method,sim=sim,blocksize=blocksize,rescale=rescale,m=m,CP=CP,rate=rate,SWD=SWD,S=S,E=E,H=H,d=d,alpha=alpha,cutout=cutout))
  
}
