#' summary.MFT
#'
#' Summary method for class 'mft'.
#'
#' @param object object of class MFT
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
#' summary(mft)
#' 
#' 
#' @seealso \code{\link{MFT.rate}, \link{MFT.variance}, \link{MFT.mean}, \link{MFT.peaks}, \link{plot.MFT}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Marietta Kirchner, Julia Schiemann, Jochen Roeper, Ralph Neininger and Gaby Schneider (2014).
#' A multiple filter test for the detection of rate changes in renewal processes with varying variance. The Annals of Applied Statistics 8(4): 2027-67
#' <doi:10.1214/14-AOAS782>
#'
#' 
#' @rdname summary.MFT
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


summary.MFT<-function(object,...)
{
  MFT<-object
  S<-MFT$S; E<-MFT$E; CP<-MFT$CP; H<-MFT$H; alpha<-MFT$alpha; sim<-MFT$sim; M<-MFT$M; rescale<-MFT$rescale;
  if (MFT$type=="mean")
  {
    X<-MFT$tech.var$X; Q<-MFT$Q; Tt<-MFT$Tt; means<-MFT$expectation
    cat("","\n")
    cat("MFT table",sep="\n"); cat("","\n")	
    cat(paste("Tt = ",Tt," and H = {",paste(H,collapse=", "),"}",sep="")); cat("","\n")	
    
    if(M>Q)
    {cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(Q,2),sep=""),sep="\n")}	
    else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(Q,2),sep=""),sep="\n")}; cat("","\n")
    
    if(MFT$method=="asymptotic"){cat("(Threshold Q derived by simulation) \n\n")}
    if(MFT$method=="fixed"){cat("(Threshold Q set by user) \n\n")}
    
    if(MFT$perform.CPD){
      cat("CPD was performed: ")
      if(dim(CP)[1]==0){cat("No change points detected")}
      if(dim(CP)[1]==1){cat(paste(dim(CP)[1],"change point detected at "))} 
      if(dim(CP)[1]>1){cat(paste(dim(CP)[1],"change points detected at "))}
      if(dim(CP)[1]>0){cat(paste(sort(CP[,1])),sep=", ")}; cat("","\n"); cat("                   ")
      if(length(means)==1){cat(paste("The estimated mean is",signif(means,2)) )} else{cat("The estimated means are ")
        cat(paste(signif(means,2)),sep=", ")}
    }
    # End-if-perform.CPD	
    if(!MFT$perform.CPD){cat("CPD was not performed")}
   
  }
  if (MFT$type=="rate")
  {
    rescale<-MFT$rescale; Q<-MFT$Q; rate<-MFT$rate; Tt<-MFT$Tt; d<-MFT$d
    
    cat("","\n")
    cat("MFT table",sep="\n"); cat("","\n")	
    cat(paste("Tt = ",Tt,", d = ",d," and H = {",paste(H,collapse=", "),"}",sep="")); cat("","\n")	
    
    if(M>Q){cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(Q,2),sep=""),sep="\n")}	else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(Q,2),sep=""),sep="\n")}#; cat("","\n")
    
    if(MFT$method=="bootstrap"){cat("(Threshold Q derived by bootstrapping) \n\n")}
    if(MFT$method=="asymptotic"){cat("(Threshold Q derived by simulation) \n\n")}
    if(MFT$method=="fixed"){cat("(Threshold Q set by user) \n\n")}
    
    
    if(MFT$perform.CPD){
      cat("CPD was performed: ")
      if(dim(CP)[1]==0){cat("No change points detected")}
      if(dim(CP)[1]==1){cat(paste(dim(CP)[1],"change point detected at "))} 
      if(dim(CP)[1]>1){cat(paste(dim(CP)[1],"change points detected at "))}
      if(dim(CP)[1]>0){cat(paste(sort(CP[,1])),sep=", ")}; cat("","\n")
      if(length(rate)==1){cat(paste("The estimated rate is",signif(rate,2) ))} else{cat("The estimated rates are ")
        cat(paste(signif(rate,2)),sep=", ")}
    }#end-if-perform.CPD	
    if(!MFT$perform.CPD){cat("CPD was not performed")}
    
 }
  if (MFT$type=="variance")
  {
    tech.var<-MFT$tech.var; varQ<-MFT$varQ;   step<-MFT$d; d<-MFT$d; G_ht<-MFT$tech.var$G_ht; Tt<-MFT$Tt; Phi<-MFT$tech.var$Phi
    
    cat("","\n")
    cat("MFT variances table",sep="\n"); cat("","\n")
    cat(paste("Tt = ",Tt, ", d = ",round(step,2)," and H = {",paste(round(H,2),collapse=", "),"}",sep="")); cat("","\n")
    if(M>varQ){cat(paste("Stationarity was rejected: M = ",round(M,2)," > Q = ", Q = round(varQ,2),sep=""),sep="\n")}	else{cat(paste("Stationarity was not rejected: M = ",round(M,2)," < Q = ", Q = round(varQ,2),sep=""),sep="\n")}; cat("","\n")
    
    if(MFT$method=="asymptotic"){cat("(Threshold Q derived by simulation) \n\n")}
    if(MFT$method=="fixed"){cat("(Threshold Q set by user) \n\n")}
    
    if(MFT$perform.CPD){
      cat("CPD was performed: ")
      if(M<=varQ){cat("No change points detected"); cat("","\n")}
      if (M>varQ) {
      if((dim(CP)[1]==1)){cat(paste(dim(CP)[1],"change point detected at "))} 
      if((dim(CP)[1]>=2)){cat(paste(dim(CP)[1],"change points detected at "))}
      if((dim(CP)[1]>0)){cat(paste(round(CP[,1])),sep=", ")}; cat("","\n")
      }
      if(M<varQ){cat(paste("The estimated variance is",signif(MFT$var,3) ));  cat("","\n")} else{cat("The estimated variances are ")
        if(MFT$perform.CPD) {cat(paste(signif(c(CP[,2],CP[dim(CP)[1],3]),3)),sep=", ")}
        cat("","\n")
      }
    }#end-if-perform.CPD	
    if(!MFT$perform.CPD){cat("CPD was not performed")}
     # cat("","\n"); cat(paste("The estimated variance is",signif(MFT$var,3) ))
    
    
    
  }
  if (MFT$type=="peaks")
  {
    rescale<-MFT$rescale; Q<-MFT$Q; rate<-MFT$rate; Tt<-MFT$Tt; d<-MFT$d
    cat("", "\n")
    cat("MFT table", sep = "\n")
    cat("", "\n")
    cat(paste("Tt = ", Tt," and H = {", paste(H, collapse = ", "), "}", sep = ""))
    cat("", "\n")
    if (M > Q) {
      cat(paste("Stationarity was rejected: M = ", round(M,  2), " > Q = ", Q = round(Q, 2), sep = ""), sep = "\n")
    }
    else {
      cat(paste("Stationarity was not rejected: M = ", 
                round(M, 2), " < Q = ", Q = round(Q, 2), sep = ""), sep = "\n")
    }
    if (MFT$method == "bootstrap") {
      cat("(Threshold Q derived by bootstrapping) \n\n")
    }
    if (MFT$method == "asymptotic") {
      cat("(Threshold Q derived by simulation) \n\n")
    }
    if (MFT$method == "fixed") {
      cat("(Threshold Q set by user) \n\n")
    }
    if (MFT$perform.CPD) {
      cat("Peak detection was performed: ")
      if (dim(CP)[1] == 0) {
        cat("No peaks detected")
      }
      if (dim(CP)[1] == 1) {
        cat(paste(dim(CP)[1], "peak detected at "))
      }
      if (dim(CP)[1] > 1) {
        cat(paste(dim(CP)[1], "peaks detected at "))
      }
      if (dim(CP)[1] > 0) {
        cat(paste(sort(CP[, 1])), sep = ", ")
      }
    }
    if (!MFT$perform.CPD) {
      cat("Peak detection was not performed")
    }
  }  
}