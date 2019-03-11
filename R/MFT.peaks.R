#' MFT.peaks
#'
#' The multiple filter test for peak detection in time series or sequences of random variables
#'
#' @param x numeric vector, input sequence of random variables
#' @param autoset.H	logical, automatic choice of window size H
#' @param S	numeric, start of time interval, default: NULL, if NULL then 1 is chosen
#' @param E numeric, end of time interval, default: NULL, if NULL then length(X) is chosen, needs E > S
#' @param H vector, window set H, the smallest element must >= 3 be and the largest =< (T/2). H is automatically set if autoset.H = TRUE
#' @param alpha numeric, in (0,1), significance level
#' @param method either "asymptotic", "bootstrap" or "fixed", defines how threshold Q is derived, default: "asymptotic", If "asymptotic": Q is derived by simulation of limit process L (Gaussian process); possible set number of simulations (sim), If "bootstrap": Q is derived by (Block)-Bootstrapping; possibly set number of simulations (sim) and blocksize (blocksize), If "fixed": Q may be set manually (Q)
#' @param sim integer, > 0, No of simulations of limit process (for approximation of Q), default = 10000
#' @param Q	numeric, rejection threshold, default: Q is simulated according to sim and alpha
#' @param blocksize  NA or integer >= 1, if method == 'bootstrap', blocksize determines the size of blocks (number of life times) for bootstrapping
#                If NA: blocksize is autoset to ceiling((length(n)^(1/4))), while n denotes the number of life times of input
#' @param two.sided logical, if TRUE a two sided test is performed and also negative peaks are considered in peak detection
#' @param perform.CPD logical, if TRUE change point detection algorithm is performed
#' @param print.output logical, if TRUE results are printed to the console	




#' @return invisible
#' \item{M}{test statistic}
#' \item{Q}{rejection threshold}
#' \item{method}{how threshold Q was derived, see 'Arguments' for detailed description}
#' \item{sim}{number of simulations of the limit process (approximation of Q)}
#' \item{blocksize}{size of blocks (number of life times) for bootstrapping (approximation of Q)}
#' \item{CP}{set of change points estmated by the multiple filter algorithm, increasingly ordered in time}
#' \item{S}{start of time interval}
#' \item{E}{end of time interval}
#' \item{Tt}{length of time interval}
#' \item{H}{window set}
#' \item{alpha}{significance level}
#' \item{two.sided}{logigal, if TRUE also negative peaks are considered}
#' \item{perform.CPD}{logical, if TRUE change point detection algorithm was performed}
#' \item{tech.var}{list of technical variables with processes x and D_ht}
#' \item{type}{type of MFT which was performed: "peaks"}
#'
#' @examples 
#' # Normal distributed sequence with 2 peaks
#' set.seed(12)
#' m <- c(rep(0,30),seq(0,3,length.out = 100),seq(3,0,length.out = 80),rep(0,10),
#'        seq(0,6,length.out = 50),seq(6,0,length.out = 50),rep(0,30))
#' x <- rnorm(length(m),m)
#' mft <- MFT.peaks(x)
#' plot(mft)
#' # Set additional parameters (window set)
#' mft <- MFT.peaks(x,autoset.H = FALSE, H =c(30,60,90))
#' plot(mft)
#' 
#' @seealso \code{\link{MFT.filterdata}, \link{plot.MFT}, \link{summary.MFT}, \link{MFT.mean}, \link{MFT.rate}, \link{MFT.variance}}
#' @author Michael Messer, Stefan Albert, Solveig Plomer and Gaby Schneider
#' @references 
#' Michael Messer, Hendrik Backhaus, Albrecht Stroh and Gaby Schneider (2019+). Peak detection in times series 
#' 
#' @rdname MFT.peaks
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


MFT.peaks <- function(x,autoset.H=TRUE,S=NULL,E=NULL,H=NULL,alpha = 0.05,method = "asymptotic", sim = 10000, Q = NA, blocksize = NA,two.sided=FALSE,perform.CPD=TRUE,print.output=TRUE)
{  
  if (is.null(S) & is.null(E)) {
    S <- 1
    E <- length(x)
  }
  if (is.null(S) & !is.null(E)) {
    S <- 1
  }
  if (!is.null(S) & is.null(E)) {
    E <- length(x)
  }
  if (E - S <= 0 | S %% 1 != 0 | E %% 1 != 0) {
    stop("Invalid choice of S and E: need to be integers with S < E")
  }
  
  x  <- x[S:E]
  Tt <- length(x)
  
  if (!method %in% c("asymptotic", "fixed", "bootstrap")) {
    stop("Invalid choice of method: method must be 'asymptotic', 'fixed' or 'bootstrap'")
  }
  if (method == "fixed" & !is.numeric(Q)) {
    stop("Invalid choice of Q: Q must be positive real number")
  }
  
  if (autoset.H) {
    if(Tt/2<50){stop("Automatic choice of windows failed. Time series is too short. Set H manually")}
    H <- seq(50,min(max(Tt/2),200),25)
  }
  
  if (!autoset.H) {
    if (is.null(H)) {
      stop("If autoset.H is FALSE, the vector of window sizes H must be set")
    }
    if (!all(c(diff(c(2, H, (Tt/2)))) > 0) | !all(H%%1 == 
                                                  0)) {
      stop("Invalid choice of window set: H needs to be an increasing ordered vector of integers, with min(H) >= 3 and max(H) <= Tt/2")
    }
  }
  if (!sim%%1 == 0 | sim <= 0) {
    stop("Invalid choice of number of simulations: sim must be a positive integer")
  }
  if (sim < 10000) {
    cat("Warning: Number of simulations might be too low",sep = "\n")
  }
  
  if (!method %in% c("asymptotic", "bootstrap", "fixed")) {
    "Invalid choice of method: method must be 'asymptotic', 'bootstrap' or 'fixed'"
  }
  if (method == "fixed" & !is.numeric(Q)) {
    stop("Invalid choice of Q: Q must be positive real number")
  }
  if (method == "fixed" & !(is.numeric(Q) & Q > 0)) {
    cat("Warning: non-positive threshold might be inappropriate. Possibly choose Q > 0", 
        sep = "\n")
  }
  if (method != "fixed" & !is.na(Q)) {
    cat("Warning: Q is derived by simulation / bootstrapping. In order to set Q manually, choose: method = 'fixed'", 
        sep = "\n")
  }
  
  if (method == "bootstrap" & (is.na(blocksize) | !is.numeric(blocksize))) {
    stop("Invalid choice of blocksize for Bootstapping: blocksize must be positive integer")
  }
  
  if (method == "bootstrap" & (is.na(blocksize) | !(is.numeric(blocksize) & blocksize%%1 == 0 & blocksize >= 1))) {
    stop("Invalid choice of blocksize for bootstapping: blocksize must be positive integer")
  }
  
  if (method == "bootstrap" & sim > 1000) {
    cat("Warning: High comupational time. Possibly reduce the number of bootstrap-simulations: sim", 
        sep = "\n")
  }
  
  
  if (alpha * (1 - alpha) <= 0) {
    stop("Invalid choice of significance level: alpha must be in (0,1)")
  }
  
  
  
  if (method == "asymptotic") {
    sim.parameter <- function(sim = sim, Tt = Tt, H = H, 
                              alpha = alpha,two.sided=two.sided) {
      maxh <- function(h = h, GP1 = GP1, GP2 = GP2, Tt = Tt, two.sided = two.sided) {
        GP1l   <-  GP1[1:(round((Tt - 2 * h), 0) + 1)]
        GP1m   <-  GP1[(h + 1):(round((Tt - h), 0) + 1)]
        GP1r   <-  GP1[(2 * h + 1):(round(Tt, 0) + 1)]
        GP2l   <-  GP2[1:(round((Tt - 2 * h), 0) + 1)]
        GP2m   <-  GP2[(h + 1):(round((Tt - h), 0) + 1)]
        GP2r   <-  GP2[(2 * h + 1):(round(Tt, 0) + 1)]
        tau    <-  h : (Tt-h)
        left  <- (GP1m - GP1l) - (tau-h+h/2) * (GP2m - GP2l)
        right <- (GP1r - GP1m) - (tau+h/2)   * (GP2r - GP2m)
        L <- sqrt(6/h^3) * (left - right)
        if(two.sided){L <- abs(L)}
        return(max(L))
      }
      
      sim.maxH <- function(Tt = Tt, H = H, two.sided = two.sided) {
        # Simulate Gaussian process
        increments <- rnorm(round(Tt, 0) , sd = sqrt(1))
        k <- seq(0,Tt,1)
        GP1 <- (cumsum(k*c(0,increments)))
        GP2 <- cumsum(c(0,increments))
        return(vapply(as.matrix(H), FUN = maxh, GP1 = GP1, GP2 = GP2, two.sided = two.sided,
                      Tt = Tt, numeric(length(1))))
      }
      
      sim.maxmax <- function(sim = sim, Tt = Tt, H = H, 
                             alpha = alpha, two.sided = two.sided) {
        re <- matrix(replicate(sim, sim.maxH(Tt = Tt, 
                                             H = H, two.sided = two.sided)), nrow = sim, byrow = TRUE)
        maxmax <- apply((t(re)), MARGIN = 2, max)
        Q <- quantile(maxmax, 1 - alpha)
        return(Q)
      }
      
      sim.maxmax(sim = sim, Tt = Tt, H = H, alpha = alpha, two.sided = two.sided)
    }
    Q <- sim.parameter(sim = sim, Tt = Tt, H = H, alpha = alpha, two.sided = two.sided)
  }
  
  # 1. Calculation of all G_h,t and test
  G_function_slope <- function(t,x,h)
  {
    xl <- x[(t-h+1):t]; xr <- x[(t+1):(t+h)];
    
    beta_r   <- (12/(h^2-1)) *  ( (1/h)*sum( ((t+1):(t+h))*xr) - (t+(h+1)/2)*mean(xr) )
    beta_l   <- (12/(h^2-1)) *  ( (1/h)*sum( ((t-h+1):t)*xl) - (t-h+(h+1)/2)*mean(xl) )
    
    alpha_r <- mean(xr) - beta_r * (t+(h+1)/2)
    alpha_l <- mean(xl) - beta_l * (t-h+(h+1)/2)
    
    sigma_sq_r <- 1/(h-2) * sum( (xr-(alpha_r+beta_r*((t+1):(t+h))))^2)
    sigma_sq_l <- 1/(h-2) * sum( (xl-(alpha_l+beta_l*((t-h+1):t)))^2)
    
    den_squared <- (12/(h*(h^2-1))) * (sigma_sq_r + sigma_sq_l)
    
    (beta_l - beta_r) / sqrt(den_squared)
  }
  
  G_h_fun <- function(h,x)
  {
    tau_h <- h:(length(x)-h)
    apply(as.matrix(tau_h),MARGIN=1,FUN=G_function_slope,x=x,h=h)
  }
  
  G_list <- lapply(as.matrix(H), FUN = G_h_fun, x = x)
  Mh <- vapply(G_list, max, numeric(1))
  M <- max(Mh)
  names(M) <- "Test statistic M"
  
  if (method == "bootstrap") {
    fun_M_boot <- function(x = x, blocksize = blocksize, H = H) {
      startblock <- sample(1:(length(x) - blocksize + 
                                1), ceiling(length(x)/blocksize), replace = TRUE)
      index_x_boot <- as.vector(apply(as.matrix(startblock), 
                                      MARGIN = 1, function(startvalue, blocksize) {
                                        startvalue:(startvalue + blocksize - 1)
                                      }, blocksize))[1:length(x)]
      x_boot    <- x[index_x_boot]
      G_ht_boot <- lapply(as.matrix(H), FUN = G_h_fun, x = x_boot)
      max(vapply(G_ht_boot, max, numeric(1)))
    }
    M_boot <- replicate(sim, fun_M_boot(x = x, blocksize = blocksize,H = H))
    Q <- quantile(M_boot, 1 - ifelse(two.sided,alpha/2,alpha))
  }
  
  # 2. Single filter algorithm 
  SFA <- function(G,Tt=Tt,Q=Q,two.sided=two.sided)
  {
    h       <- (Tt-length(G)+1)/2
    hat_cp  <- numeric(0)   # Getting updated
    Gnew    <- G            # Getting updated
    total_deleted_timeg <- c(); #total_deleted_time <- c()
    if(two.sided){Gnew <- abs(Gnew)}
    
    while(max(Gnew) > Q)
    {
      # Calculate stats (new CP and deleted time) for 'current CP'
      hat_cg <- which.max(Gnew); hat_c <- hat_cg+h-1; # maximizer
      deleted_timeg <- (hat_cg-(h-1)) : (hat_cg+(h-1))
      deleted_timeg <- deleted_timeg[deleted_timeg>=1 & deleted_timeg<=length(G)] # restrict to timedomain of G
      # update upper quantities
      total_deleted_timeg <- c(total_deleted_timeg,deleted_timeg); # total time to be deleted
      Gnew[deleted_timeg] <- 0
      hat_cp <- c(hat_cp,hat_c)
    } # end-while
    hat_cp
  } # end-SFA
  
  SWD <- lapply(G_list, FUN = SFA, Tt = Tt, Q = Q,two.sided=two.sided); SWD   
  
  # 3. Multi filter algorithm 
  delete <- function(element, h2, c1) 
  {
    hit <- any(ifelse(element - h2 < c1 & c1 < element + 
                        h2, TRUE, FALSE))
    return(hit)
  }
  
  outside <- function(c1, c2, h2) 
  {
    tot_hit <- apply(as.matrix(c2), MARGIN = 1, FUN = delete, 
                     h2 = h2, c1 = c1)
    return(tot_hit)
  }
  
  c <- SWD[[1]]
  h_val <- rep(H[1], length(SWD[[1]]))
  if (length(H) > 1) 
  {
    for (k in 2:(length(H))) 
    {
      if (!is.na(SWD[[k]][1])) 
      {
        put_out <- outside(c, SWD[[k]], H[k])
        c <- c(c, SWD[[k]][put_out == FALSE])
        h_val <- c(h_val, rep(H[k], length(c) - length(h_val)))
      }
    }
  }
  CP <- cbind(c, h_val)
  if (dim(CP)[1] > 1) {
    CP <- CP[order(CP[, 1]), ]
  }
  colnames(CP) <- c("peak", "h_value")
  if (dim(CP)[1] > 0) 
  {
    rownames(CP) <- as.character(1:dim(CP)[1])
  }
  
  R_ht   <- G_list
  
  
  
  
  SWD    <- lapply(SWD,FUN=function(vec){vec + (S-1)})
  CP[,1] <- CP[,1] + (S-1)
  
  # 4. Output 
  
  if (print.output) {
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
    if (method == "bootstrap") {
      cat("(Threshold Q derived by bootstrapping) \n\n")
    }
    if (method == "asymptotic") {
      cat("(Threshold Q derived by simulation) \n\n")
    }
    if (method == "fixed") {
      cat("(Threshold Q set by user) \n\n")
    }
    if (perform.CPD) {
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
    if (!perform.CPD) {
      cat("Peak detection was not performed")
    }
  }
  if (perform.CPD == FALSE) {
    CP <- NA
    SWD <- NA
  }
  
  tech.var <- list(D_ht=R_ht,x=x)
  mft <- list(M=M,Q=Q,method=method,sim=sim,blocksize=blocksize,CP=CP,S=S,E=E,Tt=Tt,H=H,alpha=alpha,two.sided=two.sided,
              perform.CPD=perform.CPD,tech.var=tech.var,type="peaks")
  class(mft) <- "MFT"
  invisible(mft)
}#end-MFT.peaks  
