#' @title nvi_dglm: DGLM for multivariate binominal data. 
#' @description The function \code{DGLM_Multivariate_Binominal} fits a Dynamic Generalized
#' Linear Model (DGLM) to the data input. See "Bono, C., Cornou, C., Lundbye-Christensen, S., 
#' & Kristensen, A. R. (2013). Dynamic production monitoring in pig herds II. 
#' Modeling and monitoring farrowing rate at herd level. Livestock Science, 155(1), 92-102.
#' Doi: 10.1016/j.livsci.2013.03.026" for details.
#' @details This code is not yet implemented with the "formula option", ref. 
#' \link[NVI_dynamic_linear_models]{nvi_dlm}. The code might be merged to one
#' general function covering both DLM and DGLM as much of the framework is the same.
#' @param Data [\code{data.frame}]. Data frame with three elements:\cr
#' \code{t}: Vector with times/ dates.\cr
#' \code{n}: Matrix, where each row represents the number at risk for different
#' levels at associated time step.\cr
#' \code{y}: Matrix, where each row represents the number of incidents for different
#' levels at associated time step.
#' @param GG [\code{numeric matrix}] System matrix G. Assumed constant for all time elements.
#' @param FF [\code{numeric matrix}] Design matrix F. Assumed constant for all time elements.
#' @param m0 [\code{numeric vector}] Prior mean for (latent) parameter vector \eqn{\theta}.
#' @param C0 [\code{numeric matrix}] Prior variance for (latent) parameter vector \eqn{\theta}.
#' @param W [\code{numeric matrix}] System variance.
#' @param delta [\code{numeric scalar}]. Discount factor. Only used if \code{W} is \code{NULL}. 
#' @return A list with elements at, Rt,mt,Ct,ft,Qt,pt,etaHat,vHat,Qst
#'
#' @author Lars Erik Gangsei\cr
#' lars.erik.gangsei@@vetinst.no\cr
#' +47 950 61 231
#' 
#' @import dplyr
#' 
#' @examples
#' \dontrun{
#'
#' }
#'
#'
#' @export
nvi_dglm <- function(Data,GG,FF,m0, C0, W = NULL,delta = 0.75)
  {
  
  Data <- dplyr::arrange(Data,t)
  # m: length of theta (parametervector), n: The number of time steps
  mm <- dim(m0)[1]
  nn <- dim(Data)[1]
  qq <- dim(FF)[2]
  
  # Lists of matrices - initially they are empty, but every time step the new
  # matrices are added
  at.list <- list()
  Rt.list <- list()
  mt.list <- list()
  Ct.list <- list()
  ft.list <- list()
  Qt.list <- list()
  pt.list <- list()
  etaHat.list <- list()
  vHat.list <- list()
  Qst.list <- list()
  fst.list <- list()
  ft.all <- NULL
  pt.all <- NULL
  
  oldM <- m0
  oldC <- C0
  
  
  ## Run DGLM
  for (tt in 1:nn) { #n=1
    at <- GG %*% oldM              
    Rt <- GG %*%oldC%*%t(GG) + W          #same, no Gt
    
    # Save at and Rt
    at.list[[tt]] <- round(at,6)
    Rt.list[[tt]] <- round(Rt,6)
    
    
    yt <- as.matrix(Data$y[tt,])  #observation vector
    Nt <- as.matrix(Data$n[tt,]) #animals at risk
    
    # Remove observations with no observations of individuals at risk
    idxOK <- which(Nt>0)
    
    yt <- yt[idxOK,]
    Nt <- Nt[idxOK,]
    Ft <- as.matrix(FF[,idxOK])
    qqt <- dim(Ft)[2]
    
    # "Kalman filter" - Using Gt as an identity matrix
    ft <- t(Ft) %*% at                        # One-step Forecast probability
    #rownames(ft) <- rownames(at) # Give names to rows to be easier
    
    Qt <- t(Ft) %*% Rt %*% Ft                 # One-step Forecast variance
    Qt <- (Qt + t(Qt))/2                      # Make sure Qt is symmetrical
    
    pt <- 1/(exp(-ft)+1)                      # probability of being newly detected, hyperthermic, mild CS and severe CS
    
    pt <- matrix(pmax(1.0e-8, pt))            #If lower than 1.0e-8, replace by 1.0e-8
    pt <- matrix(pmin(1 - 1.0e-8, pt))        #If higher than 1 - 1.0e-8, replace by 1 - 1.0e-8
    #rownames(pt) <- rownames(ft)
    
    etaHat <- ft + (yt - Nt*pt)/(Nt*pt*(1-pt))  # logit of the probability
    vHatvec <- as.vector(1/(Nt*pt*(1-pt)))
    if(length(vHatvec)==1){vHat <- vHatvec}else{vHat <- diag(vHatvec)}
    
    #if (length(vHatVector) == 1) {
    # vHat = vHatVector 
    #}else {
    #  vect = vHatVector[,1]
    #  vHat = diag(vect)
    #}
    
    Qst <- solve(solve(Qt) + solve(vHat))
    fst <- Qst %*% (solve(Qt) %*% ft + solve(vHat) %*% etaHat) #estimated value of eta 
    
    mt <- at + Rt %*% Ft %*% solve(Qt) %*% (fst - ft)
    Ct <- Rt - Rt %*% Ft %*% (diag(qqt) - solve(Qt) %*% Qst) %*% solve(Qt) %*% t(Ft) %*% Rt
    Ct <- (Ct + t(Ct))/2
    
    at <- mt
    Rt <- Ct
    
    # Combined ft and pt of all iterations
    ft.all <- rbind(ft.all, ft)
    pt.all <- rbind(pt.all, pt)
    
    
    oldM <- mt
    oldC <- Ct
  
    # Save the results in lists
    mt.list[[tt]] <- round(mt,6)
    Ct.list[[tt]] <- round(Ct,6)
    
    ft.list[[tt]] <- round(ft,6)
    Qt.list[[tt]] <- round(Qt,6)
    pt.list[[tt]] <- round(pt,6)
    etaHat.list[[tt]] <- round(etaHat,6)
    vHat.list[[tt]] <- round(vHat,6)
    Qst.list[[tt]] <- round(Qst,6)
    fst.list[[tt]] <- round(fst,6)
    
  }
  
  return(list(      at=at.list,
                    Rt=Rt.list,
                    mt=mt.list,
                    Ct=Ct.list,
                    ft=ft.list,
                    Qt=Qt.list,
                    pt=pt.list,
                    etaHat=etaHat.list,
                    vHat=vHat.list,
                    Qst=Qst.list))
}
