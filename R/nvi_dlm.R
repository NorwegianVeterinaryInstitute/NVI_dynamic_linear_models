#' @title nvi_dlm: DLM for multivariate normal data. 
#' @description The function \code{nvi_dlm} fits a Dynamic Linear Model (DLM) 
#' to the data input. Either the three arguments \code{systemform}, 
#' \code{errorform}, and \code{varcomps} are used (sent directly to 
#' \link[NVI_dynamic_linear_models]{nvi_dlm_design}) for constructing
#' design matrices or the design matrices (\code{G,F,W and V}) are given as 
#' direct arguments.
#' @details The code is under development. The overall goal is to get a more 
#' generic framework in the "formulas", in line with formulas used in 
#' \link[stats]{lm} functions and functions like \link[lme4]{lmer}.  
#' @param systemform [\code{formula}]. A object of class "formula", see 
#' \link[stats]{formula}. If not \code{NULL} passed directly to 
#' function \link[NVI_dynamic_linear_models]{nvi_dlm_design}, together 
#' with \code{errorform} and \code{varcomps}. See 
#' \link[NVI_dynamic_linear_models]{nvi_dlm_design} for details.
#' @param errorform [\code{formula}]. A object of class "formula", see 
#' \link[stats]{formula}. If not \code{NULL} passed directly to 
#' function \link[NVI_dynamic_linear_models]{nvi_dlm_design}, 
#' together with \code{systemform} and \code{varcomps}. See 
#' \link[NVI_dynamic_linear_models]{nvi_dlm_design} for details.
#' @param varcomps[\code{list}]. If not \code{NULL} passed directly to 
#' function \link[NVI_dynamic_linear_models]{nvi_dlm_design}, together 
#' with \code{systemform} and \code{errorform}. 
#' See \link[NVI_dynamic_linear_models]{nvi_dlm_design} 
#' for details. 
#' @param data [\code{data.frame}]. List with at least two elements:\cr
#' \code{t}: Vector with times/ dates.\cr
#' \code{y}: Matrix, of size \code{N x q} where each row represents the 
#' observations for \code{q} different sites/ levels at associated time step. 
#' (\code{N} timesteps in total). Might be given name in accordance with
#' \code{systemform}. 
#' If the "formula" option is used, \code{data} must also contain elements 
#' defining \code{nested factors}(of length \code{q}), \code{wavelengths} defining
#' harmonic waves, \code{distance matrix} (of size \code{q x 2}) associated with elements included
#' in \code{systemform} and \code{errorform} respectively.    
#' @param GG [\code{numeric matrix}] System matrix G. 
#' Assumed constant for all time elements. 
#' Not used if "formula option is applied"
#' @param FF [\code{numeric matrix}] Design matrix F. 
#' Assumed constant for all time elements. 
#' Not used if "formula option is applied"
#' @param m0 [\code{numeric vector}] Prior mean for (latent) 
#' parameter vector \eqn{\theta}.
#' @param C0 [\code{numeric matrix}] Prior variance for (latent) 
#' parameter vector \eqn{\theta}.
#' @param W [\code{numeric matrix}] System variance.
#' Assumed constant for all time elements. 
#' Not used if "formula option is applied"
#' @param V [\code{numeric matrix}] Error variance.
#' Assumed constant for all time elements. 
#' Not used if "formula option is applied"
#' @param delta [\code{numeric scalar}]. Discount factor. 
#' If not \code{NULL} discount factor is applied and input argument 
#' \code{W} (system variance) is ignored, also when "formula option" is
#' applied. \code{delta} might be given as a matrix of size \code{p x p}
#' @param smoother_run [\code{logical}]. If the DLM smoother is to 
#' be run
#' @param isnullw [\code{logical}]. If the system variance is to be returned as 
#' \code{NULL}. Must be set till \code{TRUE} if discount factor method is to 
#' be applied in combination with "formula definition". 
#' @return A list with elements:\cr
#' \item at = prior means (\code{N x p}) for parameter vector \eqn{\theta},\cr
#' \item Rt = array of size (\code{N x p x p}), with prior variances for parameter vector \eqn{\theta},\cr
#' \item mt = prior means (\code{N x p}) for parameter vector \eqn{\theta},\cr
#' \item Ct = array of size (\code{N x p x p}), with posterior variances for parameter vector \eqn{\theta},\cr
#' \item ft = forecasts for observations matrix (\code{N x q}),\cr
#' \item QT = array of size (\code{N x q x q}), with forecast variances.\cr
#' \item At = array of size (\code{N x p x q}), with AT matrices...,\cr
#' \item et = forecast errors (\code{N x p}),
#' \tem log_like = marginal log liikelihoods for each timestep (length \code{N}).\cr
#' \item Qtsmooth = NULL, "smoothed forcast variances", only used if \code{smoother_run = TRUE}.\cr
#' \item atsmooth = NULL, "smoothed forcasts", only used if \code{smoother_run = TRUE}.\cr
#' \item Rtsmooth = NULL, "smoothed R matrices", only used if \code{smoother_run = TRUE}.
#'
#' @author Lars Erik Gangsei\cr
#' lars.erik.gangsei@@vetinst.no\cr
#' +47 950 61 231
#' 
#' @import dplyr
#' @import mvtnorm
#' 
#' @examples
#' \dontrun{
#'
#' }
#'
#'
#' @export
nvi_dlm <- function(systemform = NULL,errorform = NULL,
                    varcomps = list(main = 1,trend = 10^-3,
                                    error = list(main = 1,a = 0.1,b = 0.1),
                                    waves = 10^-3,a = NULL,b = NULL,
                                    dist = list(kappa = 10,nu = 1/5,sigma_d2=1)), 
                    data,GG=NULL,FF=NULL,m0=NULL, 
                    C0=NULL, W = NULL,V=NULL, 
                    delta = 0.75,smoother_run = FALSE,isnullw = FALSE)
  {
  if(!is.null(systemform))
  {
    design_matrices <- nvi_dlm_design(systemform = systemform,
                                      errorform = errorform, 
                                      varcomps = varcomps,
                                      data = data,isnullw = isnullw) 
    W <- design_matrices$WW
    V <- design_matrices$VV
    FF <- design_matrices$FF
    GG <- design_matrices$GG
    
    resp_name <- as.character(systemform)[2]
    
    if(is.element('t',names(data)))
    {data <- data.frame(y = data[[resp_name]],t = data$t)}else{
      data <- data.frame(y = data[[resp_name]],t = 1:dim(data[[resp_name]])[1])   
      }
  }
  
  
  data <- dplyr::arrange(data,t)
  # pp: length of theta (parametervector), n: The number of time steps
  pp <- dim(GG)[1]
  if(is.null(m0))
    {
    m0 <- as.matrix(rep(0,pp))
    }else
    {
      if(dim(m0)[1]!=pp)
      {
      m0 <- matrix(m0,pp,1)
      warning('Inconsitent size of m0')
      }
    }
  if(is.null(C0))
  {
    C0 <- 10^5*diag(pp)
  }else
  {
    if(dim(C0)[1]!=pp)
    {
      stop('wrong size of C0 (initial prior variance)')
    }
  }
  
 
  nn <- dim(data)[1]
  qq <- dim(FF)[2]
  
  # If discount factor is scalar, make it a matrix.
  if(!is.null(delta)){
  if(!is.matrix(delta)){delta<- delta*diag(pp)}
    delta_mat <- diag(pp)+(diag(pp)-as.matrix(delta))%*%solve(as.matrix(delta))
  }
  
  
  # Lists of matrices - initially they are empty, but every time step the new
  # matrices are added
  Res <- list(
    at = matrix(NA,nn,pp),
    Rt = array(NA,dim = c(nn,pp,pp)),
    mt = matrix(NA,nn,pp),
    Ct = array(NA,dim = c(nn,pp,pp)),
    ft = matrix(NA,nn,qq),
    QT = array(NA,dim = c(nn,qq,qq)),
    At = array(NA,dim = c(nn,pp,qq)),
    et = matrix(NA,nn,qq),
    Qtsmooth = NULL,
    atsmooth = NULL,
    Rtsmooth = NULL,
    log_like <- rep(NA,nn))
  
  oldM <- m0
  oldC <- C0
  
  
  ## Run DLM
  for (tt in 1:nn) { #n=1
    starttime <- Sys.time()
    
    at <- GG %*% oldM   
    if(!is.null(W))
      {
      Rt <- as.matrix(GG %*%oldC%*%t(GG) + W)
    }else{
      Rt <- as.matrix(GG %*%oldC%*%t(GG)%*%solve(delta_mat))
      }
    
    # Save at and Rt
    Res$at[tt,] <- round(as.vector(at),6)
    Res$Rt[tt,,] <- round(Rt,6)
    
    yt <- as.matrix(data$y[tt,])  #observation vector
    
    # Remove observations with no observations of individuals at risk
    idxOK <- which(!is.na(yt))
    
    yt <- as.matrix(yt[idxOK,])
    Ft <- as.matrix(FF[,idxOK])
    qqt <- dim(Ft)[2]
    
    # One step forecast 
    #print(tt)
    ft <- t(Ft)%*%at                                        # Forecast
    Qt <- t(Ft)%*%Rt%*%Ft + V[idxOK,idxOK]                  # Forecast variance
    
    Res$ft[tt,idxOK] <- as.vector(ft)
    Res$QT[tt,idxOK,idxOK] <- Qt
    
    # Posterior for theta at time t
    
    Qt_inv <- eigen(Qt,symmetric = TRUE)
    Qt_inv$vectors <- Qt_inv$vectors[,Qt_inv$values>10^(-99)]
    Qt_values <- Qt_inv$values[Qt_inv$values>10^(-99)]
    
    Qt_inv <- Qt_inv$vectors%*%diag(1/Qt_values)%*%t(Qt_inv$vectors)
    endtime1 <- Sys.time()
    
    At <- Rt%*%Ft%*%Qt_inv                          # Adaptive coefficient
    et <- matrix(NA,qq,1)
    et[idxOK,] <- as.vector(yt - ft)                # Forecast error
    mt <- at + At%*%et[idxOK,]                      # Filtered (posterior) mean
    
    # Get the marginal log-likelihood.
    Res$log_like[tt] <- -(1/2)*(dim(Qt)[1]*log(2*pi)+
                                  sum(log(Qt_values))+
                                  t(et[idxOK,])%*%Qt_inv%*%et[idxOK,])
    
    
    endtime2 <- Sys.time()
    
    Ct <- Rt - At%*%Qt%*%t(At)                      # Filtered (posterior) variance
   
    Res$mt[tt,] <- as.vector(mt)
    Res$Ct[tt,,] <- Ct
    Res$At[tt,,idxOK] <- At
    Res$et[tt,] <- as.vector(et)
    
    oldM <- mt
    oldC <- Ct
    
    endtime3 <- Sys.time()
    #print(starttime-c(endtime1,endtime2,endtime3))
  
  }
  if(smoother_run == TRUE)
  {
   # Res$log_like <- rep(NA,nn)
    Res$Qtsmooth <- array(NA,dim = c(nn,qq,qq))
    Res$Ctsmooth <- array(NA,dim = c(nn,pp,pp))
    Res$Rtsmooth <- array(NA,dim = c(nn,pp,pp))
    Res$atsmooth <- matrix(NA,nn,pp)
    
    for(tt in seq(nn,1,by = -1))
    {
    
    #print(tt)
      if(tt == nn)
      {
        at <- mt                # Smoothed mean
        Ct <- Ct                # Smoothed system variance
        Rt <- Rt
      }else{
      Rtp1_inv <- eigen(Res$Rt[tt+1,,],symmetric = TRUE)
      Rtp1_inv$vectors <- Rtp1_inv$vectors[,Rtp1_inv$values>10^(-99)]
      Rtp1_inv$values <- Rtp1_inv$values[Rtp1_inv$values>10^(-99)]
      Rtp1_inv <- Rtp1_inv$vectors%*%diag(1/Rtp1_inv$values)%*%t(Rtp1_inv$vectors)
      
      Bt <- Res$Ct[tt,,]%*%t(GG)%*%Rtp1_inv
      Ct <- as.matrix(Res$Ct[tt,,] + Bt%*%(Res$Ctsmooth[tt+1,,] - Res$Rt[tt+1,,])%*%t(Bt)) # Smoothed system variance 
      at <- Res$mt[tt,] + Bt%*%(Res$atsmooth[tt+1,] - Res$at[tt+1,])       # Smoothed mean
      }
      Res$Ctsmooth[tt,,] <- Ct
      Res$Rtsmooth[tt,,] <- Rt
      Res$atsmooth[tt,] <- as.vector(at)
    }
  }
  
  return(Res)
}
