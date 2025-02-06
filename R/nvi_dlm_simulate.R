#' @title nvi_dlm_simulate: DLM for multivariate normal data. 
#' @description The function \code{DGLM_simulate} simulate data from a Dynamic 
#' Linear Model (DLM). 
#' @details ...
#' @param mu0 [\code{numeric vector}] Expected value for \eqn{\theta[0]}
#' @param GG [\code{numeric matrix}] System matrix G \code{size (q x q)}. 
#' Assumed constant for all time elements.
#' @param FF [\code{numeric matrix}] Design matrix F, \code{size (q x p)}. 
#' Assumed constant for all time elements.
#' @param W [\code{numeric matrix}] \code{size (q x q)}. System variance.
#' @param V [\code{numeric matrix}] \code{size (p x p)}. Error variance.
#' @param n [\code{integer}]. 
#' @return A data frame with elements:\cr
#' \item t [\code{numeric vector}] representing time steps (\code{1:N}).
#' \item y [\code{numeric matrix}] \code{size (N x p)} representing 
#' simulated observations. 
#' \item theta [\code{numeric matrix}] \code{size (N x q)} representing 
#' simulated observations expected values, i.e. hidden data. 
#' 
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
nvi_dlm_simulate <- function(mu0,GG,FF,W,V,N = 100)
  {
  qq <- dim(GG)[1]
  pp <- dim(FF)[2]
  Res <- data.frame(t = 1:N,
                    y = I(matrix(NA,N,pp)),
                    theta = I(matrix(NA,N,qq)))
  
  # Simulate errors for theta (system errors)
  ww <- mvtnorm::rmvnorm(N,mean = rep(0,qq),sigma = W)
  
  # Simulate errors for y (observation errors)
  vv <- mvtnorm::rmvnorm(N,mean = rep(0,pp),sigma = V)
  
  for(tt in 1:N)
  {
    if(tt == 1){Res$theta[1,] <- mu0 + ww[1,]}else{
      Res$theta[tt,] <- GG %*% Res$theta[tt-1,] + ww[tt,]
    }
    Res$y[tt,] <- t(FF) %*% Res$theta[tt,] + vv[tt,]
  }
  
 
  return(Res)
}
