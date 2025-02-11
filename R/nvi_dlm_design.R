#' @title nvi_dlm_design: Construct System matrix (G) and Design matrix (F) 
#' 
#' @description The function \code{nvi_dlm_design} constructs the 
#' \code{G,F, W and V} design matrices for DLM analyses. It is (by now) 
#' possible to include a multivariate response with common trend and 
#' common seasonality modeled as harmonic waves. Furthermore, for the system 
#' equation nested factors and matèrn covariance might be included. 
#' @details ...
#' @param systemform [\code{formula}]. A object of class "formula", see
#' \link[stats]{formula}. The \code{formula} must be of the 
#' following form:\cr
#' \code{y ~ 1 + harmonic(wave_name) + nested(a \%in\% b) + matern(coordinates)}, 
#' where the different elements are:\cr
#' \item \code{y} (response). Name of the response variable. Might be any valid
#' variable name. Has to be represented in \code{data} as a \code{n x q} matrix.\cr
#' \item \code{1} or \code{0}. If \code{1} one common trend component is added 
#' to the DLM model.\cr 
#' \item harmonic(wave_name). Might be neglected. If part of the formula, 
#' harmonic waves with wave length defined in a vector (any length >=1) given 
#' in \code{data}.The name of the variable in data is the "\code{wave_name}" in the 
#' \code{harmonic(wave_name)} term. One common variance component defined in
#' the input parameter \code{varcomps} with name given by the argument 
#' ("\code{wave_name}") in the formula term \code{harmonic(wave_name)} constitutes 
#' the elements on the main diagonal of associated block of system variance 
#' \code{W}. If the argument is given as "\code{1|wave_name}", somewhat similar to 
#' the \code{lme4} package, there will be unique wave parameters for each response
#' variable. 
#' \item nested(a\%in\%b). Might be neglected. If part of the formula, two nested
#' factors, where factor "\code{a}" is nested in factor "\code{b}" might be added
#' to the model. The factors with names \code{a} and \code{b}, any valid names 
#' might be chosen as long as they are part of the argument in formula term 
#' \code{nested(a\%in\%b)} must be present in  \code{data}, both with length 
#' \code{p}, i.e. the number of observations per time step. Two variance 
#' components defined in the input parameter \code{varcomps} with names given 
#' by the factors of the argument in formula ("\code{nested(a\%in\%b)}") constitutes 
#' the elements on the main diagonals of associated blocks of system variance 
#' \code{W}. 
#' \item matern(coordinates). Might be neglected. If part of the formula, a matèrn 
#' covariance structure is added to the system covariance, \code{W}. A variable 
#' with name defined by the argument in formula term \code{matern(coordinates)} must
#' be present in data. As by now this variable has to be a matrix of size 
#' \code{p x 2} with longitude and latitude as it is passed directly to 
#' \link[sp]{:pDists} with argument "\code{longlat=TRUE}". The distance matrix
#' (\code{p x p}) resulting from this command is passed to command
#' \link[rSPDE]{matern.covariance}, with parameters 
#' \code{kappa = varcomps[[coordinates]]$kappa}, \code{nu = varcomps[[coordinates]]$nu}, 
#' and \code{sigma = varcomps$main}), where \code{coordinates} 
#' is the argument from formula term \code{matern(coordinates)}.
#' @param errorform [\code{formula}]. A object of class "formula", see 
#' \link[stats]{formula}. The \code{formula} must be of the 
#' following form:\cr
#' \code{.~ 1 +b + a +`...` + matern(coordinates_err)}, 
#' where the different elements are:\cr
#' \item \code{1} or \code{0}. If \code{1} one common trend component is added 
#' to the DLM model. Might be combined with factors. \cr
#' \item \code{a}, \code{b}.... Might be neglected. Factors for which blocks within the error 
#' covariance (\code{V}) shall have equal covariance elements identified by
#' \code{varcomps$error[[a]]}, \code{varcomps$error[[b]]} etc. 
#' \item matern(coordinates_err). Might be neglected. If part of the formula, a matèrn 
#' covariance structure is added to the error covariance, \code{V}. A variable 
#' with name defined by the argument in formula term \code{matern(coordinates_err)} must
#' be present in data. As by now this variable has to be a matrix of size 
#' \code{p x 2} with longitude and latitude as it is passed directly to 
#' \link[sp]{:pDists} with argument "\code{longlat=TRUE}". The distance matrix
#' (\code{p x p}) resulting from this command is passed to command
#' \link[rSPDE]{matern.covariance}, with additional parameters 
#' \code{kappa = varcomps[[coordinates_err]]$kappa}, 
#' \code{nu = varcomps[[coordinates_err]]$nu}, 
#' and \code{sigma = varcomps$error$main}). Separate or common coordinate matrices
#' and parameters for \code{nu} and \code{kappa} might be applied to the error 
#' and system equations/ formulas. 
#' @param varcomps[\code{list}] with elements: \cr
#' \item main. Must be included with exact name "main". Common system variance.
#' \item trend. Must be included with exact name "trend" if trend is included 
#' in system variance, i.e. \code{systemform = formula(. ~ 1 + ...+)}
#' \item error [\code{list}] with elements \code{main} (compulsory if common 
#' variance component is added in error variance, i.e. 
#' \code{errorform = formula(. ~ 1 + ...+)}). The main error might be a vector 
#' of length \code{p} which gives the main diagonal of the error variance matrix
#' (\code{V}). Furthermore,  elements \code{a_err}, \code{b_err},...
#' giving variances associated with blocks defined by factors of errorform 
#' elements named \code{a_err}, \code{b_err}, ... 
#' \item waves, one scalar giving variance component in system variance 
#' (\code{W}) associated with harmonic waves. 
#' \item a, b, scalars defining variance component for the blocks in the 
#' system variance (\code{W}) associated with nested factors is system equation. 
#' Names must be part of the argument in element \code{nested(a \%in\% b)} in 
#' \code{systemform}.
#' \item a_err, b_err, ..., scalars defining variance component for the blocks in the 
#' error variance (\code{V}) associated with factors defined as elements in 
#' elements of the \code{errorform}.
#' \item coordinates = [\code{list}] with elements \code{kappa} and 
#' \code{nu} used for the matèrn covariance in system variance (\code{W}). 
#' \item coordinates_err = [\code{list}] with elements \code{kappa_err} and 
#' \code{nu_err} used for the matèrn covariance in error variance (\code{V}).
#' @param data [\code{list}] with elements:\cr
#' \item q, [\code{integer}]. Giving the dimension of data if not defined in 
#' \code{systemform} (respnce) with associated data matrix \code{y}
#' \item y [\code{numeric matrix}] of size \code{n x q}, with data.
#' \item a, b, a_err,b_err,... factors of length \code{q} used in \code{systemform}
#'  and/ or \code{errorform}.
#'  \item waves [\code{numeric vector}] with wave lengths for harmonic waves.
#'  \item coordinates and/ or (might be the same) coordinates_err 
#'  [\code{numeric matrices}] giving longitude and latitude if maternal elements 
#'  are to be incorporated in system- (\code{W}) or error variance (\code{V}) 
#'  matrices.
#' @param isnullw [\code{logical}]. If the system variance is to be returned as 
#' \code{NULL}. Might be used if discount factor method is to be applied. 
#' @return A list with the following elements:\cr 
#' \item GG: System design (square) matrix G with dimension [q + trend == TRUE + 2 x number of waves]
#' \item FF: Observation design matrix with dimension [q + trend == TRUE + 2 x number of waves] x q 
#' \item W: System covariance matrix with dimension [q + trend == TRUE + 2 x number of waves]
#' \item V: Error covariance matrix with dimension q x q.  
#' @author Lars Erik Gangsei\cr
#' lars.erik.gangsei@@vetinst.no\cr
#' +47 950 61 231
#' 
#' @import Matrix
#' 
#' @examples
#' \dontrun{
#'
#' }
#'
#'
#' @export
nvi_dlm_design <- function(systemform = formula(.~ 0 + harmonic(waves) + 
                                                nested(a%in%b) + 
                                                matern(coordinates)),
                           errorform = formula(.~ 1 +b + a +`...` +
                                                 matern(coordinates_err)), 
                                   varcomps = list(main = 1,
                                                 trend = 10^-3,
                                                 error = list(main = 1,
                                                              a_err = 0.1,
                                                              b_err = 0.1),
                                                 waves = 10^-3,
                                                 a = NULL,
                                                 b = NULL,
                                                 coordinates = list(kappa = 10,
                                                               nu = 1/5),
                                                 coordinates_err = list(kappa = 10,
                                                                    nu = 1/5)),
                                   data = list(q = 1,
                                               y = NULL,
                                               a = NULL,
                                               b = NULL,
                                               a_err = NULL,
                                               b_err = NULL,
                                               waves = c(2*pi/365,2*2*pi/365),
                                               coordinates = NULL,
                                               coordinates_err = NULL),
                                  isnullw = FALSE)
  {
  # Find the terms in the formula, this is clumsy coding
  systemterms <- attributes(terms(systemform))$term.labels
  
  wave_name <- systemterms[grep('harmonic',systemterms)]
  if(length(wave_name)==0){wave_name <- NULL}else
  {wave_name <- gsub('harmonic|\\(|\\)| ','',wave_name)}
  
  if(grepl('1\\|',wave_name))
  {uniq_waves <- TRUE
  wave_name <- gsub('1\\|','',wave_name)
  }else{
    uniq_waves <- FALSE}
  
  
  nested_names <- systemterms[grep('nested',systemterms)]
  if(length(nested_names)==0){nested_names <- NULL}else
    {nested_names <- gsub(' ','',
                          unlist(strsplit(gsub('nested|\\(|\\)|','',
                                     nested_names),'%in%')))}
  
  matern_name <- systemterms[grep('matern',systemterms)]
  if(length(matern_name)==0){matern_name <- NULL}else
  {matern_name <- gsub('matern|\\(|\\)','',matern_name)}
  
  # Find the name of response if it sould be used to find qq
  response_name <- as.character(systemform)[2]
  
  trend <- attributes(terms(systemform))$intercept ==1
  

  # Find the number of locations (qq)
  if(!is.null(nested_names))
  {qq <- length(data[[nested_names[1]]])
  }else{
    if(!is.null(matern_name)){qq <-dim(data[[matern_name]])[1]}
      else{
        if(response_name!="."){
          qq <- dim(as.matrix(data[[response_name]]))[2]
        }else{qq <- max(data$q,1)}}}
  
  
  
  # Check if nested is part of the formula
  if(is.null(nested_names))
  {
    GG <- FF <- diag(qq)
    WW <- varcomps$main * diag(qq)
    colnames(GG) <- rownames(GG) <- colnames(FF) <- rownames(FF) <- 
      colnames(WW) <- rownames(WW) <- paste('main',1:qq,sep = '')
  }else{
      GG <- as.matrix(1)
      FF <- t(as.matrix(rep(1,qq)))
      WW <- as.matrix(varcomps$main)
      colnames(GG) <- rownames(GG) <- rownames(FF) <- 
        colnames(WW) <- rownames(WW) <- 'main_effect'
  }
 
  
  # Set the second block, i.e. trend if TRUE
  if(trend==TRUE)
  {
    GG <- cbind(rbind(GG,0),1)
    FF <- rbind(FF,0)
    WW <- as.matrix(Matrix::bdiag(WW,varcomps$trend))
    colnames(GG) <- rownames(GG) <- rownames(FF) <- 
      colnames(WW) <- rownames(WW) <- c(colnames(GG)[-dim(GG)[1]],'trend')
    
  } 
  
    # Set the third block, i.e. harmonics
    if(!is.null(wave_name))
    {
      t_names <- colnames(GG)
      # The number of waves 
      nn_h <- length(data[[wave_name]])
      for(hh in 1:nn_h)
      {
      wave_mat <- matrix(
        c(cos(data[[wave_name]][hh]),-sin(data[[wave_name]][hh]),
          sin(data[[wave_name]][hh]),cos(data[[wave_name]][hh])),
        byrow=FALSE,2,2)
      if(uniq_waves==TRUE){wave_mat <- Matrix::kronecker(diag(qq),wave_mat)}
      GG <- as.matrix(Matrix::bdiag(GG,wave_mat))
      }
      
      WW <- Matrix::bdiag(WW,varcomps[[wave_name]]*diag(2*nn_h*ifelse(uniq_waves==TRUE,qq,1)))
      
      FF_t <- as.matrix(rep(c(1,0),nn_h))
      if(uniq_waves==TRUE){FF_t <- Matrix::kronecker(diag(qq),FF_t)}else{
        FF_t <- Matrix::kronecker(t(as.matrix(rep(1,qq))),FF_t)
      }
      FF <- rbind(FF,FF_t)
      h_names <- paste('harmonic_',rep(rep(1:nn_h,each = 2)),
              rep(1:2,nn_h),sep = '')
      if(uniq_waves==TRUE){h_names <- paste('main_',rep(1:qq,each = length(h_names)),
                                            rep(h_names,qq),sep = '')}
      
      colnames(GG) <- rownames(GG) <- rownames(FF) <- 
        colnames(WW) <- rownames(WW) <- c(t_names,h_names)
    }
  
  
      # Set the factors, this code is not to good, so far only for two nested
    # factors. LEG 17.01.2025
    if(!is.null(nested_names))
    {
      t_names <- colnames(GG)
      if(length(intersect(nested_names,names(data)))!=2)
      {stop('Factor element only valid for two nested factors')}
        nested_fac <- as.factor(data[[nested_names[1]]])
        main_fac <- as.factor(data[[nested_names[2]]])
        
        nn_f <- c(nlevels(main_fac),nlevels(nested_fac))
        GG <- Matrix::bdiag(GG,diag(nn_f[2]-1))
        WW <- Matrix::bdiag(WW,varcomps[[nested_names[2]]]*diag(nn_f[1]-1),
                            varcomps[[nested_names[1]]]*diag(nn_f[2]-nn_f[1]))
        
        FF <- rbind(FF,t(model.matrix(lm(1:qq~main_fac,
                              contrast = list(main_fac = 'contr.sum'),
                              ))[,-1]))
        lev1 <- levels(main_fac)
        FFb <- cbind(diag(sum(main_fac==lev1[1])-1),-1)

        for(ii in 2:length(lev1))
        {
          FFb <- Matrix::bdiag(FFb,
                cbind(diag(sum(main_fac==lev1[ii])-1),-1))
        }
        FF <- rbind(FF,as.matrix(FFb))
        idx_dist <- which(rowSums(as.matrix(FFb)==1)==1)
        
        colnames(GG) <- rownames(GG) <- 
          colnames(WW) <- rownames(WW) <- c(t_names,
            paste(nested_names[2],1:(nn_f[1]-1),sep = ''),
            as.character(nested_fac)[idx_dist])
        rownames(FF) <- c(t_names,
                          paste(nested_names[2],1:(nn_f[1]-1),sep = ''),
                          as.character(nested_fac)[idx_dist])
      }
  
  # Set the matern covariance 
  if(!is.null(matern_name))
    {
    if(!is.null(nested_fac))
    {
      idx_ww <- dim(WW)[1]+seq(-nn_f[2]+nn_f[1]+1,0,by = 1)
    }
    else{
      idx_ww <- 1:qq#c(1,(dim(WW)[1]-qq+2):dim(WW)[1])
      idx_dist <- 1:qq
    }
      dist <- sp::spDists(data[[matern_name]][idx_dist,],longlat=TRUE)
      tmp_names <- colnames(WW)[idx_ww]
      WW[idx_ww,idx_ww] <- rSPDE::matern.covariance(dist,
                                              kappa = varcomps[[matern_name]]$kappa, 
                                              nu = varcomps[[matern_name]]$nu,
                                              sigma = varcomps$main) 
      
      colnames(WW)[idx_ww] <- rownames(WW)[idx_ww] <- tmp_names
    }
  
  
  ## Set the error variance. 
  # Find the terms in the formula, this is clumsy coding
  errorterms <- attributes(terms(errorform))$term.labels
  
  mainerror <- attributes(terms(errorform))$intercept ==1
  
  factor_terms <- errorterms[-grep('matern',errorterms)]
  
  matern_name <- errorterms[grep('matern',errorterms)]
  
  if(length(matern_name)==0){matern_name <- NULL}else
  {matern_name <- gsub('matern|\\(|\\)','',matern_name)}
  
  
  if((mainerror==TRUE)&&(is.null(matern_name)))
    {if(length(varcomps$error$main)==1){VV <-  varcomps$error$main * diag(qq)
    }else{VV <- diag(varcomps$error$main)}}else
  {VV <- matrix(0,qq,qq)}
  
  if(length(factor_terms)>0)
  {
    for(ff in factor_terms)
    {
    T_mat <- model.matrix(lm(1:qq~0+as.factor(data[[ff]]),data = data[ff]))
    VV <- VV + varcomps$error[[ff]] * (T_mat%*%t(T_mat)>0)
    }
  }
  
  if(!is.null(matern_name))
  {
    
  dist <- sp::spDists(data[[matern_name]],longlat=TRUE)
  VV <- rSPDE::matern.covariance(dist,
                kappa = varcomps[[matern_name]]$kappa, 
                nu = varcomps[[matern_name]]$nu,
                sigma = varcomps$error$main) 
  VV[is.na(VV)] <- 0
  }
  
  if(isnullw==TRUE){WW <- NULL}else{WW <- as.matrix(WW)}
  return(list(GG=as.matrix(GG),FF = as.matrix(FF),WW = WW,VV = VV))
}
