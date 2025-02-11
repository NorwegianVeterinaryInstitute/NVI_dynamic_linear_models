# NVI_dynamic_linear_models
Repository with R-package for running Dynamic Linear Models.

# Installation of the NVI.dynamic.linear.models package
The package might be installed by the following commands: 

```
require(devtools)
devtools::install_github('https://github.com/NorwegianVeterinaryInstitute/NVI_dynamic_linear_models')
library(NVI.dynamic.linear.models)
```

# A code excample
The code in the following paragraphs will simulate and evaluate 
```
rm(list = ls())

```
## Simulate data
### Formula set up for model
Set up a system model with no trend, and harmonic waves defined by the argument
`week_wave` in data. Furthermore, a matèrn co variance structure is added to the
error/ observation equation, where the argument `long_lat` gives geographical 
position of the different "sites" represented in the data.

```
system_form <- formula(. ~ 0 + harmonic(week_wave))
error_form <- formula(.~ 1 + matern(long_lat))
```

### Wave length, latitudes and longnitudes.
Set the wave length so that there is one wave per year with weekly data (52 weeks). Simulate latitudes and longitudes. 

```
set.seed(15)
qq <- 20 #Number of locations, second dimention of response.
NN <- 200 # Number of time points/ weeks, first dimention of response.
sim_inp <- list(week_wave = c(2*pi/52),
                long_lat = I(cbind(rnorm(qq,20,5),rnorm(qq,62,8))),
                q = qq)
```

### Variance components
Set the variance components defining the system variance (`W`) and observation variance (`V`). These variance components are not totally "free", but in the present case we let there be four free parameters defined in one variable `sigma_vec`. We assume common variance components for all system elements, i.e. main effects and harmonic waves, defined by the first element in `sigma_vec`, furthermore a main error variance (main diagonal of `V`). The last element is the proportion between $\nu$ and $\kappa$ parameters passed to the Matern Covariance function.  

```
sigma_vec <- c(10^-2,0.9,10)

sigma_list <- list(main = sigma_vec[1], trend = NULL,
                   error = list(main = sigma_vec[2]), 
                   week_wave = sigma_vec[1],
                   long_lat = list(kappa = 0.001*sigma_vec[3], 
                                   nu = sigma_vec[3]))
```

### Design matrices
Based on input above the formula `nvi_dlm_design()` make design matrices:

```
design_mat <- nvi_dlm_design(systemform = system_form,
                             errorform = error_form,
                             varcomps = sigma_list,
                             data = sim_inp)

GG <- design_mat$GG
FF <- design_mat$FF
WW <- design_mat$WW
VV <- design_mat$VV

image(WW)
image(VV)
image(GG)
image(FF)
```

### Simulation
Data simulation via function `nvi_dlm_simulate`:

```
sim_data <- nvi_dlm_simulate(mu0 = as.matrix(c(-5,1,1,rep(-5,qq-1))),
                             GG = GG,FF = FF,W =WW,V = VV,N = NN)

sim_inp$y <- I(sim_data$y)

matplot(sim_inp$y)
```

### Maximum likelihood estimation of variance components
The (approximate) marginal log-likelihood is one of the outputs from `nvi_dlm()`. This is utilized via the `nlm()`function to find maximum likelihood estimates for the three variance components.  

```
ml_func <- function(sigma_vec,systemform,errorform,data)
{
  # Make sure that all variance components are positive
  if(any(sigma_vec<=0)){res <- Inf}else{
  #sigma_vec <- sqrt(sigma_vec^2)
  
  # Make the input to variance components a function of first argument
  sigma_list <- list(main = sigma_vec[1], trend = NULL,
                   error = list(main = sigma_vec[2]), 
                   week_wave = sigma_vec[1],
                   long_lat = list(kappa = 0.001*sigma_vec[3], 
                                   nu = sigma_vec[3]))
  
  # Fit DLM and calculate log likelihood
  res <- nvi_dlm(systemform = systemform,errorform = errorform,
                 varcomps = sigma_list,data = data)
  res <- -sum(res$log_like)
  #res <- sum(res$et^2,na.rm=TRUE)
  }
  # Return the negative log_likelihood
  return(res)
}


# Use nlm function to find maximul likelihood estimates
ml_eval <- nlm(f = ml_func,p = c(1,1,25),
               systemform = formula(y ~ 0 + harmonic(week_wave)),
               errorform = error_form,data = sim_inp,
               hessian = TRUE,steptol = 10^-3)
```


### Approximate confidence intervals for variance components
By using the `hessian = TRUE` argument in the `nlm()`function we get the observed fisher information, which makes it possible to find approximate confidence intervals for the variance components.

```
Estimates <- data.frame(oracle = sigma_vec,
                        estimate = ml_eval$estimate,
                        LL = ml_eval$estimate +
                        qnorm(0.025)*sqrt(diag(solve(ml_eval$hessian))),
                        UL = ml_eval$estimate +
                        qnorm(0.975)*sqrt(diag(solve(ml_eval$hessian))))
                        
plot(Estimates$estimate,pch = 19,col = 'red',main = 'True vs. estimated',
xlab = 'Variance components 1- 3',ylim = c(min(Estimates$LL),max(Estimates$UL))) 

points(Estimates$oracle,pch = 20,col = 'blue')
points(Estimates$LL,pch = '-',cex = 3)
points(Estimates$UL,pch = '-',cex = 3)


```