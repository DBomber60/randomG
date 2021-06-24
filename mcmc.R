##### SAMPLER #####
source('sample.N.R')
set.seed(1)
nIter = 2000

# array of parameter samples for each intermediate variable
# theta (1)/ gamma (p)/ beta (p)/ sigsq (1)

# input: X[time, N, ncov]; A[N, t]; Y (1x1); p (covariates in xmod); pY (covariates in Y model)

gf.mcmc = function(A, X, Y, p, pY, nu_0, nu_1, nIter) {
  
  ncov = dim(X)[3] # number of intermediate covariates
  t = dim(A)[2] # number of time points
  
  # parameter arrays (assumes normal models)
  params = array(.1, dim = c(t-1, ncov, nIter, 2 * p + 2) ) # intermediate variables parameters
  paramsY = array(.1, dim = c(nIter, 2 * pY + 2))
  
  ##### initialize all parameters #####
  for (time in 1:(t-1)) { # for each time point (where we model the intermediate covariates)
    for (cv in 1:ncov) { # for each covariate
      m.init = lm(X[time + 1, , cv] ~ A[,time] + X[time, , cv]) 
      params[time, cv, 1, 1] = 0.5 # initial theta
      params[time, cv, 1, 2:(p+1)] = rep(1,p) # initial gamma
      params[time, cv, 1, (p+2):(2*p+1)] = coefficients(m.init) # initial beta
      params[time, cv, 1, 2*p+2] = sigma(m.init)^2 # initial sigsq
    }
  }
  
  paramsY[1, 1] = 0.5
  paramsY[1, 2:(pY+1)] = rep(1, pY)
  ymod = lm(Y ~ X[3,,] + A[,3])
  paramsY[1, (pY+2):(2*pY+1)] = coefficients(ymod)
  paramsY[1, 2*pY+2] = sigma(ymod)^2
  
  ##### MCMC #####
  
  for(it in 2:nIter) {
    
    for (time in 1:(t-1)) { # for each time point (where we model the intermediate covariates)
      
      for (cv in 1:ncov) { # for each covariate
        new.params = sample.new(theta.old = params[time, cv, it-1, 1],
                                gamma.old = params[time, cv, it-1, 2:(p+1)],
                                beta.old = params[time, cv, it-1, (p+2):(2*p+1)],
                                sigsq.old = params[time, cv, it-1, 2*p+2],
                                design = cbind(1, A[,time], X[time, , cv]),
                                resp = X[time + 1, , cv])
        
        params[time, cv, it, 1] = new.params$theta.new
        params[time, cv, it, 2:(p+1)] = new.params$gamma.new
        params[time, cv, it, (p+2):(2*p+1)] = new.params$beta.new
        params[time, cv, it, 2*p+2] = new.params$sigsq.new
      }
      
    }
    
    # then sample new value for the Y's
    new.params = sample.new(theta.old = paramsY[it-1, 1],
                            gamma.old = paramsY[it-1, 2:(pY+1)],
                            beta.old = paramsY[it-1, (pY+2):(2*pY+1)],
                            sigsq.old = paramsY[it-1, 2*pY+2],
                            design = cbind(1, X[3,,], A[,3]),
                            resp = Y)
    
    paramsY[it, 1] = new.params$theta.new
    paramsY[it, 2:(pY+1)] = new.params$gamma.new
    paramsY[it, (pY+2):(2*pY+1)] = new.params$beta.new
    paramsY[it, 2*pY+2] = new.params$sigsq.new
  }
  return(list(params = params, paramsY = paramsY))
}
