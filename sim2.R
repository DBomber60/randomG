# simulation analysis for > 2 treatment times
# let's start out with 2 treatment times and get things working

# compute the true effect

library(mvtnorm)
library(tidyverse)
source('sample.N.R')
set.seed(1)

##### SIMULATE DATA #####
N = 1000 # number of subjects
ncov = 3 # 3 covariates

t = 3 # 3 time points
#u = rnorm(N, sd = 0.75)
s <- matrix(0.25, ncov, ncov) # sigma X
diag(s) <- 1.0 #
xdiffsd = 1 # used to be 1/4

su <- matrix(0.07, ncov, ncov) # sigma U
diag(su) <- 0.1
u = rmvnorm(N, mean = rep(0, ncov), sigma = su)

# hold time-varying covariates in multi-dimensional array
X = array(0, dim = c(3, N, ncov)) # 2 time points
A = matrix(0, nrow = N, ncol = t)

X[1,,] = rmvnorm(N, rep(0, ncov), su)

A[,1] = rbinom(N, 1, prob = plogis(1/3 * rowSums(X[1,,]))) # suppose the first treatment is (unconditionally) randomized

X[2,,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[1,,] + 
  0.6 * cbind(-A[,1], 0, 0) + u # treatment has effect on only the first column 

A[,2] = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X[2,,]) )  ) # 

X[3,,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[2,,] + 
  0.3 * cbind(-A[,2], 0, 0) + u # treatment has effect on only the first column 

A[,3] = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X[3,,]) )  ) # suppose the first treatment is (unconditionally) randomized

Y = 0.8 * X[3,,1] + rnorm(N, sd = 0.2) + u[,1] - 0.6 * A[,3]
#X_2 = rmvnorm(ntot, rep(0, ncov), (xdiffsd^2) * s) + 0.8 * X_1 + 0.6 * cbind(-A_0, 0, 0) + u

# for each covariate, model the inclusion


##### SAMPLER #####

set.seed(1)
nIter = 5000

nu_1 = 5
nu_0 = 0.02

p = 3 # covariate model design matrix dimension
pY = p + 2 # outcome model design matrix dimension: all covariates plus treatment/ intercept

# array of parameter samples for each intermediate variable
# theta (1)/ gamma (p)/ beta (p)/ sigsq (1)
params = array(.1, dim = c(t-1, ncov, nIter, 2 * p + 2) ) # intermediate variables
paramsY = array(.1, dim = c(nIter, 2 * pY + 2))


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
  pY = 5
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


##### posterior predictive sampler #####


# important indices. params array [time, covariate, row, column]
betas = 5:7 # beta indices for covariate models
sigsq = 8 # sigsq indices for covariate model
# design = cbind(1, A[,time], X[time, , cv]) # COVARIATE MODELS
# design = cbind(1, X[3,,], A[,3]) # OUTCOME MODEL

burn = 100
params.r = params[,,(burn+1):nIter,] # remove 100 rows for burn-in
paramsY.r = paramsY[(burn+1):nIter,]


# TODO
# posterior predictive density for Y^{111}
# posterior predictive density for Y^{000}


# important index
edge.index = 4 # index for edge indicator in posterior sample
ndraws = 3000
ypp111 = array(NA, dim = c(N, ndraws) )

for (draw in 1:ndraws) {
  X.pp = array(0, dim = c(t-1, N, ncov))
  for(time in 1:(t-1)) {
    for (covs in 1:ncov) {
      # mean function
      mu_hat = cbind(1, 1, X[time,,covs]) %*% params.r[time, covs, draw, betas]
      X.pp[time, ,covs] = rnorm(N, mean = mu_hat, sd = sqrt(params.r[time, covs, draw, sigsq]) )
      
      if (params.r[time, covs, draw, edge.index] == 0) { X.pp[time, ,covs] = X[time+1, ,covs] }
    }
  }
  mu_hat.y = cbind(1, X.pp[2,,], 1) %*% paramsY.r[draw, 7:11]
  ypp111[,draw] = rnorm(N, mean = mu_hat.y, sd = sqrt(paramsY.r[draw, 12]))
}


ypp000 = array(NA, dim = c(N, ndraws) )

for (draw in 1:ndraws) {
  X.pp = array(0, dim = c(t-1, N, ncov))
  for(time in 1:(t-1)) {
    for (covs in 1:ncov) {
      # mean function
      mu_hat = cbind(1, 0, X[time,,covs]) %*% params.r[time, covs, draw, betas]
      X.pp[time, ,covs] = rnorm(N, mean = mu_hat, sd = sqrt(params.r[time, covs, draw, sigsq]) )
      
      if (params.r[time, covs, draw, edge.index] == 0) { X.pp[time, ,covs] = X[time+1, ,covs] }
    }
  }
  mu_hat.y = cbind(1, X.pp[2,,], 0) %*% paramsY.r[draw, 7:11]
  ypp000[,draw] = rnorm(N, mean = mu_hat.y, sd = sqrt(paramsY.r[draw, 12]))
}

colSums(params[2,1,,])

j = colSums(ypp111 - ypp000)/N
quantile(j)
hist(colSums(ypp111 - ypp000)/N)

# 0.25: -.57 0.75: -.54




