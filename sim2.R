# simulation analysis for > 2 treatment times
# let's start out with 2 treatment times and get things working

library(mvtnorm)
source('sample.N.R')
set.seed(1)

##### SIMULATE DATA #####
N = 1000 # number of subjects
ncov = 3 # 3 covariates

t = 3 # 3 time points
u = rnorm(N, sd = 0.5)
s <- matrix(0.25, ncov, ncov) # sigma X
diag(s) <- 1.0 #
xdiffsd = 0.25

su <- matrix(0.05, ncov, ncov) # sigma U
diag(su) <- 0.1
u = rmvnorm(N, mean = rep(0, ncov), sigma = su)

# hold time-varying covariates in multi-dimensional array
X = array(0, dim = c(3, N, ncov)) # 2 time points
A = matrix(0, nrow = N, ncol = t)

X[1,,] = rmvnorm(N, rep(0, ncov), su)

A[,1] = rbinom(N, 1, prob = plogis(1/3 * rowSums(X[1,,]))) # suppose the first treatment is (unconditionally) randomized

X[2,,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[1,,] + 
  0.6 * cbind(-A[,1], 0, 0) + u # treatment has effect on only the first column 

A[,2] = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X[2,,]) )  ) # suppose the first treatment is (unconditionally) randomized

X[3,,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[2,,] + 
  0.6 * cbind(-A[,2], 0, 0) + u # treatment has effect on only the first column 

A[,3] = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X[3,,]) )  ) # suppose the first treatment is (unconditionally) randomized

Y = 0.8 * X[3,,1] + rnorm(N, sd = 0.2) + u[,1] - 0.6 * A[,3]
#X_2 = rmvnorm(ntot, rep(0, ncov), (xdiffsd^2) * s) + 0.8 * X_1 + 0.6 * cbind(-A_0, 0, 0) + u

# for each covariate, model the inclusion


##### SAMPLER #####

set.seed(1)
nIter = 1000

nu_1 = 5
nu_0 = 0.02

p = 3 # covariate model design matrix dimension
pY = 5 # outcome model design matrix dimension

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


# posterior predictive sampler

dat$y11 = NA
dat$y00 = NA
dat$t1 = ifelse(dat$A1 == 1 & dat$A2 == 1, 1, 0)
dat$t0 = ifelse(dat$A1 == 0 & dat$A2 == 0, 1, 0)

# first, let's do this with the correct G
row = 100
ndraws = 500

# important indices. params array [time, covariate, row, column]
betas = 5:7 # beta indices for covariate models
sigsq = 8 # sigsq indices for covariate model
# design = cbind(1, A[,time], X[time, , cv]) # COVARIATE MODELS
# design = cbind(1, X[3,,], A[,3]) # OUTCOME MODEL

burn = 100
params.r = params[,,(burn+1):nIter,] # remove 100 rows for burn-in
paramsY.r = paramsY[(burn+1):nIter,]
X.pp = array(0, dim = c(t-1, N, ncov))
ndraws = 10

ypp = array(NA, dim = c(N, ndraws) )

for (draw in 1:ndraws) {
  for(time in 1:(t-1)) {
    for (covs in 1:ncov) {
      # mean function
      mu_hat = cbind(1, A[,time], X[time,,covs]) %*% params.r[time, covs, iter, betas]
      X.pp[time, ,covs] = rnorm(N, mean = mu_hat, sd = sqrt(params.r[time, covs, iter, sigsq]) )
    }
  }
  mu_hat.y = cbind(1, X.pp[2,,], A[,3]) %*% paramsY[iter + burn, 7:11]
  ypp[,draw] = rnorm(N, mean = mu_hat.y, sd = sqrt(paramsY[iter + burn, 12]))
}

cor(ypp[,2], Y)

diffs = array(0, dim = ndraws)

for(draw in 1:ndraws) {
  
  if ( params[1, row + draw, 3] == 1 ) { # if there is an edge from treatment to covariate, impute
    
    mx1 = cbind(dat$X1, 1) %*% params[1, row + draw, c(4,5)]  # XB (treated at time 1)
    X2draw1 = rnorm(n, mean = mx1, sd = params[1, row + draw, 6]) # fixed treatment to 1
    X2draw1 = ifelse(dat$A1 == 1, dat$X2, X2draw1)
    
    
    mx0 = cbind(dat$X1, 0) %*% params[1, row + draw, c(4,5)]  # XB (treated at time 1)
    X2draw0 = rnorm(n, mean = mx0, sd = params[1, row + draw, 6]) # fixed treatment to 1
    X2draw0 = ifelse(dat$A1 == 0, dat$X2, X2draw0)
    
    
  } else { 
    
    X2draw1 = dat$X2
    X2draw0 = dat$X2
    
  }
  
  if ( params[2, row + draw, 3] == 1 ) { # if there is an edge from treatment to covariate, impute
    
    mz1 = cbind(dat$Z1, 1) %*% params[2, row + draw, c(4,5)]  # XB (treated at time 1)
    Z2draw1 = rnorm(n, mean = mz1, sd = params[2, row + draw, 6]) # fixed treatment to 1
    Z2draw1 = ifelse(dat$A1 == 1, dat$Z2, Z2draw1)
    
    
    mz0 = cbind(dat$Z1, 0) %*% params[1, row + draw, c(4,5)]  # XB (treated at time 1)
    Z2draw0 = rnorm(n, mean = mz0, sd = params[1, row + draw, 6]) # fixed treatment to 1
    Z2draw0 = ifelse(dat$A1 == 0, dat$Z2, Z2draw0)
    
    
  } else { 
    
    Z2draw1 = dat$Z2
    Z2draw0 = dat$Z2
    
  }
  
  
  
  
  
  my1 = cbind(X2draw1, Z2draw1, 1, 1) %*% paramsY[(row+draw),c(6:9)] # always treated
  Ydraw1 = rnorm(n, mean = my1, sd = paramsY[(row+draw), 10])
  
  my0 = cbind(X2draw0, Z2draw0, 1, 1) %*% paramsY[(row+draw),c(6:9)] # always treated
  Ydraw0 = rnorm(n, mean = my0, sd = paramsY[(row+draw), 10])
  
  
  dat$y11 = ifelse(dat$t1 == 1, dat$Y, Ydraw1)
  dat$y00 = ifelse(dat$t0 == 1, dat$Y, Ydraw0)
  
  diffs[draw] = sum(dat$y11 - dat$y00)/n
  
  
}

#hist(diffs)
d.random = data.frame(diffs = diffs, g = "ran" )


