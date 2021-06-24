library(mvtnorm)
library(tidyverse)
source('sample.N.R')
source('mcmc.R')
set.seed(1)

# simulation analysis for > 2 treatment times
# let's start out with 2 treatment times and get things working
# compute the true effect
# let's do 6 covariates and have treatment effect go through two 



##### SIMULATE DATA #####
N = 1000 # number of subjects
ncov = 6 # 3 covariates

t = 3 # 3 treatment time points
s <- matrix(0.25, ncov, ncov) # sigma X
diag(s) <- 1.0 #
xdiffsd = .25 # used to be 1/4

su <- matrix(0.05, ncov, ncov) # sigma U
diag(su) <- 0.1
u = rmvnorm(N, mean = rep(0, ncov), sigma = su)

# hold time-varying covariates in multi-dimensional array
X = array(0, dim = c(t, N, ncov)) # 2 time points
A = matrix(0, nrow = N, ncol = t) # treatment matrix; each column is a treatment time

X[1,,] = rmvnorm(N, rep(0, ncov), sigma = su)

A[,1] = rbinom(N, 1, prob = plogis(1/3 * rowSums(X[1,,]))) # 

X[2,,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[1,,] + 
  - 0.6 * cbind(A[,1], 0, 0, 0, 0, 0) + u     

A[,2] = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X[2,,]) )  ) # 

X[3,,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[2,,] + 
  - 0.6 * cbind(A[,2], 0, 0, 0, 0, 0) + u # treatment has effect on only the first column 

A[,3] = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X[3,,]) )  ) # suppose the first treatment is (unconditionally) randomized

Y = 0.8 * X[3,,1] + rnorm(N, sd = 0.2) + u[,1] - 0.6 * A[,3]

### sample posterior parameter values
nu_1 = 2
nu_0 = 0.02

p = 3 # covariate model design matrix dimension (intercept/ last measurement/ last treatment)
pY = ncov + 2 # outcome model design matrix dimension: all covariates plus treatment/ intercept
nIter = 1000


# sensitivity to nu+0
nu_0 = c(.001, .02, .05)
plotdf = array(0, dim = c(length(nu_0), 3))

for (m in 1:length(nu_0)) {

samples = gf.mcmc(A, X, Y, p, pY, nu_0[m], nu_1, nIter)

# look at graph characteristics
p = 2 # forget about intercept
edgesX = array(NA, dim = c(t-1, nIter, ncov * p ))
# gather the gammas and label them

for(time in 1:(t-1)) {
  for(c in 1:ncov) {
      edgesX[time, ,((c-1)*p+1): ((c-1)*p+p)] = samples$params[time, c, , 3:4]
  }
}

truev = rep( c(1,1, rep( c(0,1), 5)), 2)
tester = data.frame( cbind( edgesX[1,,], edgesX[2,,]  ) )
j = tester %>% group_by_all() %>% summarise(n = n()) %>% arrange(desc(n)) %>% ungroup() %>% mutate(prop = n/sum(n))

fp = apply(j[,-c(25,26)], 1, function(x) x %*% (x - truev) ) # FALSE POSITIVES
fn = apply(j[,-c(25,26)], 1, function(x) truev %*% (truev - x) ) # FALSE NEGATIVES

wp = as.numeric( fp %*% j$prop )
wn = as.numeric( fn %*% j$prop )

plotdf[m,] = c(nu_0[m], wp, wn)

}
# function to compute FP/ FN based on this

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
edge.index = 3 # index for edge indicator in posterior sample
ndraws = 1000
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

library(latex2exp)
j = data.frame(j)
ggplot(j, aes(j)) + geom_density() + theme_minimal() + 
  xlab(TeX('$\\phi$')) + ggtitle(TeX('Posterior Density of $\\phi$'))

#ggsave('post.png', dpi = 320)
# 0.25: -.57 0.75: -.54


##### TRUTH #####

y111.true = array(N, dim = c(N, ndraws))
X.pp.true = array(0, dim = c(2, N, ncov))

for(i in 1:ndraws) {
  
  X.pp.true[1, ,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[1,,] - 
     0.6 * cbind(rep(1,N), 0, 0) + u
  
  X.pp.true[2, ,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X.pp.true[1,,] - 
    0.3 * cbind(rep(1,N), 0, 0) + u
  
  
  Y = 0.8 * X.pp.true[2,,1] + rnorm(N, sd = 0.2) + u[,1] - 0.6 * rep(1,N)
  
  y111.true[,i] = Y
}

y000.true = array(N, dim = c(N, ndraws))
X.pp.true = array(0, dim = c(2, N, ncov))

for(i in 1:ndraws) {
  
  X.pp.true[1, ,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X[1,,] + u
  
  X.pp.true[2, ,] = rmvnorm(N, rep(0, ncov), (xdiffsd^2) * s) + 0.7 * X.pp.true[1,,] +  u
  
  Y = 0.8 * X.pp.true[2,,1] + rnorm(N, sd = 0.2) + u[,1] 
  
  y000.true[,i] = Y
}

j = colSums(y111.true - y000.true)/N
quantile(j)
hist(colSums(y111.true - y000.true)/N)



