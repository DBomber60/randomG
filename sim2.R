# simulation analysis for > 2 treatment times
# let's start out with 2 treatment times and get things working

library(mvtnorm)


##### SIMULATE DATA #####
N = 1000 # number of subjects
ncov = 3 # 3 covariates

u = rnorm(N, sd = 0.5)
s <- matrix(0.25, ncov, ncov) # sigma X
diag(s) <- 1.0 #
xdiffsd = 0.25

su <- matrix(0.05, ncov, ncov) # sigma U
diag(su) <- 0.1
u = rmvnorm(N, mean = rep(0, ncov), sigma = su)

A_0 = rbinom(N, 1, prob = 0.5) # suppose the first treatment is (unconditionally) randomized

X_1 = rmvnorm(ntot, rep(0, ncov), (xdiffsd^2) * s) + 0.6 * cbind(-A_0, 0, 0) + u # treatment has effect on only the first column 

A_1 = rbinom(N, 1, prob = plogis( 1/3 * rowSums(X_1) )  ) # suppose the first treatment is (unconditionally) randomized

Y = 0.8 * X_1[,1] * rnorm(N) + u[,1] - 0.6 * A_1
#X_2 = rmvnorm(ntot, rep(0, ncov), (xdiffsd^2) * s) + 0.8 * X_1 + 0.6 * cbind(-A_0, 0, 0) + u

# for each covariate, model the inclusion


##### SAMPLER #####
