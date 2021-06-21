# Gibbs sampler for new parameter values in the case of a normal response
# input: old parameter values (theta, gamma, beta, sigsq), data (design, response)
# output: list of new sampled values for each parameter

sample.new = function(theta.old, gamma.old, beta.old, sigsq.old, design, resp) {
  
  p = ncol(design)
  # sample new theta
  # assumption: U(0,1) prior on theta/ binomial distribution on gamma (beta - binomial)
  
  theta.new = rbeta(1, shape1 = sum(gamma.old) + 1  ,shape2 = p - sum(gamma.old) + 1)
  
  # sample new gamma 
  # assumption: independent elements, fixed nu_0, nu_1
 
  p1 = dnorm(beta.old, sd = sqrt(nu_1)) * theta.new
  p0 = dnorm(beta.old, sd = sqrt(nu_0)) * (1-theta.new)
  gamma.new = rbinom(p, 1, prob = p1/(p1 + p0))  

  # sample new beta
  # assume scale mixture prior and normal likelihood for response with sigma = sigma^2 * I
  M = (t(design) %*% design)*sigsq.old^-1 + diag(ifelse(gamma.new == 1, 1/nu_1, 1/nu_0)) 
  Minv = solve(M)
  meanvec = sigsq.old^-1 * Minv %*% t(design) %*% resp
  beta.new = rmvnorm(1, mean = meanvec, sigma = Minv)
  
  # sample new sigsq
  # assumptions: IG(.5, .5) prior on sigsq
  ss = sum ( (resp - design %*% array(beta.new, dim=p) )^2 )
  n = length(resp)
  sigsq.new = rgamma(1, (n+1)/2, (ss+1)/2)
  
  
  return(list(theta.new=theta.new, 
              gamma.new=gamma.new, 
              beta.new=beta.new,
              sigsq.new=sigsq.new))
  
  
}