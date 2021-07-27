library(tidyverse)
set.seed(2)
sig.u = 0.5
sig.y = 0.2

# confounding by l
n= 1000
A0 = rnorm(n, mean = 5)
u = rnorm(n, sd = sig.u ) # underlying health status
L = rbinom(n, 1, prob = plogis(3 * u + 1.5 * (A0 - 5))) # assume that high is good
A1 = rbinom(n, 1, prob = plogis(1 - 2 * L)) # higher p of treatment if L is 0
Y = rnorm(n, mean = 0.7 * u, sd = sig.y)

df = data.frame(A0, L, A1, Y)
summary(lm(Y ~ A0 + A1 + L, data = df))
sig.u^2 / (sqrt(sig.u^2 * (sig.y^2 + sig.u^2)))
cor(u, Y)
ggplot(df, aes(x=A0, y = Y)) + geom_point(aes(color=as.factor(L)), alpha = 0.5) + 
  geom_smooth(aes(color = as.factor(L)), method = lm) + 
  theme_classic()
