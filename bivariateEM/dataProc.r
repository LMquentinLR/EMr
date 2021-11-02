# Library imports

library(mvtnorm); library(pgmm); library(mvtnorm); library(mclust); library(ggplot2)
source(func.r)

# distribution parameters

mu = matrix(c(5, -1), nrow=1, ncol=2)
sigma = matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2)

# data generation

n = 100
X = rmvnorm(n, mu, sigma)

# introduction of NAs

ratio = 0.3

missing_idx.mcar <- sample.int(n, ratio*n)
missing_idx.x1 <- sample(missing_idx.mcar, (ratio/2)*n)
missing_idx.x2 <- setdiff(missing_idx.mcar, missing_idx.x1)

X[missing_idx.x1, 1] = NA
X[missing_idx.x2, 2] = NA

# EM initialization - We use the empirical mean and sigma
mu_init = apply(X, 2, mean, na.rm=TRUE)
sigma_init = cov(X, use="complete.obs")

# Running the EM algorithm

emAlgotithm(X, mu, sigma)