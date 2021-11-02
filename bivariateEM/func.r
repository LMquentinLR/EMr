library(mvtnorm); library(pgmm); library(mvtnorm); library(mclust); library(ggplot2)

expectationStep <- function(X, mu, Sigma) {
  ######
  #### Expectation step of a bivariate dataset with missing data
  ######
  
  n=nrow(X)
  missing_idx.x1 = which(is.na(X[,1]))
  missing_idx.x2 = which(is.na(X[,2]))
  
  #all the element in X1 are observed
  s1_vec = rep(0,n)
  s11_vec = rep(0,n)
  s2_vec = rep(0,n)
  s22_vec = rep(0,n)
  
  # For observed elements in X1
  s1_vec[setdiff(1:n, missing_idx.x1)] = X[setdiff(1:n, missing_idx.x1), 1]
  s11_vec[setdiff(1:n, missing_idx.x1)] = X[setdiff(1:n, missing_idx.x1), 1]^2
  
  # For observed elements in X2
  s2_vec[setdiff(1:n, missing_idx.x2)] = X[setdiff(1:n, missing_idx.x2), 2]
  s22_vec[setdiff(1:n, missing_idx.x2)] = X[setdiff(1:n, missing_idx.x2), 2]^2
  
  # For missing elements in X1
  s1_vec[missing_idx.x1] = mu[1] + (Sigma[1,2]/Sigma[2,2]) * (X[missing_idx.x1,2]-mu[2])
  s11_vec[missing_idx.x1] = s1_vec[missing_idx.x1]^2 + Sigma[1,1] - Sigma[2,1]^2 / Sigma[2,2]
  
  # For missing elements in X2
  s2_vec[missing_idx.x2] = mu[2] + (Sigma[2,1]/Sigma[1,1]) * (X[missing_idx.x2,1]-mu[1])
  s22_vec[missing_idx.x2] = s2_vec[missing_idx.x2]^2 + Sigma[2,2] - Sigma[1,2]^2 / Sigma[1,1]
  
  s12_vec = s1_vec * s2_vec
  
  return(list(
    n=n,
    s1=sum(s1_vec), 
    s2=sum(s2_vec), 
    s11=sum(s11_vec), 
    s22=sum(s22_vec), 
    s12 = sum(s12_vec)
  )
  )
}

maximizationStep <- function(n, s1, s2, s11, s22, s12) {
  ######
  #### Maximization step of a bivariate dataset with missing data
  ######
  mu1 = s1/n
  mu2 = s2/n
  sigma11 = s11/n-(mu1)^2
  sigma22 = s22/n-(mu2)^2
  sigma12 = s12/n-(mu1*mu2)
  return(list(
    mu = matrix(c(mu1, mu2), nrow=1, ncol=2),
    sigma = matrix(c(sigma11, sigma12, sigma12, sigma22), nrow=2, ncol=2)
  ))
} 

emAlgotithm <- function(X, initMu, initSigma, maxIteration=50) {
  ######
  #### Runs the EM algorithm for a given number of iterations
  ######
  
  # Copies the initialized values for mu and sigma
  copyMu = initMu
  copySigma = initSigma
  
  for (loop in 1:maxIteration) {
    e = expectationStep(X, copyMu, copySigma)
    m = maximizationStep(e$n, e$s1, e$s2, e$s11, e$s22, e$s12)
    copyMu = m$mu
    copySigma = m$sigma
  }
  
  return(list(mu=copyMu, sigma=copySigma))
  
}