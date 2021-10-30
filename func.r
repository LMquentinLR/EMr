# install.packages("pgmm")
# install.packages('mvtnorm')
# install.packages('mclust')
# install.packages('ggplot2')

library(pgmm); library(mvtnorm); library(mclust); library(ggplot2)


euclidianDistance <- function(x, centroid){
  ######
  #### Computes the euclidian distance between a datapoint and a centroid.
  ######
  
  # Declares necessary functions for computing the euclidian distance
  dist = function(a, b) (a-b)^2
  norm = function(a)    sqrt(sum(mapply(dist, a, centroid)))
  
  # Computes the distances
  distances = t(apply(x, 1, norm))
  
  # Return statement
  return(distances)
}

findClosestCentroids <- function(x, centroids, clusters){
  ######
  #### Computes the centroid closest to each datapoints in an input array.
  ######
  
  # Variable initialization
  n = dim(x)[1]
  distances = c()
  
  # Computes the euclidian distances
  for (cl in 1:clusters){
    cluster_dist = euclidianDistance(x, centroids[cl,])
    distances = append(distances, cluster_dist)
  }
  distances = array(distances, dim=c(n, clusters))
  
  # Finds the closest centroids
  closest = (distances == apply(distances, 1, min))
  
  # Return statement
  return(list("distances"=distances, "closest"=closest))
    
}

findClusters <- function(likelihoods) {
  ######
  #### Given an array of likelihood values (each column representing the 
  #### likelihood of belonging to a specific cluster, each row representing 
  #### a datapoint in an underlying dataset), finds the most matching cluster.
  ######
  
  # Finds the most likely distribution to which each point belongs
  clustering = apply((likelihoods == apply(likelihoods, 1, max)), 1, which)
  
  # Return statement
  return(clustering)
}

generateCentroid <- function(x) {
  ######
  #### Generates a random centroid for a given dataset, given a normal
  #### multivariate distribution with the mean and variance equal to the
  #### empirical mu and Sigma of the dataset.
  ######
  
  # Variable initialization
  centroid = c()
  
  # Randomly draws a value for each dimension of the dataset
  centroid = rmvnorm(1, colMeans(x), cov(x))
  colnames(centroid) <- NULL
  
  # return statement
  return(centroid)
}

logSumExp <- function (x) {
  ######
  #### Computes the log-sum-exponential of a vector/list of variables.
  ######
  
  max(x) + log(sum(exp(x - max(x))))
}

logLikelihood <- function(n, clusters, x, prop, mu, sigma) {
  ######
  #### Computes the gamma values of a dataset as part of an EM algorithm.
  ######
  
  # Variable initialization
  lgn = matrix(nrow=n, ncol=clusters)
  
  # Computes the gammas per clusters
  for (cl in 1:clusters){
    log_proba = matrix(dmvnorm(x, mu[cl,], sigma[[cl]], log=TRUE))
    lgn[,cl] = log(prop[cl]) + log_proba
  }
  
  # return statement
  return(lgn)
}

mostMatchingDistribution <- function(clusters, x, prop, mu, sigma){
  ######
  #### Given a set of points, find the most matching distributions for them.
  #### To be used for a test set.
  ######
  
  # Computes the likelihood of belonging to a distribution
  likelihoods = logLikelihood(dim(x)[1], clusters, x, prop, mu, sigma)
  clustering = findClusters(likelihoods)
  
  # return statement
  return(list("likelihoods" = likelihoods, "clustering" = clustering))
}

kMeans <- function(x, clusters, maxIterations=200, printMessage=T) {
  ######
  #### Implementation of a K-Means algorithm.
  ######
  
  if(printMessage){cat("Processing KM-Means with", clusters, "cluster(s).")}
  
  # Variable initialization: dimensions of the input data
  n = dim(x)[1]
  K = dim(x)[2]
  # Variable initialization: centroids
  centroids = t(lapply(1:clusters, function(i) generateCentroid(x)))
  centroids = matrix(unlist(centroids), ncol=K, byrow=TRUE)
  # Variables initialization: vector to record step values
  history_centroids = c(centroids)
  
  ###########################
  ###### KM loop BEGIN ###### 
  ###########################
  for (loop in 1:maxIterations){
    
    # Computes the current closest centroids for each datapoint
    closest_centroids = findClosestCentroids(x, centroids, clusters)
   
    # Updates the closest centroids
    centroids = c()
    for (cl in 1:clusters) {
      centroids = rbind(centroids, colMeans(x[closest_centroids$closest[,cl],]))
    }
    
    # Erases the matrix's columns and rows names inherited from the dataset
    colnames(centroids) <- NULL
    rownames(centroids) <- NULL
    
    # Records the loop's update
    history_centroids = append(history_centroids, centroids)
  }
  ###########################
  ####### KM loop END ####### 
  ###########################
  
  #return statement
  ret = list(
    "centroids" = centroids,
    "distances" = closest_centroids$distances,
    "closest_centroids" = apply(closest_centroids$closest, 1, which),
    "centroids_history" = history_centroids
   )
  return(ret)
}

expectationMaximization <- function(
  x, clusters, 
  maxIterations=200, kmInit = F,
  printMessage=T
) {
  ######
  #### Implementation of a multivariate EM algorithm.
  ######
  
  if (printMessage) {
    if (kmInit) {
      cat("Processing EM with", clusters, "cluster(s), KM initialization.")
    } else {
      cat("Processing EM with", clusters, "cluster(s), random initialization.")
    }
  }
  
  # Variables initialization: dimensions of the input dataset
  n = dim(x)[1]
  K = dim(x)[2]
  # Variables initialization: starting values for the EM parameters (means
  # and sigmas are generated randomly with a multivariate normal distribution)
  prop = rep(1/clusters, clusters)
  if (kmInit) {
    mu = kMeans(x, clusters, maxIterations, printMessage=F)$centroids
  } else {
    mu = rmvnorm(clusters, colMeans(x), cov(x))
  }
  sigma = lapply(1:clusters, function(i) cov(x))
  # Variables initialization: comparators, vectors for recording 
  #                           step values, etc.
  history_LLH = c()
  history_prop = c()
  history_mu = c()
  history_sigma = c()
  previous_LLH = -Inf
  counter = 0
  loop_counter = 0
  
  ###########################
  ###### EM loop BEGIN ###### 
  ###########################
  for (loop in 1:maxIterations){
    
    ################################################
    #### E(xpectation)-Step (with logExp trick) ####
    ################################################
    
    loop_counter = loop_counter + 1
    
    # Computes log-likelihood of the setup at the start of the loop
    LLH = logLikelihood(n, clusters, x, prop, mu, sigma)
    sum_LLH = sum(apply(LLH, 1, logSumExp))
    
    # Checks the LLH with the last recorded LLH for early stopping
    # Criterion: no update for 5 loops (with llh rounded to 4 decimals)
    if (round(sum_LLH,4) <= round(previous_LLH,4)) {
      counter = counter + 1
      if (counter >= 5) {
        if (printMessage) {cat("\nEarly stopping at loop:", loop)}
        break
      }
    } else {counter = 0; previous_LLH = sum_LLH}
    
    # Records the current LLH (plotting purposes)
    history_LLH = append(history_LLH, sum_LLH)
    
    # Computes the gamma values per datapoints and clusters with
    # the logExp trick
    gam = exp(LLH - apply(LLH, 1, logSumExp))
    
    ################################################
    ############## M(aximization)-Step #############
    ################################################
    
    # Computes the update parameters of the EM algorithm 
    # for the current loop for each cluster
    for (cl in 1:clusters) {
      
      # Computes the sum of gammas and updates the cluster's proportions
      nk = sum(gam[,cl])
      prop[cl] = nk/n
      # Updates the mean parameters
      mean_compute = function(i) gam[i,cl] * x[i,]
      mu[cl,] = Reduce("+", lapply(1:n, mean_compute))/nk
      # Updates the covariance matrix parameters
      m = mu[cl,]
      sigma_compute = function(i) gam[i,cl] * (x[i,]-m) %*% t(x[i,]-m)
      sigma[[cl]] = Reduce("+", lapply(1:n, sigma_compute))/nk
      
    }
    
    # Records the loop's updates
    history_prop = append(history_prop, prop)
    history_mu = append(history_mu, mu)
    history_sigma = append(history_sigma, sigma)
    
  }
  ###########################
  ####### EM loop END ####### 
  ###########################
  
  # Computes the log-likelihood results for the given dataset
  llh_results = mostMatchingDistribution(clusters, x, prop, mu, sigma)
  llh_sum = sum(apply(llh_results$likelihoods, 1, logSumExp))
  if (printMessage) {cat("\nEnd log-likelihood: ", llh_sum, "\n")}
  
  # Formatting the record variables
  history_prop = array(history_prop, c(loop_counter, clusters))
  history_mu = array(history_mu, c(clusters, K, loop_counter))
  history_sigma = array(history_sigma, c(clusters, K, K, loop_counter))
  
  # return statement
  ret = list(
    "N" = n,
    "n_clusters" = clusters,
    "prop" = prop,
    "means" = mu,
    "sigma" = sigma,
    "clustering" = llh_results$clustering,
    "llh_per_points" = llh_results$likelihoods,
    "llh_sum" = llh_sum,
    "prop_history" = history_prop,
    "means_history" = history_mu,
    "sigma_history" = history_sigma,
    "llh_history" = history_LLH
    )
  return(ret)
}

akaikeIC <- function(llh, clusters, K) {
  ######
  #### Implementation of the the Akaike Information Criterion such that:
  #### AIC = log-likelihood - eta(M)
  #### i.e. the final log-likelihood of a model minus the number of free
  #### scalar parameters in the model (nb of proportions (-1 as there are
  #### only cluster-1 degrees of freedom) + nb of means + nb of sigmas).
  ######
  llh - (clusters-1) + clusters*K + clusters*((K*(K+1))/2)
}

computeAIC <- function(x, max_cluster=max_clusters, print_steps=TRUE){
  ######
  #### Computes the AIC of an EM algorithm implementation.
  ######
  
  # Variable initialization
  akaike_results = c()
  
  # Loops through a cluster parameter range to compute the AIC
  for (cl in min_clusters:max_cluster){
    EM = expectationMaximization(x, clusters=cl, printMessage=F)
    akaike = akaikeIC(EM$llh_sum, cl, dim(x)[2])
    akaike_results = append(akaike_results, akaike)
    if (print_steps) {
      print(paste("Total LLH with ", cl, " clusters: ", round(akaike, 3)))
    }
  }
  
  # Prints the result
  print(paste("The best AIC result is achieved with ", 
              which.max(akaike_results)+2, 
              " clusters."))
}

bayesianIC <- function(llh, clusters, n, K) {
  ######
  #### Implementation of the Bayesian Information Criterion such that:
  #### BIC = LLH - 1/2*eta(M)*log(n)
  #### i.e. the final log-likelihood of a model minus the half of the number
  #### of free scalar parameters in the model (nb of proportions (-1 as there
  #### are only cluster-1 degrees of freedom) + nb of means + nb of sigmas)
  #### then multiplied by the log of the number of datapoints.
  ######
  
  # Variable declaration
  llh - 1/2 * ((clusters-1) + clusters*K + clusters*((K*(K+1))/2)) * log(n)
}

computeBIC <- function(x, max_cluster=max_clusters, print_steps=TRUE){
  ######
  #### Computes the BIC of an EM algorithm implementation.
  ######
  
  # Variable declaration
  bayesian_results = c()
  
  # Loops through each cluster parameter to compute the corresponding AIC
  for (cl in min_clusters:max_cluster){
    EM = expectationMaximization(x, clusters=cl, printMessage=F)
    bayesian = bayesianIC(EM$llh_sum, cl, dim(x)[1], dim(x)[2])
    bayesian_results = append(bayesian_results, bayesian)
    if (print_steps) {
      print(paste("Total LLH with ", cl, " clusters: ", round(bayesian, 3)))
    }
  }

  # Prints the result
  print(paste("The best BIC result is achieved with ", 
              which.max(bayesian_results)+2, 
              " clusters."))
}

doubleCrossValidation <- function(x_train, x_test, folds=10, max_cl=max_clusters) {
  ######
  #### Implements a double cross-validation with the resulting log-likelihood
  #### being the selection criteria.
  ######
  
  # Variable initialization
  n_train = dim(x_train)[1]
  foldAllocation = ceiling(seq_along(c(1:n_train))/(n_train/10))
  fold_indexes = split(c(1:n_train), foldAllocation)
  mean_cluster_criteria = c()
  
  # Iterates over the cluster range
  for (cl in min_clusters:max_cl){
    
    # Performs the first step of the double cross-validation: iteration
    # over the folds of the training set
    llhs = c()
    best_model = NULL
    
    # iterates over the k-folds
    for (kFold in fold_indexes){
      
      # Stores the training dataset depending on which fold is validation
      x_train_train = x_train[-kFold,]
      x_train_val = x_train[kFold,]
      
      # Computes the EM (on the training set)
      EM = expectationMaximization(x_train_train, cl, printMessage=F)
      
      # Retrieves the resulting llhs/distributions on the validation set
      dists = mostMatchingDistribution(
        EM$n_clusters, x_train_val, EM$prop, EM$means, EM$sigma
      )
      
      # Records the resulting log-likelihood
      sum_of_llhs = sum(apply(dists$likelihoods, 1, logSumExp))
      llhs = append(llhs, sum_of_llhs)
      
      # Updates the best models computed so far if needed
      if (is.null(best_model) || sum_of_llhs == min(llhs)){
        best_model <- EM
      }
    }
    
    # Computes the likelihood on the test set given the model
    # with the highest performance on the validation set
    dists = mostMatchingDistribution(
        best_model$n_clusters, 
        x_test,
        best_model$prop, 
        best_model$means, 
        best_model$sigma
        )
    
    # Records the resulting log-likelihood
    llhs = append(llhs, sum(apply(dists$likelihoods, 1, logSumExp)))
    
    # Records the mean llhs achieved with the cluster parameter
    print(paste("Mean log-likelihood achieved with ", cl, " clusters: ",
                round(mean(llhs), 4)))
    mean_cluster_criteria = append(mean_cluster_criteria, mean(llhs))
  }
  
  # Finds which cluster had the best performance
  best_cluster = which.max(mean_cluster_criteria)
  print(paste("The best result is achieved with ", 
            best_cluster+2, 
            " clusters (double CV)."))
  
  # return statement
  return(mean_cluster_criteria)
  
}

plotData <- function(
  x, clustering, mean_clusters,
  t="K-Means", xl="Alcohol", yl="Fixed Acidity"
) {
  ######
  #### Plots a cloud of point with clustering ellipses.
  #####
  
  # Formats the input data as a dataframe
  df = data.frame(x)
  mean_clusters = data.frame(mean_clusters)
  names(mean_clusters) = names(df)
  
  # Formats the clustering labels as factors for coloring
  colors = as.factor(clustering)
  
  # Declares the plot
  p = ggplot(df, aes_string(names(df)[1], names(df)[2], color=colors)) + 
    geom_point() +
    stat_ellipse(geom="polygon", aes(fill=colors), alpha=0.05) + 
    guides(fill = "none") + 
    labs(color="Wine Type", 
         title=paste("Clustering obtained via", t), 
         x=xl,y=yl) + 
    geom_point(data=mean_clusters, color="black")
  
  # return statement 
  return(p)
}

