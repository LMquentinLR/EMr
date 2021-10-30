# Library imports

library(pgmm); library(mvtnorm); library(mclust); library(ggplot2)
source(func.r)

# data imports

data(wine)
head(wine)

X = as.matrix(wine[,c(2,4)])
y = wine$Type

print(paste("Feature dimensions: ", paste(dim(X), collapse= ", ")))
print(paste("feature names: ", paste(colnames(X), collapse=", ")))

n = dim(X)[1]

training_split = round(2/3*n,0)

randomized_indexes = sample(1:nrow(X))
randomized_X = X[randomized_indexes,]
randomized_y = y[randomized_indexes]

X_train = randomized_X[1:training_split,]
X_test = randomized_X[training_split:n,]
y_train = randomized_y[1:training_split]
y_test = randomized_y[training_split:n]

# whole dataset EM + Kmeans

n_clusters = 3
min_clusters = length(unique(y))
max_clusters = 10

resultsEMwithoutKM = expectationMaximization(X, n_clusters)
resultsEMwithKM = expectationMaximization(X, n_clusters, kmInit=T)
resultsKM = kMeans(X, n_clusters)

# plots of results

plotData(X, resultsKM$closest_centroids, resultsKM$centroids)
plotData(X, resultsEMwithoutKM$clustering, resultsEMwithoutKM$means,
         t="EM (without K-Means init.)")
plotData(X, resultsEMwithKM$clustering, resultsEMwithKM$means, 
         t="EM (with K-Means init.)")
