####################################################################################
### Implementation of a simple k-means cluster algorithm  ####
# see: Algorithm 14.1 in Hastie, Tibshirani & Friedman (2017): Elements of Statistical Learning  
####################################################################################

#### Settings #### 
# set path to data set
map_data <- paste0(getwd(), "/Data/")

####################################################################################
#### K-means function ####

K_means <- function(mat, K, initial = NA){
  # mat: numeric data matrix
  # K: desired number of clusters 
  # initial: (optional) initial allocation matrix
  
  # checks, transformations 
  if (!is.matrix(mat)) X <- as.matrix(mat) else X <- mat
  if (!is.numeric(X)) warnings("Please specifiy a numeric data matrix.")
  n <- nrow(X)
  
  # initialize empty matrices
  centroids <- matrix(NA, nrow=K, ncol=ncol(X))
  dists <- matrix(NA, nrow=n, ncol=K)
  
  # set initial assignment 
  old_assgn <- rep(0,n) # no assignment
  if (any(is.na(initial))) assgn <- sample(1:K,n,replace=T) else assgn <- initial
  
  while (any(assgn != old_assgn)){
    old_assgn <- assgn
    
    for (k in 1:K){
      # calculate cluster centroid
      centroids[k,] <- apply(X[which(assgn==k), , drop=FALSE], 2, mean)
      # calculate euclidean distance
      dists[,k] <- apply(X, 1, function(x) {
        dist(rbind(centroids[k,], x), method="euclidean", p=2)
      })
    }
    # assign to closest cluster
    assgn <- apply(dists, 1, which.min)
  }
  
  # assign names for output
  if (!is.null(rownames(X))) names(assgn) <- rownames(X) else names(assgn) <- paste0("Row", 1:nrow(X))
  rownames(centroids) <- paste(1:K) # each row is a cluster
  if (!is.null(colnames(X))) colnames(centroids) <- colnames(X) else colnames(centroids) <- paste0("Metric", 1:ncol(X))
  
  return(list(assignment=assgn, centroids=centroids)) # return list
}

####################################################################################
### Example: cityweather ####
# compare K_means function to built-in Kmeans function

# load data
load(paste0(map_data, "cityweather.RData"))

# run built-in function
set.seed(42)
builtin <- kmeans(cityweather, 4)

# get initial assignment from built-in K-means function
## note that K_means can find local minima as solutions (see HTF)
## --> use initial assignment of built-in function for comparison
set.seed(42)
initial <- kmeans(cityweather, 4, iter.max=1)$cluster

# run own clustering function using initial assignment
own <- K_means(cityweather, 4, initial)

# compare
table(own$assignment, builtin$cluster)
