#### Settings #####
# packages 
if (!require("MASS")) install.packages("MASS")
if (!require("plotrix")) install.packages("plotrix")
if (!require("glmnet")) install.packages("glmnet")
require("MASS")
require("plotrix")
require("glmnet")

# load data
## set path to data set
map_data <- paste0(get_wd(), "/Data/")
load(paste0(map_data, "supermarket1996.RData"))

# prepare data
supermarket1996 <- subset(supermarket1996, select=-c(STORE,CITY,ZIP,GROCCOUP_sum,SHPINDX))
X <- as.matrix(supermarket1996)
X <- X[, -1]  # remove grocery sum 
X <- scale(X) # standarize
y <- as.vector(supermarket1996$GROCERY_sum)
y <- scale(y)

##########################
### Functions #########
MM_elasticnet <- function(X, y, alpha=0.5, lambda=1, epsilon=10^(-8)){
  ## Function that performs elastic net estimation for 0<alpha<1, 
  ## lasso for alpha=1 and ridge for alpha=0 using a majorization-minorization algorithm
  ## Inputs: standardized independent vars matrix X, dependent var vector y,
  ## constraint mixing parameter alpha, constraint weight lambda, control term epsilon
  ## Output: vector estimated beta
  
  #no of observations and variables
  n <- nrow(X)
  p <- ncol(X)
  
  #as an initial guess for beta start with a vector of 1s
  beta <- rep(1, p)   
  
  #compute D matrix with first beta
  D <- diag(1/(abs(beta)))
  
  #compute A matrix with first beta; first part never changes, so store it
  A_1 <- n^(-1)*t(X)%*%X + lambda*(1-alpha)*diag(p)
  A <- A_1  + lambda*alpha*D
  
  #also store second part of beta 
  beta_a <- n^(-1)*t(X)%*%y
  
  #just some value to pass to the loop constriant 
  L <- 0
  
  #write the constraint output into position two (L(beta(k))) 
  c <- (2*n)^(-1)*(t(y)%*%y)+0.5*lambda*alpha*sum(abs(beta)) 
  L[2] <- 0.5*t(beta)%*%A%*%beta-n^(-1)*(t(beta)%*%t(X)%*%y)+c 
  
  # update beta
  k <- 1
  while (k==1 | ((L[1]-L[2])/L[1])>epsilon){
    k <- k+1
    
    #update D matrix with beta(k-1)
    D <- diag(1/(pmax(epsilon,abs(beta)))) #compute D with "old" beta
    
    #update A matrix
    A <- A_1 + D*lambda*alpha
    
    #update c and beta
    c <- (2*n)^(-1)*(t(y)%*%y)+0.5*lambda*alpha*sum(abs(beta)) 
    beta <- solve(A, beta_a) #compute new beta vector
    
    #L(beta(k)) becomes L(beta(k-1)) for next iteration
    L[1] <- L[2]
    
    #update L with current k
    L[2] <- 0.5*t(beta)%*%A%*%beta-n^(-1)*(t(beta)%*%t(X)%*%y)+c
    
    #for checks
    print(paste("k=", k, "L(beta(k))=", L[2], "L(beta(k-1))-L(beta(k))=", L[1]-L[2]))
  }
  
  #return final beta vector
  return(beta)
}

k_crossval_elasticnet <- function(X, y, k, alpha=0.5, lambda=seq(0.1,1,0.1), epsilon=10^(-8)){ 
  ## Function that performs k_crossvalidation 
  ## for the lambda parameter of the elastic net
  ## calls another function for performing the elastic net estimation 
  ## Inputs: standardized independent matrix X, vector y,
  ## folds k, constraint mixing parameter alpha, range of lambdas, control term epsilon
  ## Output: vector beta of optimal estimates 
  
  #reshuffle data and cut into k folds
  ##first rebind so nothing goes wrong with shuffling 
  K <- cbind(y,X)
  ##shuffle
  K <- K[sample.int(nrow(X)),]
  ##and take apart gain
  X <- K[, -1]
  y <- K[, 1]
  ##indentify k-folds
  k_folds <- cut(1:nrow(X), breaks=k, labels=1:k)
  
  #create some necessary vectors 
  RMSEs <- numeric()
  lambdas <- numeric()
  
  #loop through lambdas
  for (lam in lambda){
    MSE <- 0
    #perform k-fold cross validation
    for (t in 1:k){
      # find test/train folds and subset
      testing_fold <- which(k_folds==t)
      X_test <- X[testing_fold, ]
      y_test <- y[testing_fold]
      X_train <- X[-testing_fold, ]
      y_train <- y[-testing_fold]
      
      #perform elastic net
      beta <- MM_elasticnet(X=X_train, y=y_train, alpha=alpha, lambda=lam, epsilon=epsilon)
      #beta <- as.vector(glmnet(X_train, y_train,  alpha = alpha, lambda=lam)$beta)
      
      #add MSE 
      y_h <- X_test%*%beta
      MSE <- MSE+length(testing_fold)^(-1)*(t(y_test-y_h)%*%(y_test-y_h))
    }
    #add RMSE of cross validation for given lambda to vector, add lambda to a vector
    RMSEs <- append(RMSEs, sqrt(k^(-1)*MSE))
    lambdas <- append(lambdas, lam)
  }
  #find smallest RMSE
  position <- which.min(RMSEs)
  
  #do elastic net for optimal lambda
  beta <- MM_elasticnet(X, y, alpha=alpha, lambda=lambdas[position], epsilon=epsilon)
  #beta <- glmnet(X, y,  alpha = alpha, lambda=lambdas[position])$beta
  #beta <- as.vector(beta)
  
  #print optimal lambda and return optimal beta
  cat(paste("optimal lambda:", lambdas[position], "\n"))
  cat(paste("Minimum RMSE:", RMSEs[position], "\n"))
  return(beta)
}

####################
### Analysis ######

set.seed(42)

#compare k_crossval_elasticnet output 
result.k_cross <- k_crossval_elasticnet(X,y, alpha=0.5, k=11, lambda=10^seq(-10,5,length.out = 50))
b_kcross <- round(result.k_cross, digits=2)

#to cv.glmnet output
result.elasticnet.cv <- cv.glmnet(X, y, alpha = 0.5, 
                                  lambda=10^seq(-10,5,length.out = 50), nfolds = 11)  
print(result.elasticnet.cv$lambda.min)      # Best cross validated lambda
result.elasticnet.best <- glmnet(X, y, alpha = 0.5, 
                                 lambda = result.elasticnet.cv$lambda.min)   
b_cvglmnet <- round(result.elasticnet.best$beta, digits = 2)
sqrt(min(result.elasticnet.cv$cvm))

compare <- b_cvglmnet - b_kcross
mean(compare)