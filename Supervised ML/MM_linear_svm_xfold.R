###############################################################################
### SVM Implementation with k-fold x-validation
# Implementation of Majorization-Minorization algorithm                   
# for the linear support vector machine with quadratic hinge                
# See Groenen et al (2008) for details:                                     
# http://link.springer.com/article/10.1007/s11634-008-0020-9                
###############################################################################

#### Settings #####
# set path to data set
map_data <- paste0(getwd(), "/Data/")

###############################################################################
#### SVM functions #####

svm.loss <- function(v, y, q, lambda){
  # quadratic hinge loss function
  L <- sum(ifelse(y==-1, pmax(0, q+1)^2, pmax(0, 1-q)^2))+lambda*t(v[-1])%*%v[-1]
  return(L)
}


svm.pupdate <- function(q, y, Z){
  # compute b for quadratic hinge
  b <- ifelse(y==-1, ifelse(q<=-1, q, -1), ifelse(q<=1, 1, q))  
  #compute update 
  v <- Z%*%b
  return(v)
}


svm.mm <- function(y, X, lambda=1, epsilon = 10^(-8), verbose = FALSE){ #use quadratic hinge
  # iterative MM algorithm, returns estimated (c,w) and loss for a given lambda
  # uses quadratic loss
  # calls svm.loss and svm.pupdate
  
  # check
  if (lambda <= 0) stop("Lambda has to be > 0")
  
  n <- length(y)
  if (!(sum(X[, 1])==n & var(X[, 1])==0)){ #add intercept column if it does not exist
    X <- cbind(1,X)
  } 
  p <- ncol(X)-1
  
  v <- rnorm(p+1, 10, 1)      #take some initial guess
  
  #create Z matrix
  P <- diag(p+1)
  P[1,1] <- 0
  Z <- solve(t(X)%*%X+lambda*P)%*%t(X)
  
  q <- X%*%v                 #make prediction with initial guess
  v <- svm.pupdate(q, y, Z)  #update parameters and store them
  q <- X%*%v                 #update the prediction with new v
  
  L <- numeric()
  L[2] <- svm.loss(v, y, q, lambda) #new loss 
  L[1] <- Inf                #old loss
  
  k <- 1
  while (k==1 | ((L[1]-L[2])/L[1])>epsilon) {  #
    k <- k+1
    L[1] <- L[2]             #old loss
    
    v <-  svm.pupdate(q, y, Z) #update parameters and store them
    q <- X%*%v               #update predictions with new v
    
    L[2] <- svm.loss(v, y, q, lambda) #new loss
    
    #for checks
    if (verbose) print(paste("k=", k, "L(beta(k))=", L[2], "L(beta(k-1))-L(beta(k))=", L[1]-L[2]))
    
  }
  c <- v[1]
  w <- v[-1]
  result <- list(c=c, w=w, Loss=L[2])
  return(result)
}

svm.predict <- function(X_test, c, w){ 
  n <- nrow(X_test)
  if (!(sum(X_test[, 1])==n & var(X_test[, 1])==0)){ #add intercept column if it does not exist
    X_test <- cbind(1,X_test)
  } 
  v <- c(c, w)
  q <- X_test%*%v
  y_h <- ifelse(q > 0, 1, -1) 
  result <- list(q=q, y_h=y_h)
  return(result)
}

#### K-fold x-validation #### 
svm.xvalid  <- function(y, X, k, lambdas){
  # necessary inputs
  # y = vector of dependent variable 
  # X = covariates
  # k = folds to use for k-vold x-validation
  # lambda = tuning parameter ranges  
  
  # create a grid of values
  grid <- expand.grid(
    lambdas = lambdas, 
    err = 0
  )
  
  folds <- cut(1:nrow(X), breaks=k, labels=1:k)   # identify k-folds
  indx <- sample(nrow(X))                         # shuffle indices of X
  
  for (t in 1:k){                                    # loop through folds
    test_indx <- indx[which(folds==t)]               # find indices of testing set 
    for (i in 1:nrow(grid)){                         # loop through parameters in grid
      # fit model with training fold
      svm_fit <- svm.mm(y_train, X_train, lambda=grid$lambdas[i], verbose = FALSE)
      
      # test model with testing fold
      svm_pred <- svm.predict(X_test = X[test_indx,], c = svm_fit$c, w = svm_fit$w)     
    
      # add error to grid
      grid$err[i] <- grid$err[i]+mean((svm_pred$y_h >= 0) != ifelse(y[test_indx]==1, TRUE, FALSE))
    }
  }
  grid$err <- grid$err/k          # average error across k-folds
  indx.min <- which.min(grid$err) # find minimum error entry index
  result <- list(lambda.min=grid$lambdas[indx.min], 
                 err.min = grid$err[indx.min],
                 table = grid)
  return(result)
}
###############################################################################
#### Example: heart data ##### 

# load data
heart <- read.csv(file = paste0(map_data, "heart.csv"),
                  header = TRUE, sep = ",")

# clean data
heart <- na.omit(heart)
y <- ifelse(heart$AHD=="Yes", 1, -1)
X <- model.matrix(AHD ~ ., data = heart)

# extract separate training and testing sets
set.seed(44)
train_n <- 207
indx <- sample(nrow(X), size = train_n)
X_train <- X[indx,]                       
y_train <- y[indx]                        
X_test  <- X[-indx,]                     
y_test  <- y[-indx]  

##### Compare function outputs #####
# run own function 
rslt_mm <- svm.mm(y_train, X_train, lambda=100, verbose = TRUE)
rslt_mm
# check correct classifications on test set
correct_100 <- mean(svm.predict(X_test, rslt_mm$c, rslt_mm$w)$y_h==y_test)  
correct_100

# compare to result of function in package SVMMaj
if (!require("SVMMaj")) install.packages("SVMMaj")
library("SVMMaj") 
rlst_svmmaj <- svmmaj(X_train, y_train, lambda=100, scale="none", hinge="quadratic")
# check correct classifications on test set
mean(ifelse(predict(rlst_svmmaj, X.new=X_test)>=0, 1, -1)==y_test)

# Do the functions produce the same predictions?
all(svm.predict(X_test, rslt_mm$c, rslt_mm$w)$y_h == ifelse(predict(rlst_svmmaj, X.new=X_test)>=0, 1, -1))


##### Using k-fold x-validation #### 
# Can we achieve better classifications using kfold x-validation? 
## train model over lambda parameter range using 10-fold x-validation
rslt_xvalid <- svm.xvalid(y_train, X_train, k=10, lambdas = seq(10, 1000, 10))
## run model with optimal lambda
rslt_mm <- svm.mm(y_train, X_train, lambda=rslt_xvalid$lambda.min, verbose = FALSE)
# store classifications on test set
correct_xvalid <- mean(svm.predict(X_test, rslt_mm$c, rslt_mm$w)$y_h==y_test)  

# message
print(paste("Optimal lambda is", rslt_xvalid$lambda.min, "with classification acccuracy on test set of", correct_xvalid)) 
print(paste("Classification accuracy on test set with lambda = 100 is", correct_100))



