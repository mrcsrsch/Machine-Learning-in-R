####################################################################################
# Simple function to calculate OLS estimate and standard statistics         
# that can handle formulas with interaction terms                            
####################################################################################

####################################################################################
### function ###
OLS <- function(formula, Data.frame){
  # function calculates OLS estimate, std. error, t-stat, p-value, F-stat and multiple R-squared
  # it can also handle interaction terms
  # formula: a formula in R formula style (y ~ x) 
  # Data.frame: a data frame 
  
  vars <- as.vector(attr(terms(formula), "term.labels"))   # take independent variables out of formula
  vars <- c(as.vector(all.vars(formula))[1], vars)
  interact <- vars[grep(":", vars)]                        # extract interaction terms if there are any
  
  if (length(interact)!=0){                                # if there are interaction terms then create them 
    for (i in seq_len(length(interact))){                  # loop through interactions and create vars
      dummy_vector <- sapply(strsplit(interact[i], ":"), "[", 1)
      dummy_vector[2] <- sapply(strsplit(interact[i], ":"), "[", 2)
      Data.frame[, paste(dummy_vector, collapse=":")] <- Data.frame[, dummy_vector[1]] * Data.frame[, dummy_vector[2]]
    }
  }
  K <- data.matrix(Data.frame[, vars])                         # convert final data.frame to matrix and select relevant subset
  K <- na.omit(K)                                          # drop rows with missing observations 
  y <- K[, 1]                                              # get vector dep. var. y
  X <- cbind(1, K[, -1])                                   # get matrix indp. vars X and add constant 
  X_X <- solve(t(X)%*%X)                                   # store (X'X)^-1
  
  # calculate OLS point estimate
  b <- X_X%*%t(X)%*%y                                      #(X'X)^-1X'y
  
  # calculate OLS standard errors 
  y_h <- X%*%b
  resid <- y-y_h                                           # get residuals
  n <- nrow(X)                                             # no obs. 
  p <- ncol(X)                                             # no regressors
  s2 <- t(resid)%*%resid*(1/(n-p))                         # estimate sigma^2 
  std <- sqrt(diag(X_X)*as.vector(s2))                     # calculate std. deviations
  
  # t-values and corresponding p-values
  t <- b/std
  pv <- 2*pt(-abs(t), df=n-1)
  
  # R^2
  y <- y-mean(y)
  y_h <- y_h-mean(y_h)
  R2 <- (t(y_h)%*%y_h)/((t(y)%*%y))
  
  # F-statistic and p-value
  RSS1 <- t(resid)%*%resid
  resid <- y-mean(y)
  RSS0 <- t(resid)%*%resid
  Fstat <- ((RSS0-RSS1)/(p-1))/(RSS1/(n-p))
  pFstat <- pf(Fstat, p-1, n-p, lower.tail = FALSE)
  
  # Output table
  out <- cbind(b, std, t, pv)
  rownames(out) <- c("(Intercept)", vars[-1])
  colnames(out) <- c("coef.", "std. error", "t value", "Pr(>|t|)")
  
  # print output
  cat(paste("Formula: ", format(formula)), "\n")
  cat('Linear regression output:', '\n')
  print(out, digits=5)
  cat(paste0('F-stat: ', round(Fstat,5),' on ', p-1, ' and ', n-p, ' DF, p-value: ', round(pFstat,5), '\n'))
  cat(paste0('Multiple R-squared: ', round(R2,5)))
}

####################################################################################
### example: Boston data ###

# load data
data("Boston", package = "MASS")

# run own function 
OLS(medv ~ age + tax*indus, Boston)

# compare to built-in function
summary(lm(medv ~ age + tax*indus, Boston))