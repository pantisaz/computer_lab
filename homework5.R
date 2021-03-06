# BEGIN CODE
# Linear Regression function
lin.reg <- function(y, X){
  
  return( solve( t(X)%*%X , t(X)%*%y ))
}
# Ridge Regression function
ridge.reg <- function(y, X, lambda){
  
  return( solve( t(X)%*%X+lambda*diag(ncol(X)) , t(X)%*%y )) 
}
# Lasso Regression function
lasso.reg <- function(y, X, lambda){
  
  max.iter <- 10; 
  P <- ncol(X);
  beta <- solve(t(X)%*%X,t(X)%*%y)  
  beta.prev <- beta
  for( iter in 1:max.iter ){
   for( i in 1:P ){
       y.aux <- y-X[,setdiff(1:P,i)]%*%beta[setdiff(1:P,i)] 
       x.aux <- X[,i]
       cov <- sum( y.aux*x.aux ) 
       var <- sum( x.aux*x.aux )
       beta[i] <- sign(cov/var)*max( c( abs(cov/var) - lambda/(2*var) , 0 )) 
   }
# Calculating Beta
    if( sum( (beta-beta.prev)**2 ) < 1e-6 ){ return(beta) }
    
     beta.prev <- beta 
     }
  
  return(beta) 
}
# Cross Validation function
cross.validation <- function (y, X, lambda, pen.reg){
  
  data <- cbind(y,X)
  folds <- cut(seq(1,nrow(data)),breaks=5,labels=FALSE)
  RSS <- rep(0,5)

  for(i in 1:5){
    # Segement data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- data[testIndexes, ]
    trainData <- data[-testIndexes, ]
    
    beta <- pen.reg(trainData[,1], trainData[,2:ncol(trainData)], lambda)
    
# Calculate RSS   
    y.test <- testData[,1]
    X.test <- testData[,2:ncol(testData)]
      
    RSS[i] <- t((y.test - X.test %*% beta)) %*% (y.test - X.test %*% beta)
    
  }
  return(mean (RSS))
  
}
# fin
# END CODE