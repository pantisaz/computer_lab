# BEGIN CODE
# Begin cholesky decomposition function
my.chol <- function(A) {
  l <- matrix(0, nrow(A), nrow(A))
  l[1, 1] <- sqrt(A[1, 1]) 
# Calculate elements  
  for (i in 2:nrow(A)) {
    for (j in 1:i-1) {
      l[i, j] <- 1/l[j, j] * (A[i, j] - sum(l[i, 1:j-1] * l[j, 1:j-1]))
      l[i, i] <- sqrt(A[i, i] - sum(l[i, 1:i-1]**2))
    }
  }
return(l)
}
# Begin forward solve function
my.forward.solve <- function(L, b){
  x <- rep(0, nrow(L))
  x[1] <- b[1]/L[1,1] 
 # Calculate x 
  for (i in 2:nrow(L)) {
    x[i] <- (b[i] - sum(x[1:(i-1)]*L[i,1:(i-1)])) / L[i,i]
  }
return(x)
}
# Begin back solve function
my.back.solve <- function(U, b){
  n <- nrow(U)
  x <- rep(0, n)
  
  x[n] <- b[n]/U[n,n] 
# Calculate x 
  for (i in (n-1):1) {
    x[i] <- (b[i] - sum(U[i, (i+1):n] * x[(i+1):n])) / U[i,i]
  }
return(x)
}
# Begin Solve function incorporating other functions
my.solve <- function(A, b) {
  
  L <- my.chol(A)
  y <- my.forward.solve(L, b)
  x <- my.back.solve(t(L), y)
  
return(x)
}
# fin
# END CODE