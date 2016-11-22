# BEGIN CODE

# PART 1
congruential.generator <- function(n, seed) {

m <- 2^31 - 1
a <- 48271 
c <- 1
x <- rep(0, n+1)
x[1] <- seed
# For loop for congruential generator
 for (s in 2:(n+1)){
  x[s] <- (a*x[s-1]+c) %% m 
 }
u <- x[2:(n+1)] / m
return(as.vector(u))
}
# PART 2
normal.simulator <- function(u, mu, sigma){
z <- rep(0,length(u))
# For loop for normal simulator  
 for (i in seq(1,length(u)-1,by=2)) {
  z[i] <- sqrt(-2*log(u[i+1]))*cos(2*pi * u[i])
  z[i+1] <- sqrt(-2*log(u[i+1]))*sin(2*pi * u[i])
  
 }
return(z*sigma + mu)
}
# PART 3
beta.simulator <- function(n,alpha,beta,K){
  
  y <- rep(0,n)
# For loop for beta simulator
  for( s in 1:n ){
    repeat{
      x     <- runif(1)
      u     <- runif(1)
      bool  <- u < dbeta(x,alpha,beta)/(K*dunif(x))
      if( bool ){
        y[s] <- x
        break
      }
    }
  }
  return(y)
}
# PART 4
# Pareto simulator, no for loop needed
pareto.simulator <- function(u,alpha) {
return((1-u)^(-(1/alpha)))
}

# END CODE





  