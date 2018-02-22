sample.once.naive <- function(n=2000){
  m <- n
  zeta <- 10 
  x.old1 <- runif(m,-1,1)
  x.old2 <- runif(m,-1,1)
  x.old <- cbind(x.old1, x.old2)
  
  v <- function(x){
    ans <- 0
    if (x[1] > 0 & x[2] > 0 & x[1] < 1 & x[2] < 1){ # x is in the first quadrant
      ans <- rnorm(1,10,1)}
    return(ans)
  }
  
  y <- unlist(lapply(1:m, function(x) v(x.old[x,]))) # these are the sampled V's. 
  sum(y == 0)/length(y)  # This should go to 0.75 as m goes to infinity.
  
  s <- ifelse(y>zeta, 1, 0)
  return(mean(s))
}