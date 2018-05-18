n = 2000
noise.dim = 1
noise.var = 0.002 
c = 1
m <- floor(c * floor (n^(2/3))) # first-stage samples
theta <- runif(m, 0, 2 * pi)
u <- cbind(cos(theta), sin(theta))
epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
x <- u + epsilon
x <- cbind(x,  mvrnorm(m, 
                       mu = rep(0, noise.dim), 
                       Sigma =  noise.var * diag(noise.dim)))
x.old <- x

v <- function(x) {
  return(rnorm(1, exp(x[1] + x[2]), 1))
}

## ------------------------------------------------------------------------
y <- unlist(lapply(1:m, function(x) v(x.old[x, ]))) # these are the sampled V's.
# sum(y == 0)/length(y)

## ------------------------------------------------------------------------
zeta <- 1.35 # this is chosen to be the median of one trial y. Goal is this represents P(Y>zeta)= 0.5.
s <- ifelse(y > zeta, 1, 0)

## ------------------------------------------------------------------------
size <- floor(sqrt(m))  # number of k-means clusters
# model <- kmeans(x.old[,1:2], size)
model <- kmeans(x.old, size, iter.max=5)
model1 <- kmeans(x.old, size, iter.max=5)
arr <- model$cluster
arr1 <- model1$cluster
model.2d <- kmeans(x.old[,1:2], size)
for (i in 1:size){
  print(i)
  print(which(arr1==i))
  print(which(arr==i))
}

