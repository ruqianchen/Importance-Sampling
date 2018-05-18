circle.density <- function(x.tmp = c(0,0), n.cut = 100){
  theta.tmp <- seq(0, 2*pi, by = 2*pi/n.cut)   
  u.tmp <- cbind(cos(theta.tmp), sin(theta.tmp))
  d.tmp <- sum(unlist(lapply(1:nrow(u.tmp), function(x){dmvnorm(x.tmp, mean = u.tmp[x,], sigma = 0.01 * diag(2))})))/n.cut
  return(d.tmp)
}

circle.density.df <- function(x.tmp = c(0,0), n.cut = 100){
  theta.tmp <- seq(0, 2*pi, length.out = n.cut + 1)   
  u.tmp <- cbind(cos(theta.tmp), sin(theta.tmp))
  d.tmp <- dmvnorm(x.tmp, u.tmp[1,], 0.01 * diag(2))
  for (j in 2:n.cut){
    d.tmp <- d.tmp + dmvnorm(x.tmp, u.tmp[j,], 0.01 * diag(2))
  }
  d.tmp <- d.tmp / n.cut
  return(d.tmp)
}

zeta <- 1.35 # this is chosen to be the median of one trial y. Goal is this represents P(Y>zeta)= 0.5.
v <- function(x) {
  return(rnorm(1, exp(x[1] + x[2]), 1)) 
}

sample.once.importance <- function(n = 2000, noise.dim = 1, noise.var = 0.000001, c = 1){
  out <- tryCatch(
    {
      m <- floor(c * floor (n^(2/3))) # first-stage samples
      # print(m) 
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
      model <- kmeans(x.old, size)
      cluster.size.count <-
        unlist(lapply(1:size, function(x) {
          length(which(model$cluster == x))
        }))
      
      ## ------------------------------------------------------------------------
      cluster.prob <-
        unlist(lapply(1:nrow(model$centers), function(x) {
          mean(s[which(model$cluster == x)])
        }))
      
      ## ------------------------------------------------------------------------
      df <- data.frame(model$centers)
      df$prob <- cluster.prob
      
      ## ------------------------------------------------------------------------
      sigma.new <- function(i) {
        cov(x.old[which(model$cluster == i), ])
      }
      
      ## ------------------------------------------------------------------------
      cluster.size.percentage <-
        unlist(lapply(1:size, function(x) {
          sum(s[which(model$cluster == x)]) / sum(s)
        }))
      
      ## ------------------------------------------------------------------------
      components <- 
        sample(1:size,
               prob = cluster.size.percentage,
               size = (n - m),
               replace = TRUE)
      
      ## ------------------------------------------------------------------------
      sigma.array <- lapply(1:size, function(x) unlist(sigma.new(x)))
      mean.array <- lapply(1:size, function(x) df[,1:(2 + noise.dim)][x, ])
      df$sigma <- sigma.array
      
      
      ## ------------------------------------------------------------------------
      x.new <- ldply( 
        components,
        .fun = function(x) {
          mvrnorm(1,
                  mu = unlist(mean.array[[x]]),
                  Sigma = sigma.array[[x]])
        }
      )
      
      y.new <- ldply(
        1:nrow(x.new),
        .fun = function(x)
          v(as.numeric(x.new[x, ]))
      )
      s.new <- ldply(
        1:nrow(x.new),
        .fun = function(x) {
          ifelse(y.new[x, ] > zeta, 1, 0)
        }
      )
      
      ## ------------------------------------------------------------------------
      p.old <- circle.density.df(as.matrix(x.new[,1:2]))
      if (noise.dim == 1){
        p.old <- p.old * dnorm(x.new[, (2+1):(2+noise.dim)], 0, 0.01)
      } else{
        p.old <- p.old * dmvnorm(x.new[,(2+1):(2+noise.dim)], rep(0, noise.dim), 0.01 * diag(noise.dim))
      }
      p.new <- cluster.size.percentage[1] * dmvnorm(x.new, unlist(mean.array[[1]]), sigma.array[[1]])
      for (j in 2:size) {
        p.new <- p.new + cluster.size.percentage[j] * dmvnorm(x.new, unlist(mean.array[[j]]), sigma.array[[j]])
      }
      
      ## ------------------------------------------------------------------------
      estimator.new <- mean(c(s.new[, 1] / p.new * p.old, s) )
      return(estimator.new)
    },
    error = function(cond){
      return(-999)
    }
  )
  return(out)
}


sample.once.naive <- function(n = 2000){
  m <- n
  ## ------------------------------------------------------------------------
  zeta <- 1.35
  m <- floor (n^(2/3)) # first-stage samples 
  theta <- runif(m, 0, 2 * pi)
  u <- cbind(cos(theta), sin(theta))
  epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
  x <- u + epsilon
  x.old <- x
  
  v <- function(x) {
    return(rnorm(1, exp(x[1] + x[2]), 1))
  }
  
  ## ------------------------------------------------------------------------
  y <- unlist(lapply(1:m, function(x) v(x.old[x,]))) # these are the sampled V's. 
  
  ## ------------------------------------------------------------------------
  s <- ifelse(y>zeta, 1, 0)
  return(mean(s))
}