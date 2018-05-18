
# we set noise.dim = 1 throughout
# we compute four estimators at once to fasciliate comparison while holding inputs constant for each trial
# first: importance sampling, with k means using 3D
# second: importance sampling, with k means using 2D structural
# third: naive, with sample from 3D, eval using 3D
# fourth: naive, with sample from 3D, eval using 2D



sample.once <- function(n = 2000, noise.var = 0.002, c = 1){
  out <- tryCatch(
    {
      noise.dim <- 1

      # CASE 3 NAIVE
      #################################################################################

      m <- n
      # m <- floor(c * (n^(2/3))) # first-stage samples
      # print(m) 
      theta <- runif(m, 0, 2 * pi)
      u <- cbind(cos(theta), sin(theta))
      epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
      x <- u + epsilon
      x <- cbind(x,  mvrnorm(m, 
                             mu = rep(0, noise.dim), 
                             Sigma =  noise.var * diag(noise.dim)))
      x.old <- x
      ## ------------------------------------------------------------------------
      y <- unlist(lapply(1:m, function(x) v(x.old[x,]))) # these are the sampled V's. 
  
      ## ------------------------------------------------------------------------
      s <- ifelse(y>zeta, 1, 0)
      estimator.new3 <- mean(s)


      # CASE 4 NAIVE
      #################################################################################

      ## ------------------------------------------------------------------------
      y <- unlist(lapply(1:m, function(x) v(x.old[x,1:2]))) # these are the sampled V's. 
  
      ## ------------------------------------------------------------------------
      s <- ifelse(y>zeta, 1, 0)
      estimator.new4 <- mean(s)

      # IMPORTANCE
      #################################################################################

      m <- floor(c * (n^(2/3))) 
      x.old <- x.old[1:m, ]
      

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
      
      # CASE 1
      #################################################################################

      ## ------------------------------------------------------------------------
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
      estimator.new1 <- mean(c(s.new[, 1] / p.new * p.old, s) )

      # CASE 2
      #################################################################################

      ## ------------------------------------------------------------------------
      model <- kmeans(x.old[,1:2], size)
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
        p.old <- p.old * dnorm(x.new[, (2+1):(2+noise.dim)], 0, sqrt(noise.var))
      } else{
        p.old <- p.old * dmvnorm(x.new[,(2+1):(2+noise.dim)], rep(0, noise.dim), noise.var * diag(noise.dim))
      }
      p.new <- cluster.size.percentage[1] * dmvnorm(x.new, unlist(mean.array[[1]]), sigma.array[[1]])
      for (j in 2:size) {
        p.new <- p.new + cluster.size.percentage[j] * dmvnorm(x.new, unlist(mean.array[[j]]), sigma.array[[j]])
      }
      
      ## ------------------------------------------------------------------------
      estimator.new2 <- mean(c(s.new[, 1] / p.new * p.old, s) )

      #################################################################################



      return(c(estimator.new1, estimator.new2, estimator.new3, estimator.new4))
    },
    error = function(cond){
      return(-999)
    }
  )
  return(out)
}