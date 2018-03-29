library(distr)
library(plyr)
library(mclust)
library(MASS)
library(ggplot2)
library(mvtnorm)

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

# z1 <- mvrnorm(m, mu = c( 0), Sigma = 0.01 * diag(1))
# z2 <- mvrnorm(m, mu = c( 0), Sigma = 0.01 * diag(1))
# z1.prob <- dmvnorm(cbind(z1), c(0), 0.01 * diag(1))
# z2.prob <- dmvnorm(cbind(z2), c(0), 0.01 * diag(1))
# z.prob <- dmvnorm(cbind(z1,z2), c(0,0), 0.01 * diag(2))

n <- 2000
m <- floor (n^(2/3)) # first-stage samples 
theta <- runif(m, 0, 2 * pi)
u <- cbind(cos(theta), sin(theta))
plot(u)

epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
x <- u + epsilon
plot(x)
x <- cbind(x,  mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2)))

x.old <- x

v <- function(x) {
  return(rnorm(1, exp(x[1] + x[2]), 1))
}

## ------------------------------------------------------------------------
y <-
  unlist(lapply(1:m, function(x)
    v(x.old[x, ]))) # these are the sampled V's.
sum(y == 0)/length(y)
mean(y) # 1.525384
median(y) # 1.377392

## ------------------------------------------------------------------------
zeta <- 1.35 # this is chosen to be the median of one trial y. Goal is this represents P(Y>zeta)= 0.5.
s <- ifelse(y > zeta, 1, 0)


# ## ----plot0---------------------------------------------------------------
color <- ifelse(y > zeta, "blue", "red")
plot(x.old,
     col = color,
     pch = 20,
     main = "First stage samples of x, color by values of y's.")
grid()
legend(
  'bottomright',
  c("y>10", "y<10"),
  pch = 20,
  col = c("blue", "red"),
  bty = 'y',
  cex = .65
)
# 
## ------------------------------------------------------------------------
size <- floor(sqrt(m))  # number of k-means clusters
model <- kmeans(x.old, size)
cluster.size.count <-
  unlist(lapply(1:size, function(x) {
    length(which(model$cluster == x))
  }))

# ## ----plot1--------------------------------------------------------------------
barplot(
  cluster.size.count,
  main = "Size of Clusters",
  xlab = "Cluster index",
  names.arg = 1:size,
  col = "purple"
)

# ## ----plot2.0---------------------------------------------------------------

plot(data.frame(x.old),
     col = model$cluster,
     # xlim = c(-1,1),
     # ylim = c(-1,1),
     pch = 20)


# ## ----plot2---------------------------------------------------------------

plot(x.old,
     col = model$cluster,
     # xlim = c(-1,1),
     # ylim = c(-1,1),
     pch = 20)
# points(model$centers, col = "black", pch = 2)
grid()
text(
  model$centers[, 1],
  model$centers[, 2],
  labels = as.character(1:size),
  col = "black",
  pos = c(1),
  offset = -0.16
)

# NOTE: THESE SHOULD LOOK MORE OR LESS UNIFORMLY DISTRIBUTED
plot(x.old[,2:3],
     col = model$cluster,
     # xlim = c(-1,1),
     # ylim = c(-1,1),
     pch = 20)
# points(model$centers, col = "black", pch = 2)
grid()
text(
  model$centers[, 1],
  model$centers[, 2],
  labels = as.character(1:size),
  col = "black",
  pos = c(1),
  offset = -0.16
)

# NOTE: THESE SHOULD LOOK MORE OR LESS UNIFORMLY DISTRIBUTED
plot(x.old[,3:4],
     col = model$cluster,
     # xlim = c(-1,1),
     # ylim = c(-1,1),
     pch = 20)
# points(model$centers, col = "black", pch = 2)
grid()
text(
  model$centers[, 1],
  model$centers[, 2],
  labels = as.character(1:size),
  col = "black",
  pos = c(1),
  offset = -0.16
)

## ------------------------------------------------------------------------
cluster.prob <-
  unlist(lapply(1:nrow(model$centers), function(x) {
    mean(s[which(model$cluster == x)])
  }))

# ## ----plot3---------------------------------------------------------------
barplot(
  cluster.prob,
  main = "Proportion of Y>zeta's in the Clusters",
  xlab = "Cluster index",
  names.arg = 1:size,
  col = "purple"
)

# ## ----plot4---------------------------------------------------------------
pos.prob <- which(cluster.prob > 0)
colors.grayscale <-
  paste("gray", floor((1 - cluster.prob) * 100), sep = "")
plot(
  x.old[which(model$cluster == pos.prob[1]), ],
  col = colors.grayscale[1],
  bg = colors.grayscale[1],
  xlim = c(-1.5, 1.5),
  ylim = c(-1.5, 1.5),
  main = "Clusters and Their Respective Probabilities",
  xlab = "x_1",
  ylab = "x_2",
  pch = 20
)
grid()
for (i in 2:size) {
  points(x.old[which(model$cluster == i), ],
         col = colors.grayscale[i],
         bg = colors.grayscale[i],
         pch = 20)
}
text(
  model$centers[, 1],
  model$centers[, 2],
  labels = paste(
    paste("#", as.character(1:size), ":", sep = ""),
    paste(floor(cluster.prob * 100), rep("%", size), sep = ""),
    sep = " "
  ),
  col = "red",
  pos = c(1),
  offset = -0.16
)
# text(model$centers[,1], model$centers[,2], labels=as.character(1:size), col="red", pos=c(1), offset=-0.16)
# text(model$centers[,1],
#      model$centers[,2],
#      labels=paste(floor(cluster.prob*100), rep("%", size), sep = ""),
#      col="red",
#      pos=c(4),
#      offset=1)


## ------------------------------------------------------------------------
df <- data.frame(model$centers)
# names(df) <- c("centers.old1", "centers.old2")
df$prob <- cluster.prob
# centers.new <- matrix(unlist(lapply(1:size,function(x){colMeans(x.old[which(model$cluster == x),])})),,2)
# df$centers.new1 <- centers.new[,1]
# df$centers.new2 <- centers.new[,2]
# df$centers.new2 <- df$centers.old2
# df$centers.new1 <- df$centers.old1
# df

## ------------------------------------------------------------------------
sigma.new <- function(i) {
  cov(x.old[which(model$cluster == i), ])
}

## ------------------------------------------------------------------------
cluster.size.percentage <-
  unlist(lapply(1:size, function(x) {
    sum(s[which(model$cluster == x)]) / sum(s)
  }))
# sum(cluster.size.percentage) # = 1

## ----plot5--------------------------------------------------------------------
barplot(
  cluster.size.percentage,
  main = "Probability of Clusters",
  xlab = "Cluster index",
  names.arg = 1:size,
  col = "purple"
)

## ------------------------------------------------------------------------
components <- 
  sample(1:size,
         prob = cluster.size.percentage,
         size = (n - m),
         replace = TRUE)

# ## ----plot6--------------------------------------------------------------------
# hist(
#   components,
#   col  = "purple",
#   breaks = c(0:(size + 1)),
#   prob = TRUE,
#   xaxt = 'n',
#   xlab = "Component Index"
# )
# # lines(density(components))
# axis(1,
#      at = seq(0.5, (size + 1.5), by = 1),
#      labels = c(1:(size + 2)),
#      las = 1)

## ------------------------------------------------------------------------
sigma.array <- lapply(1:size, function(x) unlist(sigma.new(x)))
mean.array <- lapply(1:size, function(x) df[,1:4][x, ])
df$sigma <- sigma.array
# df

# ## -----gmm.unused-------------------------------------------------------------------
# gmm.pdf <- function(x){ # goal: This function should integrate to 1.
#   sum(cluster.size.percentage *
#         unlist(lapply(1:size,
#                       function(j){dmvnorm(x,
#                                           mean = mean.array[[j]],
#                                           sigma = sigma.array[[j]]
#                       )
#                       }
#         )
#         )
#   )
# }

## ----plot7.visualize.gmm--------------------------------------------------------------------
plot.cluster <- function(counter){ # counter indicate cluster indices
  
  num.sample <- 300 # sample 300 times from each normal.
  
  this.cluster.points <- rmvnorm(
    num.sample, mean = unlist(mean.array[counter]),
    sigma = matrix(unlist(sigma.array[counter]),4,4) # CRQ: updated this
  ) # generate points from the normal
  
  df.tmp1 <- data.frame(this.cluster.points) # save the generated points in a data frame
  
  df.tmp1$prob <- unlist(lapply(1:num.sample, function(x)dmvnorm(this.cluster.points[x,], mean = unlist(mean.array[counter]), sigma = matrix(unlist(sigma.array[counter]),4,4)))) # compute the pdf of the normal at each sampled point
  # CRQ: updated this
  
  df.tmp2 <- data.frame(cbind(x.old[which(model$cluster==counter),],0)) # obtain the sampled points belonging to this cluster from stage-one sampling
  
  names(df.tmp2) <- names(df.tmp1)
  df.tmp1 <- rbind(df.tmp1, df.tmp2) # combine the two temporary data frames
  
  rbPal <- colorRampPalette(c('red','blue')) # specify color scale.
  
  df.tmp1$Col <- rbPal(nrow(df.tmp1))[as.numeric(cut(df.tmp1$prob, breaks = nrow(df.tmp1)))] # adds a column of color values based on the y values.
  
  plot(x.old,
       col = model$cluster,
       xlim = c(-1.2, 1.2),
       ylim = c(-1.2, 1.2),
       pch = model$cluster,
       main = paste("Sampling from Normal Distribution #", counter))
  grid()
  points(model$centers, col = 1:size, pch = 2)
  text(model$centers[,1], model$centers[,2], labels=as.character(1:size), col="dark orange", pos=c(1), offset=-0.16)
  
  points(df.tmp1$X1[1:num.sample],df.tmp1$X2[1:num.sample],pch = 20,col = alpha(df.tmp1$Col[1:num.sample],0.5))
}

for (j in 1:size){
  plot.cluster(j)
  # Simulate bivariate normal data, using code from http://blog.revolutionanalytics.com/2016/02/multivariate_data_with_r.html
  mu <- unlist(mean.array[j])  # Mean
  Sigma <-  matrix(unlist(sigma.array[j]),4,4) # Covariance matrix # CRQ: updated this
  bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma)  # from Mass package
  head(bivn)
  # Calculate kernel density estimate
  bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)
  image(bivn.kde) # from base graphics package
  contour(bivn.kde, add = TRUE)
}


## ------------------------------------------------------------------------
# # We note that mvrnorm is faster than rmvnorm
# x.new <- ldply(
#   components,
#   .fun = function(x) {
#     rmvnorm(1,
#             mean = mean.array[[x]],
#             sigma = sigma.array[[x]])
#   }
# )
x.new <- ldply( 
  components,
  .fun = function(x) {
    mvrnorm(1,
            mu = unlist(mean.array[[x]]),
            Sigma = sigma.array[[x]])
  }
)
plot(x.new)

# y.new <- ldply(
#     1:nrow(x.new),
#     .fun = function(x)
#       v(unlist(x.new[x, ]))
#   )
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
# https://stackoverflow.com/questions/38761453/confusion-on-2-dimension-kernel-density-estimation-in-r
# https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/kde2d.html
# https://stackoverflow.com/questions/30814545/r-get-joint-probabilities-from-2d-kernel-density-estimate
# density <- kde2d(x.old[, 1], x.old[, 2], n = 50)
# points <- data.frame(x = x.new[, 1], y = x.new[, 2])
# p.old2 <- interp.surface(density, points)

# ## ----plot8.gmm.kde.plot--------------------------------------------------------------------
# filled.contour(density, color.palette = colorRampPalette(c(
#   'white', 'blue',
#   'yellow', 'red',
#   'darkred'
# )))

## ------------------------------------------------------------------------
p.old <- circle.density.df(as.matrix(x.new[,1:2]))
p.old <- p.old * dmvnorm(x.new[,3:4], c(0,0), 0.01 * diag(2))

p.new <- cluster.size.percentage[1] * dmvnorm(x.new, unlist(mean.array[[1]]), sigma.array[[1]])
for (j in 2:size) {
  p.new <- p.new + cluster.size.percentage[j] * dmvnorm(x.new, unlist(mean.array[[j]]), sigma.array[[j]])
}

# ## ----plot9--------------------------------------------------------------------
# hist(
#   p.new,
#   breaks = seq(min(p.new),
#                max(p.new), l = 35),
#   col = "purple",
#   xlab = "Probabilities from GMM",
#   main = "Histogram of p.new, with n-m Probabilities"
# )

## ------------------------------------------------------------------------
estimator.new <- mean(s.new[, 1] / p.new * p.old)
estimator.new

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------





sample.once.importance <- function(n=2000){
  out <- tryCatch(
    {
      m <- floor (n^(2/3)) # first-stage samples 
      theta <- runif(m, 0, 2 * pi)
      u <- cbind(cos(theta), sin(theta))
      # plot(u)
      
      epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
      x <- u + epsilon
      # plot(x)
      x <- cbind(x,  mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2)))
      
      x.old <- x
      
      v <- function(x) {
        return(rnorm(1, exp(x[1] + x[2]), 1))
      }
      
      ## ------------------------------------------------------------------------
      y <-
        unlist(lapply(1:m, function(x)
          v(x.old[x, ]))) # these are the sampled V's.
      sum(y == 0)/length(y)
      # mean(y) # 1.525384
      # median(y) # 1.377392
      
      ## ------------------------------------------------------------------------
      zeta <- 1.35 # this is chosen to be the median of one trial y. Goal is this represents P(Y>zeta)= 0.5.
      s <- ifelse(y > zeta, 1, 0)
      
      
      # ## ----plot0---------------------------------------------------------------
      # color <- ifelse(y > zeta, "blue", "red")
      # plot(x.old,
      #      col = color,
      #      pch = 20,
      #      main = "First stage samples of x, color by values of y's.")
      # grid()
      # legend(
      #   'bottomright',
      #   c("y>10", "y<10"),
      #   pch = 20,
      #   col = c("blue", "red"),
      #   bty = 'y',
      #   cex = .65
      # )
      # 
      ## ------------------------------------------------------------------------
      size <- floor(sqrt(m))  # number of k-means clusters
      model <- kmeans(x.old, size)
      cluster.size.count <-
        unlist(lapply(1:size, function(x) {
          length(which(model$cluster == x))
        }))
      
      # ## ----plot1--------------------------------------------------------------------
      # barplot(
      #   cluster.size.count,
      #   main = "Size of Clusters",
      #   xlab = "Cluster index",
      #   names.arg = 1:size,
      #   col = "purple"
      # )
      # 
      # ## ----plot2.0---------------------------------------------------------------
      
      # plot(data.frame(x.old),
      #      col = model$cluster,
      #      # xlim = c(-1,1),
      #      # ylim = c(-1,1),
      #      pch = 20)
      # 
      # 
      # ## ----plot2---------------------------------------------------------------
      
      # plot(x.old,
      #      col = model$cluster,
      #      # xlim = c(-1,1),
      #      # ylim = c(-1,1),
      #      pch = 20)
      # # points(model$centers, col = "black", pch = 2)
      # grid()
      # text(
      #   model$centers[, 1],
      #   model$centers[, 2],
      #   labels = as.character(1:size),
      #   col = "black",
      #   pos = c(1),
      #   offset = -0.16
      # )
      # 
      # # NOTE: THESE SHOULD LOOK MORE OR LESS UNIFORMLY DISTRIBUTED
      # plot(x.old[,2:3],
      #      col = model$cluster,
      #      # xlim = c(-1,1),
      #      # ylim = c(-1,1),
      #      pch = 20)
      # # points(model$centers, col = "black", pch = 2)
      # grid()
      # text(
      #   model$centers[, 1],
      #   model$centers[, 2],
      #   labels = as.character(1:size),
      #   col = "black",
      #   pos = c(1),
      #   offset = -0.16
      # )
      # 
      # # NOTE: THESE SHOULD LOOK MORE OR LESS UNIFORMLY DISTRIBUTED
      # plot(x.old[,3:4],
      #      col = model$cluster,
      #      # xlim = c(-1,1),
      #      # ylim = c(-1,1),
      #      pch = 20)
      # # points(model$centers, col = "black", pch = 2)
      # grid()
      # text(
      #   model$centers[, 1],
      #   model$centers[, 2],
      #   labels = as.character(1:size),
      #   col = "black",
      #   pos = c(1),
      #   offset = -0.16
      # )
      # 
      ## ------------------------------------------------------------------------
      cluster.prob <-
        unlist(lapply(1:nrow(model$centers), function(x) {
          mean(s[which(model$cluster == x)])
        }))
      
      # ## ----plot3---------------------------------------------------------------
      # barplot(
      #   cluster.prob,
      #   main = "Proportion of Y>zeta's in the Clusters",
      #   xlab = "Cluster index",
      #   names.arg = 1:size,
      #   col = "purple"
      # )
      
      # ## ----plot4---------------------------------------------------------------
      # pos.prob <- which(cluster.prob > 0)
      # colors.grayscale <-
      #   paste("gray", floor((1 - cluster.prob) * 100), sep = "")
      # plot(
      #   x.old[which(model$cluster == pos.prob[1]), ],
      #   col = colors.grayscale[1],
      #   bg = colors.grayscale[1],
      #   xlim = c(-1.5, 1.5),
      #   ylim = c(-1.5, 1.5),
      #   main = "Clusters and Their Respective Probabilities",
      #   xlab = "x_1",
      #   ylab = "x_2",
      #   pch = 20
      # )
      # grid()
      # for (i in 2:size) {
      #   points(x.old[which(model$cluster == i), ],
      #          col = colors.grayscale[i],
      #          bg = colors.grayscale[i],
      #          pch = 20)
      # }
      # text(
      #   model$centers[, 1],
      #   model$centers[, 2],
      #   labels = paste(
      #     paste("#", as.character(1:size), ":", sep = ""),
      #     paste(floor(cluster.prob * 100), rep("%", size), sep = ""),
      #     sep = " "
      #   ),
      #   col = "red",
      #   pos = c(1),
      #   offset = -0.16
      # )
      # text(model$centers[,1], model$centers[,2], labels=as.character(1:size), col="red", pos=c(1), offset=-0.16)
      # text(model$centers[,1],
      #      model$centers[,2],
      #      labels=paste(floor(cluster.prob*100), rep("%", size), sep = ""),
      #      col="red",
      #      pos=c(4),
      #      offset=1)
      
      
      ## ------------------------------------------------------------------------
      df <- data.frame(model$centers)
      # names(df) <- c("centers.old1", "centers.old2")
      df$prob <- cluster.prob
      # centers.new <- matrix(unlist(lapply(1:size,function(x){colMeans(x.old[which(model$cluster == x),])})),,2)
      # df$centers.new1 <- centers.new[,1]
      # df$centers.new2 <- centers.new[,2]
      # df$centers.new2 <- df$centers.old2
      # df$centers.new1 <- df$centers.old1
      # df
      
      ## ------------------------------------------------------------------------
      sigma.new <- function(i) {
        cov(x.old[which(model$cluster == i), ])
      }
      
      ## ------------------------------------------------------------------------
      cluster.size.percentage <-
        unlist(lapply(1:size, function(x) {
          sum(s[which(model$cluster == x)]) / sum(s)
        }))
      # sum(cluster.size.percentage) # = 1
      
      ## ----plot5--------------------------------------------------------------------
      # barplot(
      #   cluster.size.percentage,
      #   main = "Probability of Clusters",
      #   xlab = "Cluster index",
      #   names.arg = 1:size,
      #   col = "purple"
      # )
      
      ## ------------------------------------------------------------------------
      components <- 
        sample(1:size,
               prob = cluster.size.percentage,
               size = (n - m),
               replace = TRUE)
      
      # ## ----plot6--------------------------------------------------------------------
      # hist(
      #   components,
      #   col  = "purple",
      #   breaks = c(0:(size + 1)),
      #   prob = TRUE,
      #   xaxt = 'n',
      #   xlab = "Component Index"
      # )
      # # lines(density(components))
      # axis(1,
      #      at = seq(0.5, (size + 1.5), by = 1),
      #      labels = c(1:(size + 2)),
      #      las = 1)
      
      ## ------------------------------------------------------------------------
      sigma.array <- lapply(1:size, function(x) unlist(sigma.new(x)))
      mean.array <- lapply(1:size, function(x) df[,1:4][x, ])
      df$sigma <- sigma.array
      # df
      
      # ## -----gmm.unused-------------------------------------------------------------------
      # gmm.pdf <- function(x){ # goal: This function should integrate to 1.
      #   sum(cluster.size.percentage *
      #         unlist(lapply(1:size,
      #                       function(j){dmvnorm(x,
      #                                           mean = mean.array[[j]],
      #                                           sigma = sigma.array[[j]]
      #                       )
      #                       }
      #         )
      #         )
      #   )
      # }
      
      ## ----plot7.visualize.gmm--------------------------------------------------------------------
      # plot.cluster <- function(counter){ # counter indicate cluster indices
      #   
      #   num.sample <- 300 # sample 300 times from each normal.
      #   
      #   this.cluster.points <- rmvnorm(
      #     num.sample, mean = unlist(mean.array[counter]),
      #     sigma = matrix(unlist(sigma.array[counter]),4,4) # CRQ: updated this
      #   ) # generate points from the normal
      #   
      #   df.tmp1 <- data.frame(this.cluster.points) # save the generated points in a data frame
      #   
      #   df.tmp1$prob <- unlist(lapply(1:num.sample, function(x)dmvnorm(this.cluster.points[x,], mean = unlist(mean.array[counter]), sigma = matrix(unlist(sigma.array[counter]),4,4)))) # compute the pdf of the normal at each sampled point
      #   # CRQ: updated this
      #   
      #   df.tmp2 <- data.frame(cbind(x.old[which(model$cluster==counter),],0)) # obtain the sampled points belonging to this cluster from stage-one sampling
      #   
      #   names(df.tmp2) <- names(df.tmp1)
      #   df.tmp1 <- rbind(df.tmp1, df.tmp2) # combine the two temporary data frames
      #   
      #   rbPal <- colorRampPalette(c('red','blue')) # specify color scale.
      #   
      #   df.tmp1$Col <- rbPal(nrow(df.tmp1))[as.numeric(cut(df.tmp1$prob, breaks = nrow(df.tmp1)))] # adds a column of color values based on the y values.
      #   
      #   plot(x.old,
      #        col = model$cluster,
      #        xlim = c(-1.2, 1.2),
      #        ylim = c(-1.2, 1.2),
      #        pch = model$cluster,
      #        main = paste("Sampling from Normal Distribution #", counter))
      #   grid()
      #   points(model$centers, col = 1:size, pch = 2)
      #   text(model$centers[,1], model$centers[,2], labels=as.character(1:size), col="dark orange", pos=c(1), offset=-0.16)
      #   
      #   points(df.tmp1$X1[1:num.sample],df.tmp1$X2[1:num.sample],pch = 20,col = alpha(df.tmp1$Col[1:num.sample],0.5))
      # }
      # 
      # for (j in 1:size){
      #   plot.cluster(j)
      #   # Simulate bivariate normal data, using code from http://blog.revolutionanalytics.com/2016/02/multivariate_data_with_r.html
      #   mu <- unlist(mean.array[j])  # Mean
      #   Sigma <-  matrix(unlist(sigma.array[j]),4,4) # Covariance matrix # CRQ: updated this
      #   bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma)  # from Mass package
      #   head(bivn)
      #   # Calculate kernel density estimate
      #   bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)
      #   image(bivn.kde) # from base graphics package
      #   contour(bivn.kde, add = TRUE)
      # }
      
      
      ## ------------------------------------------------------------------------
      # # We note that mvrnorm is faster than rmvnorm
      # x.new <- ldply(
      #   components,
      #   .fun = function(x) {
      #     rmvnorm(1,
      #             mean = mean.array[[x]],
      #             sigma = sigma.array[[x]])
      #   }
      # )
      x.new <- ldply( 
        components,
        .fun = function(x) {
          mvrnorm(1,
                  mu = unlist(mean.array[[x]]),
                  Sigma = sigma.array[[x]])
        }
      )
      # plot(x.new)
      
      # y.new <- ldply(
      #     1:nrow(x.new),
      #     .fun = function(x)
      #       v(unlist(x.new[x, ]))
      #   )
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
      # https://stackoverflow.com/questions/38761453/confusion-on-2-dimension-kernel-density-estimation-in-r
      # https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/kde2d.html
      # https://stackoverflow.com/questions/30814545/r-get-joint-probabilities-from-2d-kernel-density-estimate
      # density <- kde2d(x.old[, 1], x.old[, 2], n = 50)
      # points <- data.frame(x = x.new[, 1], y = x.new[, 2])
      # p.old2 <- interp.surface(density, points)
      
      # ## ----plot8.gmm.kde.plot--------------------------------------------------------------------
      # filled.contour(density, color.palette = colorRampPalette(c(
      #   'white', 'blue',
      #   'yellow', 'red',
      #   'darkred'
      # )))
      
      ## ------------------------------------------------------------------------
      p.old <- circle.density.df(as.matrix(x.new[,1:2]))
      p.old <- p.old * dmvnorm(x.new[,3:4], c(0,0), 0.01 * diag(2))
      
      p.new <- cluster.size.percentage[1] * dmvnorm(x.new, unlist(mean.array[[1]]), sigma.array[[1]])
      for (j in 2:size) {
        p.new <- p.new + cluster.size.percentage[j] * dmvnorm(x.new, unlist(mean.array[[j]]), sigma.array[[j]])
      }
      
      # ## ----plot9--------------------------------------------------------------------
      # hist(
      #   p.new,
      #   breaks = seq(min(p.new),
      #                max(p.new), l = 35),
      #   col = "purple",
      #   xlab = "Probabilities from GMM",
      #   main = "Histogram of p.new, with n-m Probabilities"
      # )
      
      ## ------------------------------------------------------------------------
      estimator.new <- mean(s.new[, 1] / p.new * p.old)
      # estimator.new
      return(estimator.new)
    },
    error = function(cond){
      return(-999)
    }
  )
  return(out)
}

# results1000 <- ldply(1:1000, .fun = function(x){sample.once.importance(1000)})
# mean(results1000[results1000<1 & results1000>0,1])
# results2000 <- ldply(1:1000, .fun = function(x){sample.once.importance(2000)})
# mean(results2000[results2000<1 & results2000>0,1])

results1000 <- sample.once.importance(1000)
results2000 <- sample.once.importance(2000)
results4000 <- sample.once.importance(4000)
results8000 <- sample.once.importance(8000)

for (tmp.counter in 1:9000){
  results1000 <- c(results1000, sample.once.importance(1000))
  results2000 <- c(results2000, sample.once.importance(2000))
  results4000 <- c(results4000, sample.once.importance(4000))
  results8000 <- c(results8000, sample.once.importance(8000))
  results1000.naive <- c(results1000.naive, sample.once.naive(1000))
  results2000.naive <- c(results2000.naive, sample.once.naive(2000))
  results4000.naive <- c(results4000.naive, sample.once.naive(4000))
  results8000.naive <- c(results8000.naive, sample.once.naive(8000))
}
results1000 <- matrix(unlist(results1000), ncol = 1)[,1]
results2000 <- matrix(unlist(results2000), ncol = 1)[,1]

mean(results4000[results4000<1 & results4000>0])
mean(results8000[results8000<1 & results8000>0])

results16000 <- sample.once.importance(16000)
results16000.naive <- sample.once.naive(16000)
for (tmp.counter in 1:5000){
  results16000 <- c(results16000, sample.once.importance(16000))
  results16000.naive <- c(results16000.naive, sample.once.naive(16000))
}


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

sample.once.naive <- function(n=2000){
  m <- n
  ## ------------------------------------------------------------------------
  zeta <- 1.35
  m <- floor (n^(2/3)) # first-stage samples 
  theta <- runif(m, 0, 2 * pi)
  u <- cbind(cos(theta), sin(theta))
  # plot(u)
  
  epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
  x <- u + epsilon
  # plot(x)
  # x <- cbind(x,  mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2)))
  
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

results1000.naive <- sample.once.naive(1000)
results2000.naive <- sample.once.naive(2000)
results4000.naive <- sample.once.naive(4000)
results8000.naive <- sample.once.naive(8000)

for (tmp.counter in 1:9000){
  results1000.naive <- c(results1000.naive, sample.once.naive(1000))
  results2000.naive <- c(results2000.naive, sample.once.naive(2000))
  results4000.naive <- c(results4000.naive, sample.once.naive(4000))
  results8000.naive <- c(results8000.naive, sample.once.naive(8000))
}


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

df.tmp <- data.frame(results8000[results8000<1 & results8000>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
plot(1:nrow(df.tmp), 
     mse.tmp[,1], 
     type = "l", 
     col = "firebrick",
     lty = "solid",
     main = "MSE with Different Numbers of Samples", 
     xlab = "Number of Simulations",
     ylab = "MSE")

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

df.tmp <- data.frame(results1000[results1000<1 & results1000>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
plot(1:nrow(df.tmp), 
     mse.tmp[,1], 
     type = "l", 
     col = "firebrick",
     lty = "solid",
     main = "MSE with Different Numbers of Samples", 
     xlab = "Number of Simulations",
     ylab = "MSE")
df.tmp <- data.frame(results2000[results2000<1 & results2000>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "aquamarine4",
      lty = "solid")
df.tmp <- data.frame(results4000[results4000<1 & results4000>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "darkorchid4",
      lty = "solid")
df.tmp <- data.frame(results8000[results8000<1 & results8000>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "darkgoldenrod3",
      lty = "solid")
grid()
legend('topright', 
       bg="transparent",
       c("n=1000", "n=2000", "n=4000", "n=8000"), 
       lty = c("solid"),
       col = c("firebrick", "aquamarine4", "darkorchid4", "darkgoldenrod3"), 
       bty = 'y', 
       lwd = 2,
       title = "Importance Sampling",
       title.col = "black",
       cex = 0.5)

df.tmp <- data.frame(results1000.naive[results1000.naive<1 & results1000.naive>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "firebrick1",
      lty = "dotted",
      main = "MSE with n = 2000", 
      xlab = "Number of Simulations",
      ylab = "MSE")
df.tmp <- data.frame(results2000.naive[results2000.naive<1 & results2000.naive>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "mediumaquamarine",
      lty = "dotted")
df.tmp <- data.frame(results4000.naive[results4000.naive<1 & results4000.naive>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "darkorchid1",
      lty = "dotted")
df.tmp <- data.frame(results8000.naive[results8000.naive<1 & results8000.naive>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "darkgoldenrod1",
      lty = "dotted")
legend('topleft', 
       bg="transparent",
       c("n=1000", "n=2000", "n=4000", "n=8000"), 
       lty = c("dotted"),
       col = c("firebrick1", "mediumaquamarine", "darkorchid1", "darkgoldenrod1"), 
       bty = 'y', 
       lwd = 2,
       title = "Naive Sampling",
       title.col = "black",
       cex = 0.5)
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------


df.tmp <- data.frame(results16000.naive[results16000.naive<1 & results16000.naive>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
plot(1:nrow(df.tmp), 
     mse.tmp[,1], 
     type = "l", 
     col = "darkgoldenrod1",
     lty = "dotted")
df.tmp <- data.frame(results16000[results16000<1 & results16000>0])
mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
lines(1:nrow(df.tmp), 
      mse.tmp[,1], 
      type = "l", 
      col = "darkgoldenrod1",
      lty = "dotted")

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------