library(distr)
library(plyr)
library(mclust)
library(MASS)
library(ggplot2)
library(mvtnorm)
library(fields)

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
noise.dim <- 1
noise.var <- 0.00001
n <- 2000 # total samples
m <- floor (n^(2/3)) # first-stage samples, m = 158

theta <- runif(m, 0, 2 * pi)
u <- cbind(cos(theta), sin(theta))
## ----plot u--------------------------------------------------------------------
plot(u, main = "Sampled Points on the Unit Circle", col = "purple", pch = 20)
grid()



epsilon <- mvrnorm(m, mu = c(0, 0), Sigma = 0.01 * diag(2))
x <- u + epsilon
## ----plot x--------------------------------------------------------------------
plot(x, main = "Sampled Noisy x's", col = "purple", pch = 20)
grid()




x <- cbind(x,  mvrnorm(m, 
                       mu = rep(0, noise.dim), 
                       Sigma = noise.var * diag(noise.dim)))



x.old <- x
y <-
  unlist(lapply(1:m, function(x)
    v(x.old[x, ]))) # these are the sampled V's.
sum(y == 0)/length(y)
mean(y)  
median(y) # close to the zeta = 1.35



s <- ifelse(y > zeta, 1, 0)



## ----plot0---------------------------------------------------------------
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




## ----plot0---------------------------------------------------------------
color <- ifelse(y > zeta, "blue", "red")
plot(data.frame(x.old),
     col = color,
     pch = 20,
     main = "First stage samples of x, color by values of y's.")




size <- floor(sqrt(m))  # number of k-means clusters
model <- kmeans(x.old, size)
cluster.size.count <-
  unlist(lapply(1:size, function(x) {
    length(which(model$cluster == x))
  }))

## ----plot the size of clusters--------------------------------------------------------------------
barplot(
  cluster.size.count,
  main = "Size of Clusters",
  xlab = "Cluster index",
  names.arg = 1:size,
  col = "purple"
)





## ----plot2---------------------------------------------------------------

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




## ----plot2.0---------------------------------------------------------------

plot(data.frame(x.old),
     col = model$cluster,
     # xlim = c(-1,1),
     # ylim = c(-1,1),
     pch = 20)





cluster.prob <-
  unlist(lapply(1:nrow(model$centers), function(x) {
    mean(s[which(model$cluster == x)])
  }))

## ----plot---------------------------------------------------------------
barplot(
  cluster.prob,
  main = "Proportion of Y>zeta's in the Clusters",
  xlab = "Cluster index",
  names.arg = 1:size,
  col = "purple"
)






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
for (i in 1:size) {
  points(x.old[which(model$cluster == i), ], col = "purple")
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
  offset = -0.16,
  font = 2
)






df <- data.frame(model$centers)
df$prob <- cluster.prob
## ------------------------------------------------------------------------
sigma.new <- function(i) {
  cov(x.old[which(model$cluster == i), ])
}




cluster.size.percentage <-
  unlist(lapply(1:size, function(x) {
    sum(s[which(model$cluster == x)]) / sum(s)
  })) # total percentage of s == 1 's among all first stage samples
sum(cluster.size.percentage) # = 1




## ----plot5--------------------------------------------------------------------
barplot(
  cluster.size.percentage,
  main = "Probability of Clusters",
  xlab = "Cluster index",
  names.arg = 1:size,
  col = "purple"
)




components <- 
  sample(1:size,
         prob = cluster.size.percentage,
         size = (n - m),
         replace = TRUE)
## ----plot histogram of components--------------------------------------------------------------------
hist(
  components,
  col  = "purple",
  breaks = c(0:(size + 1)),
  prob = TRUE,
  xaxt = 'n',
  xlab = "Component Index"
)
axis(1,
     at = seq(0.5, (size + 1.5), by = 1),
     labels = c(1:(size + 2)),
     las = 1)





sigma.array <- lapply(1:size, function(x) unlist(sigma.new(x)))
mean.array <- lapply(1:size, function(x) df[,1:(2 + noise.dim)][x, ])
df$sigma <- sigma.array

# DIAGNOSIS
for (i in 1:12) {print(sigma.array[[i]] < 9e-05)}
for (i in 1:12) {print(sigma.array[[i]])}

## ----plot to visualize.gmm--------------------------------------------------------------------
plot.cluster <- function(counter){ # counter indicate cluster indices
  
  num.sample <- 300 # sample 300 times from each normal.
  
  this.cluster.points <- rmvnorm(
    num.sample, mean = unlist(mean.array[counter]),
    sigma = matrix(unlist(sigma.array[counter]),(2 + noise.dim), (2 + noise.dim)) # CRQ: updated this
  ) # generate points from the normal
  
  df.tmp1 <- data.frame(this.cluster.points) # save the generated points in a data frame
  
  df.tmp1$prob <- unlist(lapply(1:num.sample, function(x)dmvnorm(this.cluster.points[x,], mean = unlist(mean.array[counter]), sigma = matrix(unlist(sigma.array[counter]),(2 + noise.dim), (2 + noise.dim))))) # compute the pdf of the normal at each sampled point
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
for (j in 1:3){
  plot.cluster(j)
  # Simulate bivariate normal data, using code from http://blog.revolutionanalytics.com/2016/02/multivariate_data_with_r.html
  mu <- unlist(mean.array[j])  # Mean
  Sigma <-  matrix(unlist(sigma.array[j]), (2 + noise.dim), (2 + noise.dim)) # Covariance matrix # CRQ: updated this
  bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma)  # from Mass package
  head(bivn)
  # Calculate kernel density estimate
  bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)
  image(bivn.kde) # from base graphics package
  contour(bivn.kde, add = TRUE)
}




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






density <- kde2d(x.new[, 1], x.new[, 2], n = 50)
points.tmp <- data.frame(x = x.new[, 1], y = x.new[, 2])
p.old2 <- interp.surface(density, points.tmp)



## ----plot gmm kde ----------------------------------------------------------------
filled.contour(density, color.palette = colorRampPalette(c(
  'white', 'blue',
  'yellow', 'red',
  'darkred'
)))





# plot(data.frame(x.new), pch = ".", col = alpha("black", 0.5))
# 
# plot(x.new[,1:2], pch = ".", col = alpha("black", 0.5))


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
estimator.new <- mean(s.new[, 1] / p.new * p.old)
estimator.new
