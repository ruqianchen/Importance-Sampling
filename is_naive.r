## ---- results='hide', message=FALSE, warning=FALSE-----------------------
library(distr)
library(plyr)
library(mclust)
library(MASS)
library(ggplot2)
library(mvtnorm)

sample.once.naive <- function(n=2000){
  m <- n
  ## ------------------------------------------------------------------------
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
  
  ## ------------------------------------------------------------------------
  y <- unlist(lapply(1:m, function(x) v(x.old[x,]))) # these are the sampled V's. 
  sum(y == 0)/length(y)  # This should go to 0.75 as m goes to infinity.
  
  ## ------------------------------------------------------------------------
  s <- ifelse(y>zeta, 1, 0)
  return(mean(s))
}


sample.once.importance <- function(n=2000){
  m <- floor (n^(2/3)) # first-stage samples 
  
  ## ------------------------------------------------------------------------
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
  
  ## ------------------------------------------------------------------------
  y <- unlist(lapply(1:m, function(x) v(x.old[x,]))) # these are the sampled V's. 
  sum(y == 0)/length(y)  # This should go to 0.75 as m goes to infinity.
  
  ## ------------------------------------------------------------------------
  s <- ifelse(y>zeta, 1, 0)
  
  # ## ----plot0---------------------------------------------------------------
  # color <- ifelse(y>=zeta,"blue","red")
  # plot(x.old, col=color, pch = 20, main = "First stage samples of x, color by values of y's.")
  # grid()
  # legend('bottomright', c("y>10", "y<10"), pch = 20, col=c("blue", "red"), bty='y', cex=.65)

  ## ------------------------------------------------------------------------
  size <- floor(sqrt(m))  # number of k-means clusters
  model <- kmeans(x.old, size)
  cluster.size.count <- unlist(lapply(1:size,function(x){length(which(model$cluster == x))}))
  
  # ## ----plot1--------------------------------------------------------------------
  # barplot(cluster.size.count, main="Size of Clusters",
  #   	xlab="Cluster index",names.arg=1:size,col = "purple")

  # ## ----plot2---------------------------------------------------------------
  # plot(x.old, col = model$cluster, pch = 20, xlim = c(-1,1), ylim = c(-1,1))
  # # points(model$centers, col = "black", pch = 2)
  # grid()
  # text(model$centers[,1], model$centers[,2], labels = as.character(1:size), col= "black", pos=c(1), offset=-0.16)

  ## ------------------------------------------------------------------------
  s <- ifelse(y>zeta, 1, 0)
  cluster.prob <- unlist(lapply(1:nrow(model$centers),function(x){mean(s[which(model$cluster == x)])}))
  
  # ## ----plot3---------------------------------------------------------------
  # barplot(cluster.prob, main="Proportion of Y>10's in the Clusters",
  #   	xlab="Cluster index",names.arg=1:size, col = "purple")

  # ## ----plot4---------------------------------------------------------------
  # pos.prob <- which(cluster.prob>0)
  # plot(x.old[which(model$cluster == pos.prob[1]),], col=1, xlim = c(-1,1),ylim = c(-1,1), main = "Clusters with Positive Probability", pch=20)
  # grid()
  # for (i in 2:length(pos.prob)){
  #   j = which(cluster.prob>0)[i]
  #   points(x.old[which(model$cluster == j),], col=i, pch = 20)
  # }
  # text(model$centers[,1], model$centers[,2], labels=as.character(1:size), col="black", pos=c(1), offset=-0.16)

  ## ------------------------------------------------------------------------
  df <- data.frame(model$centers)
  names(df) <- c("centers.old1", "centers.old2")
  df$prob <- cluster.prob
  # centers.new <- matrix(unlist(lapply(1:size,function(x){colMeans(x.old[which(model$cluster == x),])})),,2)
  # df$centers.new1 <- centers.new[,1]
  # df$centers.new2 <- centers.new[,2]
  df$centers.new2 <- df$centers.old2
  df$centers.new1 <- df$centers.old1
  # df
  
  ## ------------------------------------------------------------------------
  sigma.new <- function(i){
    cov(x.old[which(model$cluster==i),])
  }
  
  ## ------------------------------------------------------------------------
  cluster.size.percentage <- unlist(lapply(1:size, function(x){sum(s[which(model$cluster == x)])/sum(s)}))
  sum(cluster.size.percentage) # = 1
  
  # ## ----plot5--------------------------------------------------------------------
  # barplot(cluster.size.percentage, main="Probability of Clusters", xlab="Cluster index", names.arg=1:size, col = "purple")

  ## ------------------------------------------------------------------------
  components <- sample(1:size, prob=cluster.size.percentage, size=(n-m), replace=TRUE)
  
  # ## ----plot6--------------------------------------------------------------------
  # hist(components, col  = "purple", breaks = c(0:(size+1)), prob = TRUE, xaxt = 'n', xlab = "Component Index")
  # # lines(density(components))
  # axis(1, at = seq(0.5, (size+1.5), by = 1), labels = c(1:(size+2)), las = 1)

  ## ------------------------------------------------------------------------
  sigma.array <- lapply(1:size, function(x) unlist(sigma.new(x)))
  mean.array <- lapply(1:size, function(x) cbind(df$centers.new1,df$centers.new2)[x,])
  df$sigma <- sigma.array
  # df 
  
  ## -----gmm.unused-------------------------------------------------------------------
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
  
  # ## ----plot7.visualize.gmm--------------------------------------------------------------------
  # plot.cluster <- function(counter){ # counter indicate cluster indices
  # 
  #   num.sample <- 300 # sample 300 times from each normal.
  # 
  #   this.cluster.points <- rmvnorm(
  #     num.sample, mean = unlist(mean.array[counter]),
  #     sigma = matrix(unlist(sigma.array[counter]),2,2)
  #     ) # generate points from the normal
  # 
  #   df.tmp1 <- data.frame(this.cluster.points) # save the generated points in a data frame
  # 
  #   df.tmp1$prob <- unlist(lapply(1:num.sample, function(x)dmvnorm(this.cluster.points[x,], mean = unlist(mean.array[counter]), sigma = matrix(unlist(sigma.array[counter]),2,2)))) # compute the pdf of the normal at each sampled point
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
  #   plot(x.old, col = model$cluster, xlim = c(-1,1),ylim = c(-1,1), pch = model$cluster,main =paste("Sampling from Normal Distribution #", counter))
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
  #   Sigma <-  matrix(unlist(sigma.array[j]),2,2) # Covariance matrix
  #   bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )  # from Mass package
  #   head(bivn)
  #   # Calculate kernel density estimate
  #   bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)
  #   image(bivn.kde) # from base graphics package
  #   contour(bivn.kde, add = TRUE)
  # }

  
  ## ------------------------------------------------------------------------
  x.new <- ldply(
    components, 
    .fun = function(x){rmvnorm(1, mean = unlist(mean.array[x]), sigma = matrix(unlist(sigma.array[x]),2,2) )}
  )
  
  p.old <- rep(0, nrow(x.new))
  p.old[x.new[,1] > -1 & x.new[,1] < 1 & x.new[,2] > -1 & x.new[,2] < 1] <- 0.25
  
  # ## ------------------------------------------------------------------------
  # plot(x.new, col = ifelse(p.old==0.25," black", "red"))
  
  p.new <- cluster.size.percentage[1] * dmvnorm(x.new, mean.array[[1]], sigma.array[[1]])
  for (j in 2:size){
    p.new <- p.new +  cluster.size.percentage[j] * dmvnorm(x.new, mean.array[[j]], sigma.array[[j]])
  }
 
  # ## ----plot8--------------------------------------------------------------------
  # hist(p.new,
  #      breaks=seq(min(p.new),
  #                 max(p.new), l=35),
  #      col = "purple",
  #      xlab = "Probabilities from GMM",
  #      main = "Histogram of p.new, with n-m Probabilities")

 
  ## ------------------------------------------------------------------------
  y.new <- ldply(1:nrow(x.new), .fun=function(x) v(unlist(x.new[x,])))
  s.new <- ldply(1:nrow(x.new), .fun=function(x) {ifelse(y.new[x,] > zeta,1,0)})
  estimator.new <- mean(s.new[,1]/p.new*p.old)
  # estimator.new
  
  # ## ---plot9---------------------------------------------------------------------
  # plot(x.old,
  #      col = color,
  #      pch = 20,
  #      main = "First stage samples of x, color by values of y's.")
  # grid()
  # legend('bottomright', c("y>10", "y<10"), pch = 20, col=c("blue", "red"), bty='y', cex=.65)

  # ## ----plot10--------------------------------------------------------------------
  # col.new <- ifelse(s.new == 1, "blue", "red")
  # plot(x.new,
  #      col=col.new,
  #      pch = 20,
  #      cex = 0.25,
  #      xlim = c(-1,1),
  #      ylim = c(-1,1),
  #      xlab = "X1",
  #      ylab = "X2",
  #      main = "Second stage samples of x, color by values of y's.")
  # grid()
  # legend('bottomright', c("y>10", "y<10"), pch = 20, col = c("blue", "red"), bty = 'y', cex = .65)

  ## ------------------------------------------------------------------------
  return(estimator.new)
}