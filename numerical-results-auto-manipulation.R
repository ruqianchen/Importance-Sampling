library(distr)
library(plyr)
library(mclust)
library(MASS)
library(ggplot2)
library(mvtnorm)
library(fields)


source("alg-imp.r")

# 0.02, 0.002, 0.0002, 0.00002, 0.000002
dat.noise.var <- rbind(dat.noise.var, dat)
test(dat.noise.var)
plot(c(0.02, 0.002, 0.0002, 0.00002, 0.000002), test(dat.noise.var)[1:5])


# 0.017 0.019 0.021 0.023 0.025
dat.noise.var.2 <- dat 
test(dat.noise.var.2)
plot(c(0.017, 0.019, 0.021, 0.023, 0.025), test(dat.noise.var.2)[1:5])


dat.noise.var.0220 <- dat
test(dat.noise.var.0220)
plot(c( 0.0220, 0.0224, 0.0228, 0.0232, 0.0236, 0.0240), test(dat.noise.var.0220)[1:6])


dat.noise.var.00416 <- dat
test(dat.noise.var.00416)
plot(c( 0.00416, 0.00812, 0.01208, 0.01604, 0.02000), test(dat.noise.var.00416)[1:5])


dat.arithmetic.002.02 <- dat
plot(c(0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020), test(dat.arithmetic.002.02)[1:10])

dat.c.arithmetic.1.9 <- dat
test(dat.c.arithmetic.1.9)




test <- function(df){
  t <- ncol(df)
  a <- c()
  b <- c()
  for (i in 1:t){
    temp.col <- df[[i]]
    a <- c(a, mean((temp.col[temp.col > 0 & temp.col < 1]-0.5)^2))
    b <- c(b, mean((temp.col[temp.col > 0 & temp.col < 1])))
    # print(mean((temp.col[temp.col > 0 & temp.col < 1]-0.5)^2))
    # print(mean((temp.col[temp.col > 0 & temp.col < 1])))
  }
  return(c(a,b))
}



dat <- data.frame(matrix(ncol = 5, nrow = 20))
colnames(dat) <- c(0.02, 0.002, 0.0002, 0.00002, 0.000002)
for (i in 1:5){
  noise.tmp <- 2 * 10^(-i-1)
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp)
  for (tmp.counter in 1:19){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp))
  }
  dat[i] <- tmp.array
}
dat.noise.var <- dat






dat <- data.frame(matrix(ncol = 4, nrow = 20))
colnames(dat) <- c(1,2,3,4)
for (i in 1:4){
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = 0.002, c = i)
  for (tmp.counter in 1:19){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = 0.002, c = i))
  }
  dat[[i]] <- tmp.array
}

for (i in 0:5){
  c <- 1.5 + i/5
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = 0.002, c = c)
  for (tmp.counter in 1:19){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = 0.002, c = c))
  }
  dat[,as.character(c)] <- tmp.array
}

dat.backup2 <- dat





dat <- data.frame(matrix(ncol = 4, nrow = 20))
colnames(dat) <- c(1,2,3,4)
for (i in 1:4){
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = 0.002, c = i)
  for (tmp.counter in 1:19){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = 0.002, c = i))
  }
  dat[[i]] <- tmp.array
}





dat <- data.frame(matrix(ncol = 5, nrow = 2))
for (i in 1:5){
  noise.tmp <- 0.015 + 2 * 10^(-3) * i
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp)
  for (tmp.counter in 1:2){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp))
  }
  dat[i] <- tmp.array
}



total <- 2000
dat <- data.frame(matrix(ncol = 9, nrow = total))
for (i in 1:9){
  noise.tmp <- 0.002
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 0.1*i)
  for (tmp.counter in 1:(total-1)){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 0.1*i))
  }
  dat[i] <- tmp.array
}


total <- 1500
dat <- data.frame(matrix(ncol = 10, nrow = total))
for (i in 1:9){
  noise.tmp <- 0.002
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 0.1*i)
  for (tmp.counter in 1:(total-1)){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp))
  }
  dat[i] <- tmp.array
}

total <- 150
dat <- data.frame(matrix(ncol = 2, nrow = total))
for (i in 1:2){
  noise.tmp <- 0.002
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 10^(-i))
  for (tmp.counter in 1:(total-1)){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp))
  }
  dat[i] <- tmp.array
}
dat.backup <- dat

total <- 150
dat <- data.frame(matrix(ncol = 2, nrow = total))
for (i in 1:2){
  noise.tmp <- 0.002
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 20^(-i))
  for (tmp.counter in 1:(total-1)){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp))
  }
  dat[i] <- tmp.array
}


total <- 30
dat <- data.frame(matrix(ncol = 1, nrow = total))
for (i in 1:1){
  noise.tmp <- 0.002
  tmp.array <- sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 5e-04)
  for (tmp.counter in 1:(total-1)){
    tmp.array <- c(tmp.array, sample.once.importance(noise.dim = 1, noise.var = noise.tmp, c = 5e-04))
  }
  dat[i] <- tmp.array
}



tmp1 <- dat[["2e-04"]]
plot(tmp1[tmp1>0 & tmp1 < 1])


results1000 <- sample.once.importance(1000)
results2000 <- sample.once.importance(2000)
results4000 <- sample.once.importance(4000)
results8000 <- sample.once.importance(8000)
results16000 <- sample.once.importance(16000)
results1000.naive <- sample.once.naive(1000)
results2000.naive <- sample.once.naive(2000)
results4000.naive <- sample.once.naive(4000)
results8000.naive <- sample.once.naive(8000)
results16000.naive <- sample.once.naive(16000)

for (tmp.counter in 1:2000){
  results1000 <- c(results1000, sample.once.importance(1000))
  results2000 <- c(results2000, sample.once.importance(2000))
  results4000 <- c(results4000, sample.once.importance(4000))
  # results8000 <- c(results8000, sample.once.importance(8000))
  # results16000 <- c(results16000, sample.once.importance(16000))
  results1000.naive <- c(results1000.naive, sample.once.naive(1000))
  results2000.naive <- c(results2000.naive, sample.once.naive(2000))
  results4000.naive <- c(results4000.naive, sample.once.naive(4000))
  # results8000.naive <- c(results8000.naive, sample.once.naive(8000))
  # results16000.naive <- c(results16000.naive, sample.once.naive(16000))
}


# 
# results2000.1 <- sample.once.importance()
# for (tmp.counter in 1:20){
#   results2000.1 <- c(results2000.1, sample.once.importance(2000))
# }
# df.tmp <- data.frame(results2000.1[results2000.1<1 & results2000.1>0])
# mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
# lines(1:nrow(df.tmp), 
#      mse.tmp[,1], 
#      type = "l", 
#      col = "firebrick",
#      lty = "solid",
#      main = "MSE with Different Numbers of Samples", 
#      xlab = "Number of Simulations",
#      ylab = "MSE")
# 
# results2000.naive <- sample.once.naive()
# for (tmp.counter in 1:1810){
#   results2000.naive <- c(results2000.naive, sample.once.naive())
# }
# df.tmp <- data.frame(results2000.naive[results2000.naive<1 & results2000.naive>0])
# mse.tmp <- ldply(1:nrow(df.tmp), .fun = function(x){mean((df.tmp[1:x,1]-0.5)^2)})
# plot(1:nrow(df.tmp), 
#      mse.tmp[,1], 
#      type = "l", 
#      col = "firebrick",
#      lty = "dashed")
