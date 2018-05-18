#### IS
### high-D case
library(FNN)
###
g_function = function(x){
  x>1.5
}
###
mu0 = 0.4092262

##### large simulation (100 times)
N_rep = 100
V_err = 0.5
d_err = 5
err_dim = 0.01
n = 400
m = 200
k0 = 20


###
mu_naive = rep(NA, N_rep)
mu_is = rep(NA,N_rep)


###
for(i_rep in 1:N_rep){
  U = runif(m, max=2*pi)
  X0 = cos(U)
  Y0 = sin(U)
  
  D = cbind(X0,Y0, matrix(rep(0, m*d_err), ncol=d_err))+
    cbind(matrix(rnorm(2*m, sd=0.1),ncol=2), 
          matrix(rnorm(d_err*m,sd=err_dim), ncol=d_err))
  V = rnorm(m, mean=exp(D[,1]+D[,2]), sd=V_err)
  
  ### new datasets
  U = runif(5*n, max=2*pi)
  X1 = cos(U)
  Y1 = sin(U)
  D1 = cbind(X1,Y1, matrix(rep(0, length(X1)*d_err), ncol=d_err))+
    cbind(matrix(rnorm(2*length(X1), sd=0.1),ncol=2), 
          matrix(rnorm(d_err*length(X1),sd=err_dim), ncol=d_err))
  
  
  D1_fit = knn.reg(train = D, test=D1, y=g_function(V)^2, k = k0)
  # D1_fit_sum = mean(sqrt(
  #   knn.reg(train = D, test=D1, y=g_function(V)^2, k = k0)$pred))
  D1_fit_sum = mean(sqrt(D1_fit$pred))
  
  # max(D1_fit$pred)
  # min(D1_fit$pred)
  
  p_max = max(sqrt(D1_fit$pred))
  
  idx_pass = runif(length(U))< (sqrt(D1_fit$pred)/p_max)
  
  D_pass = D1[idx_pass,]
  p_pass = (sqrt(D1_fit$pred)/p_max)[idx_pass]
  
  D_use = D_pass[1:n,]
  p_use = p_pass[1:n]
  V_use = rnorm(n, mean=exp(D_use[,1]+D_use[,2]), sd=V_err)
  
  mu_is[i_rep] = mean(g_function(V_use)/p_use)*D1_fit_sum
  
  ### naive
  U = runif(n, max=2*pi)
  X2 = cos(U)
  Y2 = sin(U)
  
  D2 = cbind(X2,Y2, matrix(rep(0, n*d_err), ncol=d_err))+
    cbind(matrix(rnorm(2*n, sd=0.1),ncol=2), 
          matrix(rnorm(d_err*n,sd=err_dim), ncol=d_err))
  
  V2 = rnorm(n, mean=exp(D2[,1]+D2[,2]), sd=V_err)
  
  mu = mean(g_function(V2))
  mu_naive[i_rep] = mu
  
}


sd(mu_is)
sd(mu_naive)

sqrt(mean((mu_is - mu0)^2))
sqrt(mean((mu_naive-mu0)^2))


### define function

sample.new <- function(V_err = 0.5, d_err = 5, err_dim = 0.01, N_rep = 100, n = 400, m = 200, r = 1){
  mu_naive = rep(NA, N_rep)
  mu_is = rep(NA,N_rep)
  for(i_rep in 1:N_rep){
    U = runif(m, max=2*pi)
    X0 = r* cos(U)
    Y0 = r* sin(U)
    
    D = cbind(X0,Y0, matrix(rep(0, m*d_err), ncol=d_err))+
      cbind(matrix(rnorm(2*m, sd=0.1),ncol=2), 
            matrix(rnorm(d_err*m,sd=err_dim), ncol=d_err))
    V = rnorm(m, mean=exp(D[,1]+D[,2]), sd=V_err)
    
    ### new datasets
    U = runif(5*n, max=2*pi)
    X1 = r* cos(U)
    Y1 = r* sin(U)
    D1 = cbind(X1,Y1, matrix(rep(0, length(X1)*d_err), ncol=d_err))+
      cbind(matrix(rnorm(2*length(X1), sd=0.1),ncol=2), 
            matrix(rnorm(d_err*length(X1),sd=err_dim), ncol=d_err))
    D1_fit = knn.reg(train = D, test=D1, y=g_function(V)^2, k = k0)
    D1_fit_sum = mean(sqrt(
      knn.reg(train = D, test=D1, y=g_function(V)^2, k = k0)$pred))
    
    # max(D1_fit$pred)
    # min(D1_fit$pred)
    
    p_max = max(sqrt(D1_fit$pred))
    
    idx_pass = runif(length(U))< (sqrt(D1_fit$pred)/p_max)
    
    D_pass = D1[idx_pass,]
    p_pass = (sqrt(D1_fit$pred)/p_max)[idx_pass]
    
    D_use = D_pass[1:n,]
    p_use = p_pass[1:n]
    V_use = rnorm(n, mean=exp(D_use[,1]+D_use[,2]), sd=V_err)
    
    mu_is[i_rep] = mean(g_function(V_use)/p_use)*D1_fit_sum
    
    ### naive
    U = runif(n, max=2*pi)
    X2 = r* cos(U)
    Y2 = r* sin(U)
    
    D2 = cbind(X2,Y2, matrix(rep(0, n*d_err), ncol=d_err))+
      cbind(matrix(rnorm(2*n, sd=0.1),ncol=2), 
            matrix(rnorm(d_err*n,sd=err_dim), ncol=d_err))
    
    V2 = rnorm(n, mean=exp(D2[,1]+D2[,2]), sd=V_err)
    
    mu = mean(g_function(V2))
    mu_naive[i_rep] = mu
  }
  return(cbind(mu_is, mu_naive))
}

performance <- function(input){
  return(c(sqrt(mean((input[,1] - mu0)^2)), sqrt(mean((input[,2]-mu0)^2))))
}

### test
sample.0.5..5..0.01..1000..1000..800 <- sample.new(V_err = 0.5, d_err = 5, n = 1000, m = 800)


### varying error dimension
sample.0.5..5..0.1..1000..400..200 <- sample.new(V_err = 0.5, 
                                                  d_err = 5, 
                                                  err_dim = 0.1, 
                                                  N_rep = 1000, 
                                                  n = 400, 
                                                  m = 200
                                                  )
sample.0.5..5..0.2..1000..400..200 <- sample.new(err_dim = 0.2)
sample.0.5..5..0.3..1000..400..200 <- sample.new(err_dim = 0.3)
sample.0.5..5..0.4..1000..400..200 <- sample.new(err_dim = 0.4)
sample.0.5..5..0.5..1000..400..200 <- sample.new(err_dim = 0.5)
sample.0.5..5..0.6..1000..400..200 <- sample.new(err_dim = 0.6)
sample.0.5..5..0.7..1000..400..200 <- sample.new(err_dim = 0.7)
sample.0.5..5..0.8..1000..400..200 <- sample.new(err_dim = 0.8)
sample.0.5..5..0.9..1000..400..200 <- sample.new(err_dim = 0.9)
sample.0.5..5..1.0..1000..400..200 <- sample.new(err_dim = 1.0)
result <- rbind(
  performance(sample.0.5..5..0.1..1000..400..200),
  performance(sample.0.5..5..0.2..1000..400..200),
  performance(sample.0.5..5..0.3..1000..400..200),
  performance(sample.0.5..5..0.4..1000..400..200),
  performance(sample.0.5..5..0.5..1000..400..200),
  performance(sample.0.5..5..0.6..1000..400..200),
  performance(sample.0.5..5..0.7..1000..400..200),
  performance(sample.0.5..5..0.8..1000..400..200),
  performance(sample.0.5..5..0.9..1000..400..200),
  performance(sample.0.5..5..1.0..1000..400..200)
) 
plot(0.1*c(1:10), 
     result[,1], 
     ylim = c(0,0.03),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5 ||  err_dim varies || n = 400 || m = 200",
     xlab = "err_dim = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0",
     ylab = "MSE")
points(0.1*c(1:10), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)
ticks <- 0.1*c(1:10)
axis(side = 1, at = ticks)
abline(h=0.01*(0:12), v= 0.1*c(1:10), col="gray", lty=3)




### Varying  n while holding m = n/2
sample.0.5..5..0.01..1000..400..200 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 200)
sample.0.5..5..0.01..1000..600..300 <- sample.new(V_err = 0.5, d_err = 5, n = 600, m = 300)
sample.0.5..5..0.01..1000..800..400 <- sample.new(V_err = 0.5, d_err = 5, n = 800, m = 400)
sample.0.5..5..0.01..1000..1000..500 <- sample.new(V_err = 0.5, d_err = 5, n = 1000, m = 500)
sample.0.5..5..0.01..1000..1200..600 <- sample.new(V_err = 0.5, d_err = 5, n = 1200, m = 600)
sample.0.5..5..0.01..1000..1400..700 <- sample.new(V_err = 0.5, d_err = 5, n = 1400, m = 700)
sample.0.5..5..0.01..1000..2000..1000 <- sample.new(V_err = 0.5, d_err = 5, n = 2000, m = 1000)
sample.0.5..5..0.01..1000..3000..1500 <- sample.new(V_err = 0.5, d_err = 5, n = 3000, m = 1500)

result <- rbind(
  performance(sample.0.5..5..0.01..1000..400..200),
  performance(sample.0.5..5..0.01..1000..600..300),
  performance(sample.0.5..5..0.01..1000..800..400),
  performance(sample.0.5..5..0.01..1000..1000..500),
  performance(sample.0.5..5..0.01..1000..1200..600),
  performance(sample.0.5..5..0.01..1000..1400..700),
  performance(sample.0.5..5..0.01..1000..2000..1000),
  performance(sample.0.5..5..0.01..1000..3000..1500)
) 
plot(c(200*c(2:7),2000, 3000), 
     result[,1], 
     ylim = c(0,0.03),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5 ||  err_dim = 0.01 || n varies || m = n/2",
     xlab = "n = 400, 600, 800, 1000, 1200, 1400, 2000, 3000",
     ylab = "MSE")
points(c(200*c(2:7),2000, 3000), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### Varying m while holding n constant
sample.0.5..5..0.01..1000..400..25 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 25)
sample.0.5..5..0.01..1000..400..31 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 31)
sample.0.5..5..0.01..1000..400..37 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 37)
sample.0.5..5..0.01..1000..400..43 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 43)
sample.0.5..5..0.01..1000..400..50 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 50)
sample.0.5..5..0.01..1000..400..56 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 56)
sample.0.5..5..0.01..1000..400..62 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 62)
sample.0.5..5..0.01..1000..400..68 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 68)
sample.0.5..5..0.01..1000..400..75 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 75)
sample.0.5..5..0.01..1000..400..100 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 100)
sample.0.5..5..0.01..1000..400..125 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 125)
sample.0.5..5..0.01..1000..400..150 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 150)
sample.0.5..5..0.01..1000..400..175 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 175)
sample.0.5..5..0.01..1000..400..200 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 200)
sample.0.5..5..0.01..1000..400..225 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 225)
sample.0.5..5..0.01..1000..400..250 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 250)
sample.0.5..5..0.01..1000..400..275 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 275)
sample.0.5..5..0.01..1000..400..300 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 300)
sample.0.5..5..0.01..1000..400..325 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 325)
sample.0.5..5..0.01..1000..400..350 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 350)
sample.0.5..5..0.01..1000..400..375 <- sample.new(V_err = 0.5, d_err = 5, n = 400, m = 375)

result <- rbind(
  performance(sample.0.5..5..0.01..1000..400..25),
  performance(sample.0.5..5..0.01..1000..400..31),
  performance(sample.0.5..5..0.01..1000..400..37),
  performance(sample.0.5..5..0.01..1000..400..43),
  performance(sample.0.5..5..0.01..1000..400..50),
  performance(sample.0.5..5..0.01..1000..400..56),
  performance(sample.0.5..5..0.01..1000..400..62),
  performance(sample.0.5..5..0.01..1000..400..68),
  performance(sample.0.5..5..0.01..1000..400..75),
  performance(sample.0.5..5..0.01..1000..400..100),
  performance(sample.0.5..5..0.01..1000..400..125),
  performance(sample.0.5..5..0.01..1000..400..150),
  performance(sample.0.5..5..0.01..1000..400..175),
  performance(sample.0.5..5..0.01..1000..400..200),
  performance(sample.0.5..5..0.01..1000..400..225),
  performance(sample.0.5..5..0.01..1000..400..250),
  performance(sample.0.5..5..0.01..1000..400..275),
  performance(sample.0.5..5..0.01..1000..400..300),
  performance(sample.0.5..5..0.01..1000..400..325),
  performance(sample.0.5..5..0.01..1000..400..350),
  performance(sample.0.5..5..0.01..1000..400..375)
) 
plot(c(25,31,37,43,50,56,62,68,25*c(3:15)), 
     result[,1], 
     ylim = c(0,0.12),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5 ||  err_dim = 0.01 || n = 400 || m varies",
     xlab = "m = 25, 31, 37, 43, 50, 56, 62, 68, 75, 100 ..., 350, 375",
     ylab = "MSE")
points(c(25,31,37,43,50,56,62,68,25*c(3:15)), result[,2], col = "red", pch = 20)
legend('topright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)
ticks <- c(25*c(1:15))
axis(side = 1, at = ticks)
abline(h=0.01*(0:12), v=c(25*c(1:15)), col="gray", lty=3)
points(0.05*c(1:30), result[,2], col = "red", pch = 20)



### varying d_err
sample.0.5..5..0.01..1000 <- sample.new()
sample.0.5..10..0.01..1000 <- sample.new(d_err = 10)
sample.0.5..15..0.01..1000 <- sample.new(d_err = 15)
sample.0.5..20..0.01..1000 <- sample.new(d_err = 20)
sample.0.5..25..0.01..1000 <- sample.new(d_err = 25)
sample.0.5..30..0.01..1000 <- sample.new(d_err = 30)
sample.0.5..35..0.01..1000 <- sample.new(d_err = 35)
sample.0.5..40..0.01..1000 <- sample.new(d_err = 40)
sample.0.5..45..0.01..1000 <- sample.new(d_err = 45)
sample.0.5..50..0.01..1000 <- sample.new(d_err = 50)
sample.0.5..55..0.01..1000 <- sample.new(d_err = 55)

result <- rbind(
performance(sample.0.5..5..0.01..1000),
performance(sample.0.5..10..0.01..1000),
performance(sample.0.5..15..0.01..1000),
performance(sample.0.5..20..0.01..1000),
performance(sample.0.5..25..0.01..1000),
performance(sample.0.5..30..0.01..1000),
performance(sample.0.5..35..0.01..1000),
performance(sample.0.5..40..0.01..1000),
performance(sample.0.5..45..0.01..1000),
performance(sample.0.5..50..0.01..1000),
performance(sample.0.5..55..0.01..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.03),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err varies  ||  err_dim = 0.01",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)





### varying d_err to larger values

sample.0.5..50..0.01..1000 <- sample.new(d_err = 50)
sample.0.5..100..0.01..1000 <- sample.new(d_err = 100)
sample.0.5..150..0.01..1000 <- sample.new(d_err = 150)
sample.0.5..200..0.01..1000 <- sample.new(d_err = 200)
sample.0.5..250..0.01..1000 <- sample.new(d_err = 250)
sample.0.5..300..0.01..1000 <- sample.new(d_err = 300)
sample.0.5..350..0.01..1000 <- sample.new(d_err = 350)
sample.0.5..400..0.01..1000 <- sample.new(d_err = 400)
sample.0.5..450..0.01..1000 <- sample.new(d_err = 450)
sample.0.5..500..0.01..1000 <- sample.new(d_err = 500)
sample.0.5..550..0.01..1000 <- sample.new(d_err = 550)
sample.0.5..600..0.01..1000 <- sample.new(d_err = 600)
sample.0.5..650..0.01..1000 <- sample.new(d_err = 650)

result <- rbind(
  performance(sample.0.5..50..0.01..1000),
  performance(sample.0.5..100..0.01..1000),
  performance(sample.0.5..150..0.01..1000),
  performance(sample.0.5..200..0.01..1000),
  performance(sample.0.5..250..0.01..1000),
  performance(sample.0.5..300..0.01..1000),
  performance(sample.0.5..350..0.01..1000),
  performance(sample.0.5..400..0.01..1000),
  performance(sample.0.5..450..0.01..1000),
  performance(sample.0.5..500..0.01..1000),
  performance(sample.0.5..550..0.01..1000),
  performance(sample.0.5..600..0.01..1000),
  performance(sample.0.5..650..0.01..1000)
)

plot(50*c(1:13), 
     result[,1], 
     ylim = c(0,0.03),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err varies  ||  err_dim = 0.01",
     xlab = "d_err = 50, 100, 150, ..., 600, 650",
     ylab = "MSE")
points(50*c(1:13), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### Varying d_err while err_dim is high (err_dim = 0.5)

sample.0.5..5..0.5..1000 <- sample.new(err_dim = 0.5)
sample.0.5..10..0.5..1000 <- sample.new(d_err = 10, err_dim = 0.5)
sample.0.5..15..0.5..1000 <- sample.new(d_err = 15, err_dim = 0.5)
sample.0.5..20..0.5..1000 <- sample.new(d_err = 20, err_dim = 0.5)
sample.0.5..25..0.5..1000 <- sample.new(d_err = 25, err_dim = 0.5)
sample.0.5..30..0.5..1000 <- sample.new(d_err = 30, err_dim = 0.5)
sample.0.5..35..0.5..1000 <- sample.new(d_err = 35, err_dim = 0.5)
sample.0.5..40..0.5..1000 <- sample.new(d_err = 40, err_dim = 0.5)
sample.0.5..45..0.5..1000 <- sample.new(d_err = 45, err_dim = 0.5)
sample.0.5..50..0.5..1000 <- sample.new(d_err = 50, err_dim = 0.5)
sample.0.5..55..0.5..1000 <- sample.new(d_err = 55, err_dim = 0.5)
result <- rbind(
  performance(sample.0.5..5..0.5..1000),
  performance(sample.0.5..10..0.5..1000),
  performance(sample.0.5..15..0.5..1000),
  performance(sample.0.5..20..0.5..1000),
  performance(sample.0.5..25..0.5..1000),
  performance(sample.0.5..30..0.5..1000),
  performance(sample.0.5..35..0.5..1000),
  performance(sample.0.5..40..0.5..1000),
  performance(sample.0.5..45..0.5..1000),
  performance(sample.0.5..50..0.5..1000),
  performance(sample.0.5..55..0.5..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 0.5",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### Varying d_err while err_dim is even higher (err_dim = 1)

sample.0.5..5..1..1000 <- sample.new(err_dim = 1)
sample.0.5..10..1..1000 <- sample.new(d_err = 10, err_dim = 1)
sample.0.5..15..1..1000 <- sample.new(d_err = 15, err_dim = 1)
sample.0.5..20..1..1000 <- sample.new(d_err = 20, err_dim = 1)
sample.0.5..25..1..1000 <- sample.new(d_err = 25, err_dim = 1 )
sample.0.5..30..1..1000 <- sample.new(d_err = 30, err_dim = 1)
sample.0.5..35..1..1000 <- sample.new(d_err = 35, err_dim = 1)
sample.0.5..40..1..1000 <- sample.new(d_err = 40, err_dim = 1)
sample.0.5..45..1..1000 <- sample.new(d_err = 45, err_dim = 1)
sample.0.5..50..1..1000 <- sample.new(d_err = 50, err_dim = 1)
sample.0.5..55..1..1000 <- sample.new(d_err = 55, err_dim = 1)
result <- rbind(
  performance(sample.0.5..5..1..1000),
  performance(sample.0.5..10..1..1000),
  performance(sample.0.5..15..1..1000),
  performance(sample.0.5..20..1..1000),
  performance(sample.0.5..25..1..1000),
  performance(sample.0.5..30..1..1000),
  performance(sample.0.5..35..1..1000),
  performance(sample.0.5..40..1..1000),
  performance(sample.0.5..45..1..1000),
  performance(sample.0.5..50..1..1000),
  performance(sample.0.5..55..1..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 1",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)



### Varying d_err while err_dim is even higher (err_dim = 2)

sample.0.5..5..2..1000 <- sample.new(err_dim = 2)
sample.0.5..10..2..1000 <- sample.new(d_err = 10, err_dim = 2)
sample.0.5..15..2..1000 <- sample.new(d_err = 15, err_dim = 2)
sample.0.5..20..2..1000 <- sample.new(d_err = 20, err_dim = 2)
sample.0.5..25..2..1000 <- sample.new(d_err = 25, err_dim = 2)
sample.0.5..30..2..1000 <- sample.new(d_err = 30, err_dim = 2)
sample.0.5..35..2..1000 <- sample.new(d_err = 35, err_dim = 2)
sample.0.5..40..2..1000 <- sample.new(d_err = 40, err_dim = 2)
sample.0.5..45..2..1000 <- sample.new(d_err = 45, err_dim = 2)
sample.0.5..50..2..1000 <- sample.new(d_err = 50, err_dim = 2)
sample.0.5..55..2..1000 <- sample.new(d_err = 55, err_dim = 2)
result <- rbind(
  performance(sample.0.5..5..2..1000),
  performance(sample.0.5..10..2..1000),
  performance(sample.0.5..15..2..1000),
  performance(sample.0.5..20..2..1000),
  performance(sample.0.5..25..2..1000),
  performance(sample.0.5..30..2..1000),
  performance(sample.0.5..35..2..1000),
  performance(sample.0.5..40..2..1000),
  performance(sample.0.5..45..2..1000),
  performance(sample.0.5..50..2..1000),
  performance(sample.0.5..55..2..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 2",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### change r vs err_dim

sample.0.5..5..0.5..1000.10 <- sample.new(err_dim = 0.5, r = 10)
sample.0.5..10..0.5..1000.10 <- sample.new(d_err = 10, err_dim = 0.5, r = 10)
sample.0.5..15..0.5..1000.10 <- sample.new(d_err = 15, err_dim = 0.5, r = 10)
sample.0.5..20..0.5..1000.10 <- sample.new(d_err = 20, err_dim = 0.5, r = 10)
sample.0.5..25..0.5..1000.10 <- sample.new(d_err = 25, err_dim = 0.5, r = 10)
sample.0.5..30..0.5..1000.10 <- sample.new(d_err = 30, err_dim = 0.5, r = 10)
sample.0.5..35..0.5..1000.10 <- sample.new(d_err = 35, err_dim = 0.5, r = 10)
sample.0.5..40..0.5..1000.10 <- sample.new(d_err = 40, err_dim = 0.5, r = 10)
sample.0.5..45..0.5..1000.10 <- sample.new(d_err = 45, err_dim = 0.5, r = 10)
sample.0.5..50..0.5..1000.10 <- sample.new(d_err = 50, err_dim = 0.5, r = 10)
sample.0.5..55..0.5..1000.10 <- sample.new(d_err = 55, err_dim = 0.5, r = 10)
result <- rbind(
  performance(sample.0.5..5..0.5..1000.10),
  performance(sample.0.5..10..0.5..1000.10),
  performance(sample.0.5..15..0.5..1000.10),
  performance(sample.0.5..20..0.5..1000.10),
  performance(sample.0.5..25..0.5..1000.10),
  performance(sample.0.5..30..0.5..1000.10),
  performance(sample.0.5..35..0.5..1000.10),
  performance(sample.0.5..40..0.5..1000.10),
  performance(sample.0.5..45..0.5..1000.10),
  performance(sample.0.5..50..0.5..1000.10),
  performance(sample.0.5..55..0.5..1000.10)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 0.5 || r = 10",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)


### change r vs err_dim

sample.0.5..5..1..1000.10 <- sample.new(err_dim = 1, r = 10)
sample.0.5..10..1..1000.10 <- sample.new(d_err = 10, err_dim = 1, r = 10)
sample.0.5..15..1..1000.10 <- sample.new(d_err = 15, err_dim = 1, r = 10)
sample.0.5..20..1..1000.10 <- sample.new(d_err = 20, err_dim = 1, r = 10)
sample.0.5..25..1..1000.10 <- sample.new(d_err = 25, err_dim = 1, r = 10)
sample.0.5..30..1..1000.10 <- sample.new(d_err = 30, err_dim = 1, r = 10)
sample.0.5..35..1..1000.10 <- sample.new(d_err = 35, err_dim = 1, r = 10)
sample.0.5..40..1..1000.10 <- sample.new(d_err = 40, err_dim = 1, r = 10)
sample.0.5..45..1..1000.10 <- sample.new(d_err = 45, err_dim = 1, r = 10)
sample.0.5..50..1..1000.10 <- sample.new(d_err = 50, err_dim = 1, r = 10)
sample.0.5..55..1..1000.10 <- sample.new(d_err = 55, err_dim = 1, r = 10)
result <- rbind(
  performance(sample.0.5..5..1..1000.10),
  performance(sample.0.5..10..1..1000.10),
  performance(sample.0.5..15..1..1000.10),
  performance(sample.0.5..20..1..1000.10),
  performance(sample.0.5..25..1..1000.10),
  performance(sample.0.5..30..1..1000.10),
  performance(sample.0.5..35..1..1000.10),
  performance(sample.0.5..40..1..1000.10),
  performance(sample.0.5..45..1..1000.10),
  performance(sample.0.5..50..1..1000.10),
  performance(sample.0.5..55..1..1000.10)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 1 || r = 10",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### change r vs err_dim

sample.0.5..5..2..1000.10 <- sample.new(err_dim = 2, r = 10)
sample.0.5..10..2..1000.10 <- sample.new(d_err = 10, err_dim = 2, r = 10)
sample.0.5..15..2..1000.10 <- sample.new(d_err = 15, err_dim = 2, r = 10)
sample.0.5..20..2..1000.10 <- sample.new(d_err = 20, err_dim = 2, r = 10)
sample.0.5..25..2..1000.10 <- sample.new(d_err = 25, err_dim = 2, r = 10)
sample.0.5..30..2..1000.10 <- sample.new(d_err = 30, err_dim = 2, r = 10)
sample.0.5..35..2..1000.10 <- sample.new(d_err = 35, err_dim = 2, r = 10)
sample.0.5..40..2..1000.10 <- sample.new(d_err = 40, err_dim = 2, r = 10)
sample.0.5..45..2..1000.10 <- sample.new(d_err = 45, err_dim = 2, r = 10)
sample.0.5..50..2..1000.10 <- sample.new(d_err = 50, err_dim = 2, r = 10)
sample.0.5..55..2..1000.10 <- sample.new(d_err = 55, err_dim = 2, r = 10)
result <- rbind(
  performance(sample.0.5..5..2..1000.10),
  performance(sample.0.5..10..2..1000.10),
  performance(sample.0.5..15..2..1000.10),
  performance(sample.0.5..20..2..1000.10),
  performance(sample.0.5..25..2..1000.10),
  performance(sample.0.5..30..2..1000.10),
  performance(sample.0.5..35..2..1000.10),
  performance(sample.0.5..40..2..1000.10),
  performance(sample.0.5..45..2..1000.10),
  performance(sample.0.5..50..2..1000.10),
  performance(sample.0.5..55..2..1000.10)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 2 || r = 10",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### change r vs err_dim

sample.0.5..5..3..1000.10 <- sample.new(err_dim = 3, r = 10)
sample.0.5..10..3..1000.10 <- sample.new(d_err = 10, err_dim = 3, r = 10)
sample.0.5..15..3..1000.10 <- sample.new(d_err = 15, err_dim = 3, r = 10)
sample.0.5..20..3..1000.10 <- sample.new(d_err = 20, err_dim = 3, r = 10)
sample.0.5..25..3..1000.10 <- sample.new(d_err = 25, err_dim = 3, r = 10)
sample.0.5..30..3..1000.10 <- sample.new(d_err = 30, err_dim = 3, r = 10)
sample.0.5..35..3..1000.10 <- sample.new(d_err = 35, err_dim = 3, r = 10)
sample.0.5..40..3..1000.10 <- sample.new(d_err = 40, err_dim = 3, r = 10)
sample.0.5..45..3..1000.10 <- sample.new(d_err = 45, err_dim = 3, r = 10)
sample.0.5..50..3..1000.10 <- sample.new(d_err = 50, err_dim = 3, r = 10)
sample.0.5..55..3..1000.10 <- sample.new(d_err = 55, err_dim = 3, r = 10)
result <- rbind(
  performance(sample.0.5..5..3..1000.10),
  performance(sample.0.5..10..3..1000.10),
  performance(sample.0.5..15..3..1000.10),
  performance(sample.0.5..20..3..1000.10),
  performance(sample.0.5..25..3..1000.10),
  performance(sample.0.5..30..3..1000.10),
  performance(sample.0.5..35..3..1000.10),
  performance(sample.0.5..40..3..1000.10),
  performance(sample.0.5..45..3..1000.10),
  performance(sample.0.5..50..3..1000.10),
  performance(sample.0.5..55..3..1000.10)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5 ||  d_err varies  ||  err_dim = 3 || r = 10",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### Varying d_err while V_err is 1

sample.1..5..0.01..1000 <- sample.new(V_err = 1)
sample.1..10..0.01..1000 <- sample.new(d_err = 10, V_err = 1)
sample.1..15..0.01..1000 <- sample.new(d_err = 15, V_err = 1)
sample.1..20..0.01..1000 <- sample.new(d_err = 20, V_err = 1)
sample.1..25..0.01..1000 <- sample.new(d_err = 25, V_err = 1)
sample.1..30..0.01..1000 <- sample.new(d_err = 30, V_err = 1)
sample.1..35..0.01..1000 <- sample.new(d_err = 35, V_err = 1)
sample.1..40..0.01..1000 <- sample.new(d_err = 40, V_err = 1)
sample.1..45..0.01..1000 <- sample.new(d_err = 45, V_err = 1)
sample.1..50..0.01..1000 <- sample.new(d_err = 50, V_err = 1)
sample.1..55..0.01..1000 <- sample.new(d_err = 55, V_err = 1)

result <- rbind(
  performance(sample.1..5..0.01..1000),
  performance(sample.1..10..0.01..1000),
  performance(sample.1..15..0.01..1000),
  performance(sample.1..20..0.01..1000),
  performance(sample.1..25..0.01..1000),
  performance(sample.1..30..0.01..1000),
  performance(sample.1..35..0.01..1000),
  performance(sample.1..40..0.01..1000),
  performance(sample.1..45..0.01..1000),
  performance(sample.1..50..0.01..1000),
  performance(sample.1..55..0.01..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 1  ||  d_err varies  ||  err_dim = 0.01",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### Varying d_err while V_err is 2

sample.2..5..0.01..1000 <- sample.new(V_err = 2)
sample.2..10..0.01..1000 <- sample.new(d_err = 10, V_err = 2)
sample.2..15..0.01..1000 <- sample.new(d_err = 15, V_err = 2)
sample.2..20..0.01..1000 <- sample.new(d_err = 20, V_err = 2)
sample.2..25..0.01..1000 <- sample.new(d_err = 25, V_err = 2)
sample.2..30..0.01..1000 <- sample.new(d_err = 30, V_err = 2)
sample.2..35..0.01..1000 <- sample.new(d_err = 35, V_err = 2)
sample.2..40..0.01..1000 <- sample.new(d_err = 40, V_err = 2)
sample.2..45..0.01..1000 <- sample.new(d_err = 45, V_err = 2)
sample.2..50..0.01..1000 <- sample.new(d_err = 50, V_err = 2)
sample.2..55..0.01..1000 <- sample.new(d_err = 55, V_err = 2)


result <- rbind(
  performance(sample.2..5..0.01..1000),
  performance(sample.2..10..0.01..1000),
  performance(sample.2..15..0.01..1000),
  performance(sample.2..20..0.01..1000),
  performance(sample.2..25..0.01..1000),
  performance(sample.2..30..0.01..1000),
  performance(sample.2..35..0.01..1000),
  performance(sample.2..40..0.01..1000),
  performance(sample.2..45..0.01..1000),
  performance(sample.2..50..0.01..1000),
  performance(sample.2..55..0.01..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 2  ||  d_err varies  ||  err_dim = 0.01",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### changing r
sample.0.5..5..1..1000.1 <- sample.new(err_dim = 1, r = 1)
sample.0.5..5..1..1000.2 <- sample.new(err_dim = 1, r = 2)
sample.0.5..5..1..1000.3 <- sample.new(err_dim = 1, r = 3)
sample.0.5..5..1..1000.4 <- sample.new(err_dim = 1, r = 4)
sample.0.5..5..1..1000.5 <- sample.new(err_dim = 1, r = 5)
sample.0.5..5..1..1000.6 <- sample.new(err_dim = 1, r = 6)
sample.0.5..5..1..1000.7 <- sample.new(err_dim = 1, r = 7)
sample.0.5..5..1..1000.8 <- sample.new(err_dim = 1, r = 8)
sample.0.5..5..1..1000.9 <- sample.new(err_dim = 1, r = 9)
sample.0.5..5..1..1000.10 <- sample.new(err_dim = 1, r = 10)
sample.0.5..5..1..1000.11 <- sample.new(err_dim = 1, r = 11)
sample.0.5..5..1..1000.12 <- sample.new(err_dim = 1, r = 12)
sample.0.5..5..1..1000.15 <- sample.new(err_dim = 1, r = 15)
sample.0.5..5..1..1000.20 <- sample.new(err_dim = 1, r = 20)
sample.0.5..5..1..1000.25 <- sample.new(err_dim = 1, r = 25)
sample.0.5..5..1..1000.30 <- sample.new(err_dim = 1, r = 30)
sample.0.5..5..1..1000.35 <- sample.new(err_dim = 1, r = 35)
sample.0.5..5..1..1000.40 <- sample.new(err_dim = 1, r = 40)
sample.0.5..5..1..1000.45 <- sample.new(err_dim = 1, r = 45)
sample.0.5..5..1..1000.50 <- sample.new(err_dim = 1, r = 50)
sample.0.5..5..1..1000.55 <- sample.new(err_dim = 1, r = 55)
sample.0.5..5..1..1000.60 <- sample.new(err_dim = 1, r = 60)
result <- rbind(
  performance(sample.0.5..5..1..1000.1),
  performance(sample.0.5..5..1..1000.2),
  performance(sample.0.5..5..1..1000.3),
  performance(sample.0.5..5..1..1000.4),
  performance(sample.0.5..5..1..1000.5),
  performance(sample.0.5..5..1..1000.6),
  performance(sample.0.5..5..1..1000.7),
  performance(sample.0.5..5..1..1000.8),
  performance(sample.0.5..5..1..1000.9),
  performance(sample.0.5..5..1..1000.10),
  performance(sample.0.5..5..1..1000.11),
  performance(sample.0.5..5..1..1000.12),
  performance(sample.0.5..5..1..1000.15),
  performance(sample.0.5..5..1..1000.20),
  performance(sample.0.5..5..1..1000.25),
  performance(sample.0.5..5..1..1000.30),
  performance(sample.0.5..5..1..1000.35),
  performance(sample.0.5..5..1..1000.40),
  performance(sample.0.5..5..1..1000.45),
  performance(sample.0.5..5..1..1000.50),
  performance(sample.0.5..5..1..1000.55),
  performance(sample.0.5..5..1..1000.60)
) 

plot(c(c(1:12), 5*c(3:12)), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 2  ||  d_err = 5  ||  err_dim = 1 || r varies",
     xlab = "r = 1, 2, ..., 12",
     ylab = "MSE")
points(c(c(1:12), 5*c(3:12)), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)



### zooming in r

for (i in seq(1, 5, by = 0.5)){
  assign(paste("sample.0.5..5..1..1000.", i,sep=""),sample.new(err_dim = 1, r = i))
}
result <- performance(sample.0.5..5..1..1000.1)
for (i in seq(1.5,5,by=0.5)){
  result <- rbind(result, performance( get(paste("sample.0.5..5..1..1000.", i,sep=""))  ))
}
plot(seq(1, 5, by = 0.5), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 10.5  ||  d_err = 5  ||  err_dim = 1 || r varies",
     xlab = "r = 1, 1.5, 2, ..., 5",
     ylab = "MSE")
points(seq(1, 5, by = 0.5), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)
  


### combining above two

result <- rbind(
  performance(sample.0.5..5..1..1000.1),
  performance(sample.0.5..5..1..1000.2),
  performance(sample.0.5..5..1..1000.3),
  performance(sample.0.5..5..1..1000.4),
  performance(sample.0.5..5..1..1000.5),
  performance(sample.0.5..5..1..1000.6),
  performance(sample.0.5..5..1..1000.7),
  performance(sample.0.5..5..1..1000.8),
  performance(sample.0.5..5..1..1000.9),
  performance(sample.0.5..5..1..1000.10),
  performance(sample.0.5..5..1..1000.11),
  performance(sample.0.5..5..1..1000.12),
  performance(sample.0.5..5..1..1000.15),
  performance(sample.0.5..5..1..1000.20),
  performance(sample.0.5..5..1..1000.25),
  performance(sample.0.5..5..1..1000.30),
  performance(sample.0.5..5..1..1000.35),
  performance(sample.0.5..5..1..1000.40),
  performance(sample.0.5..5..1..1000.45),
  performance(sample.0.5..5..1..1000.50),
  performance(sample.0.5..5..1..1000.55),
  performance(sample.0.5..5..1..1000.60)
) 

plot(c(c(1:12), 5*c(3:12)), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5  ||  err_dim = 1 || r varies",
     xlab = "r = 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 11, 12",
     ylab = "MSE")
points(c(c(1:12), 5*c(3:12)), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)



result <- performance(sample.0.5..5..1..1000.1)
for (i in seq(1.5,5,by=0.5)){
  result <- rbind(result, performance( get(paste("sample.0.5..5..1..1000.", i,sep=""))  ))
}
points(seq(1, 5, by = 0.5), result[,1], pch = 20)
points(seq(1, 5, by = 0.5), result[,2], col = "red", pch = 20)



### err_dim = 0.5, r varies

for (i in seq(1, 12, by = 0.5)){
  assign(paste("sample.0.5..5..0.5..1000..", i,sep=""),sample.new(err_dim = 0.5, r = i))
}
result <- performance(sample.0.5..5..0.5..1000..1)
for (i in seq(1.5,12,by=0.5)){
  result <- rbind(result, performance( get(paste("sample.0.5..5..0.5..1000..", i,sep=""))  ))
}
plot(seq(1, 12, by = 0.5), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5  ||  err_dim = 0.5 || r varies",
     xlab = "r = 1, 1.5, 2, ..., 5",
     ylab = "MSE")
points(seq(1, 12, by = 0.5), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)





result <- rbind(
  performance(sample.0.5..5..1..1000.1),
  performance(sample.0.5..5..1..1000.2),
  performance(sample.0.5..5..1..1000.3),
  performance(sample.0.5..5..1..1000.4),
  performance(sample.0.5..5..1..1000.5),
  performance(sample.0.5..5..1..1000.6),
  performance(sample.0.5..5..1..1000.7),
  performance(sample.0.5..5..1..1000.8),
  performance(sample.0.5..5..1..1000.9),
  performance(sample.0.5..5..1..1000.10),
  performance(sample.0.5..5..1..1000.11),
  performance(sample.0.5..5..1..1000.12)
) 

points(c(c(1:12)), 
     result[,1], 
     pch = 20,
     col = "purple")
points(c(c(1:12)), result[,2], col = "orange", pch = 20)

result <- performance(sample.0.5..5..1..1000.1)
for (i in seq(1.5,5,by=0.5)){
  result <- rbind(result, performance( get(paste("sample.0.5..5..1..1000.", i,sep=""))  ))
}
points(seq(1, 5, by = 0.5), result[,1], col = "purple", pch = 20)
points(seq(1, 5, by = 0.5), result[,2], col = "orange", pch = 20)


### d_err = 5, while err_dim varies

for (i in seq(0.01, 1, by = 0.01)){
  assign(paste("sample.0.5..5..", i, "..1000..1",sep=""), sample.new(err_dim = i))
}
for (i in seq(1.1, 2, by = 0.1)){
  assign(paste("sample.0.5..5..", i, "..1000..1",sep=""), sample.new(err_dim = i))
}
result <- performance(sample.0.5..5..0.01..1000..1)
for (i in c(seq(0.02,1,by=0.01),seq(1.1, 2, by = 0.1)) ){
  result <- rbind(result, performance( get(paste("sample.0.5..5..", i, "..1000..1",sep=""))  ))
}
plot(c(seq(0.01,1,by=0.01),seq(1.1, 2, by = 0.1)), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = c(rep(".",100), rep("o", 10)),
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.01, 0.02, ..., 0.99, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0",
     ylab = "MSE")
points(c(seq(0.01,1,by=0.01),seq(1.1, 2, by = 0.1)), result[,2], col = "red", pch= c(rep(".",100), rep("o", 10)))
lines(c(seq(0.01,1,by=0.01),seq(1.1, 2, by = 0.1)), result[,1])
lines(c(seq(0.01,1,by=0.01),seq(1.1, 2, by = 0.1)), result[,2], col = "red")
ticks <- 0.1*c(1:19)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:19), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### changing d_err to be higher than five, and varying err_dim

for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..10..", i, "..1000..1",sep=""), sample.new(d_err = 10, err_dim = i))
}
result <- performance(sample.0.5..10..0.1..1000..1)
for (i in seq(0.2, 2, by = 0.1) ){
  result <- rbind(result, performance( get(paste("sample.0.5..10..", i, "..1000..1",sep=""))  ))
}
plot(seq(0.1, 2, by = 0.1), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch =  rep("o", 20),
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 10  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.1, 0.2, ..., 1.9, 2.0",
     ylab = "MSE")
points(seq(0.1, 2, by = 0.1), result[,2], col = "red", pch = rep("o", 20))
lines(seq(0.1, 2, by = 0.1), result[,1])
lines(seq(0.1, 2, by = 0.1), result[,2], col = "red")
ticks <- 0.1*c(1:19)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:19), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

### changing d_err to be 20, and varying err_dim

for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..20..", i, "..1000..1",sep=""), sample.new(d_err = 20, err_dim = i))
}
result <- performance(sample.0.5..20..0.1..1000..1)
for (i in seq(0.2, 2, by = 0.1) ){
  result <- rbind(result, performance( get(paste("sample.0.5..20..", i, "..1000..1",sep=""))  ))
}
plot(seq(0.1, 2, by = 0.1), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch =  rep("o", 20),
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 20  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.1, 0.2, ..., 1.9, 2.0",
     ylab = "MSE")
points(seq(0.1, 2, by = 0.1), result[,2], col = "red", pch = rep("o", 20))
lines(seq(0.1, 2, by = 0.1), result[,1])
lines(seq(0.1, 2, by = 0.1), result[,2], col = "red")
ticks <- 0.1*c(1:19)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:19), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

###

for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..", tmp.count, "..", i, "..1000..1",sep=""), sample.new(d_err = 30, err_dim = i))
}
tmp.count = 30
result <- performance(get(paste("sample.0.5..", tmp.count, "..", "0.1", "..1000..1",sep=""))[1:20,])
for (i in seq(0.2, 2, by = 0.1) ){
  result <- rbind(result, performance(get(paste("sample.0.5..", tmp.count, "..", i, "..1000..1",sep=""))[1:200,]  ))
}
plot(seq(0.1, 2, by = 0.1), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch =  rep("o", 20),
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 30  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.1, 0.2, ..., 1.9, 2.0",
     ylab = "MSE")
points(seq(0.1, 2, by = 0.1), result[,2], col = "red", pch = rep("o", 20))
lines(seq(0.1, 2, by = 0.1), result[,1])
lines(seq(0.1, 2, by = 0.1), result[,2], col = "red")
ticks <- 0.1*c(1:19)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:19), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)


for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..40..", i, "..20..1",sep=""), sample.new(d_err = 40, err_dim = i, N_rep = 100))
}
tmp.count = 40
result <- performance(get(paste("sample.0.5..", tmp.count, "..", "0.1", "..20..1",sep="")))
for (i in seq(0.2, 2, by = 0.1) ){
  result <- rbind(result, performance(get(paste("sample.0.5..", tmp.count, "..", i, "..20..1",sep="")) ))
}
plot(seq(0.1, 2, by = 0.1), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch =  rep("o", 20),
     main = "N_rep = 100  ||  V_err = 0.5  ||  d_err = 40  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.1, 0.2, ..., 1.9, 2.0",
     ylab = "MSE")
points(seq(0.1, 2, by = 0.1), result[,2], col = "red", pch = rep("o", 20))
lines(seq(0.1, 2, by = 0.1), result[,1])
lines(seq(0.1, 2, by = 0.1), result[,2], col = "red")
ticks <- 0.1*c(1:19)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:19), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)




for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..50..", i, "..100..1",sep=""), sample.new(d_err = 50, err_dim = i, N_rep = 100))
}
tmp.count = 50
result <- performance(get(paste("sample.0.5..", tmp.count, "..", "0.1", "..100..1",sep="")))
for (i in seq(0.2, 2, by = 0.1) ){
  result <- rbind(result, performance(get(paste("sample.0.5..", tmp.count, "..", i, "..100..1",sep="")) ))
}
plot(seq(0.1, 2, by = 0.1), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch =  rep("o", 20),
     main = "N_rep = 100  ||  V_err = 0.5  ||  d_err = 50  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.1, 0.2, ..., 1.9, 2.0",
     ylab = "MSE")
points(seq(0.1, 2, by = 0.1), result[,2], col = "red", pch = rep("o", 20))
lines(seq(0.1, 2, by = 0.1), result[,1])
lines(seq(0.1, 2, by = 0.1), result[,2], col = "red")
ticks <- 0.1*c(1:19)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:19), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)


for (i in seq(0.1, 0.8, by = 0.1)){
  assign(paste("sample.0.5..60..", i, "..100..1",sep=""), sample.new(d_err = 60, err_dim = i, N_rep = 100))
}
tmp.count = 60
result <- performance(get(paste("sample.0.5..", tmp.count, "..", "0.1", "..100..1",sep="")))
for (i in seq(0.2, 0.8, by = 0.1) ){
  result <- rbind(result, performance(get(paste("sample.0.5..", tmp.count, "..", i, "..100..1",sep="")) ))
}

plot(seq(0.1, 0.8, by = 0.1), 
     result[,1], 
     ylim = c(0,0.1),
     xlim = c(0, 2.0),
     grid(),
     pch =  rep("o", 20),
     main = "N_rep = 100  ||  V_err = 0.5  ||  d_err = 60  ||  err_dim varies || r = 1",
     xlab = "err_dim = 0.1, 0.2, ..., 0.7. 0.8",
     ylab = "MSE")
points(seq(0.1, 0.8, by = 0.1), result[,2], col = "red", pch = rep("o", 20))
lines(seq(0.1, 0.8, by = 0.1), result[,1])
lines(seq(0.1, 0.8, by = 0.1), result[,2], col = "red")
ticks <- 0.1*c(1:7)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:7), col="gray", lty=3)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)


for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..70..", i, "..1000..1",sep=""), sample.new(d_err = 70, err_dim = i))
}
for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..90..", i, "..1000..1",sep=""), sample.new(d_err = 90, err_dim = i))
}
for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..110..", i, "..1000..1",sep=""), sample.new(d_err = 110, err_dim = i))
}
for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..130..", i, "..1000..1",sep=""), sample.new(d_err = 130, err_dim = i))
}
for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..150..", i, "..1000..1",sep=""), sample.new(d_err = 150, err_dim = i))
}
for (i in seq(0.1, 2, by = 0.1)){
  assign(paste("sample.0.5..170..", i, "..1000..1",sep=""), sample.new(d_err = 170, err_dim = i))
}


### Varying d_err while V_err is 10.5

sample.10.5..5..0.01..1000 <- sample.new(V_err = 10.5)
sample.10.5..10..0.01..1000 <- sample.new(d_err = 10, V_err = 10.5)
sample.10.5..15..0.01..1000 <- sample.new(d_err = 15, V_err = 10.5)
sample.10.5..20..0.01..1000 <- sample.new(d_err = 20, V_err = 10.5)
sample.10.5..25..0.01..1000 <- sample.new(d_err = 25, V_err = 10.5)
sample.10.5..30..0.01..1000 <- sample.new(d_err = 30, V_err = 10.5)
sample.10.5..35..0.01..1000 <- sample.new(d_err = 35, V_err = 10.5)
sample.10.5..40..0.01..1000 <- sample.new(d_err = 40, V_err = 10.5)
sample.10.5..45..0.01..1000 <- sample.new(d_err = 45, V_err = 10.5)
sample.10.5..50..0.01..1000 <- sample.new(d_err = 50, V_err = 10.5)
sample.10.5..55..0.01..1000 <- sample.new(d_err = 55, V_err = 10.5)

result <- rbind(
  performance(sample.10.5..5..0.01..1000),
  performance(sample.10.5..10..0.01..1000),
  performance(sample.10.5..15..0.01..1000),
  performance(sample.10.5..20..0.01..1000),
  performance(sample.10.5..25..0.01..1000),
  performance(sample.10.5..30..0.01..1000),
  performance(sample.10.5..35..0.01..1000),
  performance(sample.10.5..40..0.01..1000),
  performance(sample.10.5..45..0.01..1000),
  performance(sample.10.5..50..0.01..1000),
  performance(sample.10.5..55..0.01..1000)
)
plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 10.5  ||  d_err varies  ||  err_dim = 0.01",
     xlab = "d_err = 5, 10, 15, ..., 55",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)


### Varying err_dim

sample.0.5..5..0.0001..1000 <- sample.new(err_dim = 0.0001)
sample.0.5..5..0.01..1000 <- sample.new(err_dim = 0.01)
sample.0.5..5..0.02..1000 <- sample.new(err_dim = 0.02)
sample.0.5..5..0.03..1000 <- sample.new(err_dim = 0.03)
sample.0.5..5..0.04..1000 <- sample.new(err_dim = 0.04)
sample.0.5..5..0.05..1000 <- sample.new(err_dim = 0.05)
sample.0.5..5..0.06..1000 <- sample.new(err_dim = 0.06)
sample.0.5..5..0.07..1000 <- sample.new(err_dim = 0.07)
sample.0.5..5..0.08..1000 <- sample.new(err_dim = 0.08)
sample.0.5..5..0.09..1000 <- sample.new(err_dim = 0.09)
sample.0.5..5..0.10..1000 <- sample.new(err_dim = 0.10)
result <- rbind(
  performance(sample.0.5..5..0.0001..1000),
  performance(sample.0.5..5..0.01..1000),
  performance(sample.0.5..5..0.02..1000),
  performance(sample.0.5..5..0.03..1000),
  performance(sample.0.5..5..0.04..1000),
  performance(sample.0.5..5..0.05..1000),
  performance(sample.0.5..5..0.06..1000),
  performance(sample.0.5..5..0.07..1000),
  performance(sample.0.5..5..0.08..1000),
  performance(sample.0.5..5..0.09..1000),
  performance(sample.0.5..5..0.10..1000)
)

plot(c(0.0001, 0.01*c(1:10)), 
     result[,1], 
     ylim = c(0,0.03),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err = 0.5  ||  d_err = 5  ||  err_dim varies",
     xlab = "err_dim  = 0.0001, 0.01, 0.02, ..., 0.08, 0.09, 0.10",
     ylab = "MSE")
points(c(0.0001, 0.01*c(1:10)), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)



### Varying v_err


sample.0.5..5..0.01..1000 <- sample.new(V_err = 0.5)
sample.1.5..5..0.01..1000 <- sample.new(V_err = 1.5)
sample.2.5..5..0.01..1000 <- sample.new(V_err = 2.5)
sample.3.5..5..0.01..1000 <- sample.new(V_err = 3.5)
sample.4.5..5..0.01..1000 <- sample.new(V_err = 4.5)
sample.5.5..5..0.01..1000 <- sample.new(V_err = 5.5)
sample.6.5..5..0.01..1000 <- sample.new(V_err = 6.5)
sample.7.5..5..0.01..1000 <- sample.new(V_err = 7.5)
sample.8.5..5..0.01..1000 <- sample.new(V_err = 8.5)
sample.9.5..5..0.01..1000 <- sample.new(V_err = 9.5)
sample.10.5..5..0.01..1000 <- sample.new(V_err = 10.5)

result <- rbind(
  performance(sample.0.5..5..0.01..1000),
  performance(sample.1.5..5..0.01..1000),
  performance(sample.2.5..5..0.01..1000),
  performance(sample.3.5..5..0.01..1000),
  performance(sample.4.5..5..0.01..1000),
  performance(sample.5.5..5..0.01..1000),
  performance(sample.6.5..5..0.01..1000),
  performance(sample.7.5..5..0.01..1000),
  performance(sample.8.5..5..0.01..1000),
  performance(sample.9.5..5..0.01..1000),
  performance(sample.10.5..5..0.01..1000)
)


plot(0.5+1*c(0:10), 
     result[,1], 
     ylim = c(0,0.1),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err varies  ||  d_err = 5  ||  err_dim = 0.01",
     xlab = "V_err = 0.5, 1.5, 2.5, ..., 9.5, 10.5",
     ylab = "MSE")
points(0.5+1*c(0:10), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)



### Varying V_err
sample.0.05..5..0.01..1000 <- sample.new(V_err = 0.05)
sample.0.10..5..0.01..1000 <- sample.new(V_err = 0.10)
sample.0.15..5..0.01..1000 <- sample.new(V_err = 0.15)
sample.0.20..5..0.01..1000 <- sample.new(V_err = 0.20)
sample.0.25..5..0.01..1000 <- sample.new(V_err = 0.25)
sample.0.30..5..0.01..1000 <- sample.new(V_err = 0.30)
sample.0.35..5..0.01..1000 <- sample.new(V_err = 0.35)
sample.0.40..5..0.01..1000 <- sample.new(V_err = 0.40)
sample.0.45..5..0.01..1000 <- sample.new(V_err = 0.45)
sample.0.50..5..0.01..1000 <- sample.new(V_err = 0.50)
sample.0.55..5..0.01..1000 <- sample.new(V_err = 0.55)
sample.0.60..5..0.01..1000 <- sample.new(V_err = 0.60)
sample.0.65..5..0.01..1000 <- sample.new(V_err = 0.65)
sample.0.70..5..0.01..1000 <- sample.new(V_err = 0.70)
sample.0.75..5..0.01..1000 <- sample.new(V_err = 0.75)
sample.0.80..5..0.01..1000 <- sample.new(V_err = 0.80)
sample.0.85..5..0.01..1000 <- sample.new(V_err = 0.85)
sample.0.90..5..0.01..1000 <- sample.new(V_err = 0.90)
sample.0.95..5..0.01..1000 <- sample.new(V_err = 0.95)
sample.1.00..5..0.01..1000 <- sample.new(V_err = 1.0)
sample.1.05..5..0.01..1000 <- sample.new(V_err = 1.05)
sample.1.10..5..0.01..1000 <- sample.new(V_err = 1.1)
sample.1.15..5..0.01..1000 <- sample.new(V_err = 1.15)
sample.1.20..5..0.01..1000 <- sample.new(V_err = 1.2)
sample.1.25..5..0.01..1000 <- sample.new(V_err = 1.25)
sample.1.30..5..0.01..1000 <- sample.new(V_err = 1.3)
sample.1.35..5..0.01..1000 <- sample.new(V_err = 1.35)
sample.1.40..5..0.01..1000 <- sample.new(V_err = 1.4)
sample.1.45..5..0.01..1000 <- sample.new(V_err = 1.45)
sample.1.50..5..0.01..1000 <- sample.new(V_err = 1.5)

result <- rbind(
  performance(sample.0.05..5..0.01..1000),
  performance(sample.0.10..5..0.01..1000),
  performance(sample.0.15..5..0.01..1000),
  performance(sample.0.20..5..0.01..1000),
  performance(sample.0.25..5..0.01..1000),
  performance(sample.0.30..5..0.01..1000),
  performance(sample.0.35..5..0.01..1000),
  performance(sample.0.40..5..0.01..1000),
  performance(sample.0.45..5..0.01..1000),
  performance(sample.0.50..5..0.01..1000),
  performance(sample.0.55..5..0.01..1000),
  performance(sample.0.60..5..0.01..1000),
  performance(sample.0.65..5..0.01..1000),
  performance(sample.0.70..5..0.01..1000),
  performance(sample.0.75..5..0.01..1000),
  performance(sample.0.80..5..0.01..1000),
  performance(sample.0.85..5..0.01..1000),
  performance(sample.0.90..5..0.01..1000),
  performance(sample.0.95..5..0.01..1000),
  performance(sample.1.00..5..0.01..1000),
  performance(sample.1.05..5..0.01..1000),
  performance(sample.1.10..5..0.01..1000),
  performance(sample.1.15..5..0.01..1000),
  performance(sample.1.20..5..0.01..1000),
  performance(sample.1.25..5..0.01..1000),
  performance(sample.1.30..5..0.01..1000),
  performance(sample.1.35..5..0.01..1000),
  performance(sample.1.40..5..0.01..1000),
  performance(sample.1.45..5..0.01..1000),
  performance(sample.1.50..5..0.01..1000)
)

plot(0.05*c(1:30), 
     result[,1], 
     ylim = c(0,0.1),
     # xlim = c(0, 57),
     # grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err varies  ||  d_err = 5  ||  err_dim = 0.01",
     xlab = "V_err = 0.05, 0.10, ..., 1.45, 1.50",
     ylab = "MSE")
ticks <- 0.1*c(1:15)
axis(side = 1, at = ticks)
abline(h=0.01*(0:10), v=0.1*c(1:15), col="gray", lty=3)
points(0.05*c(1:30), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)


### Varying V_err from 5 through 55
sample.5..5..0.01..1000 <- sample.new(V_err = 5)
sample.10..5..0.01..1000 <- sample.new(V_err = 10)
sample.15..5..0.01..1000 <- sample.new(V_err = 15)
sample.20..5..0.01..1000 <- sample.new(V_err = 20)
sample.25..5..0.01..1000 <- sample.new(V_err = 25)
sample.30..5..0.01..1000 <- sample.new(V_err = 30)
sample.35..5..0.01..1000 <- sample.new(V_err = 35)
sample.40..5..0.01..1000 <- sample.new(V_err = 40)
sample.45..5..0.01..1000 <- sample.new(V_err = 45)
sample.50..5..0.01..1000 <- sample.new(V_err = 50)
sample.55..5..0.01..1000 <- sample.new(V_err = 55)

result <- rbind(
  performance(sample.5..5..0.01..1000),
  performance(sample.10..5..0.01..1000),
  performance(sample.15..5..0.01..1000),
  performance(sample.20..5..0.01..1000),
  performance(sample.25..5..0.01..1000),
  performance(sample.30..5..0.01..1000),
  performance(sample.35..5..0.01..1000),
  performance(sample.40..5..0.01..1000),
  performance(sample.45..5..0.01..1000),
  performance(sample.50..5..0.01..1000),
  performance(sample.55..5..0.01..1000)
)

plot(5*c(1:11), 
     result[,1], 
     ylim = c(0,0.1),
     xlim = c(0, 57),
     grid(),
     pch = 20,
     main = "N_rep = 1000  ||  V_err varies  ||  d_err = 5  ||  err_dim = 0.01",
     xlab = "V_err = 0.05 to 1.5 by 0.05, 0.5 to 10.5 by 1, 5 to 55 by 5",
     ylab = "MSE")
points(5*c(1:11), result[,2], col = "red", pch = 20)
legend('bottomright',
       bg="transparent",
       c("naive",  "importance"),
       # lty = c("solid",  "solid"),
       col = c("red",   "black"),
       pch = 20,
       # bty = 'y',
       # lwd = 2,
       # title = "Naive Sampling",
       # title.col = "black",
       cex = 1)

result <- rbind(
  performance(sample.0.5..5..0.01..1000),
  performance(sample.1.5..5..0.01..1000),
  performance(sample.2.5..5..0.01..1000),
  performance(sample.3.5..5..0.01..1000),
  performance(sample.4.5..5..0.01..1000),
  performance(sample.5.5..5..0.01..1000),
  performance(sample.6.5..5..0.01..1000),
  performance(sample.7.5..5..0.01..1000),
  performance(sample.8.5..5..0.01..1000),
  performance(sample.9.5..5..0.01..1000),
  performance(sample.10.5..5..0.01..1000)
)
points(0.5+1*c(0:10), result[,1], pch = 20)
points(0.5+1*c(0:10), result[,2], col = "red", pch = 20)



result <- rbind(
  performance(sample.0.05..5..0.01..1000),
  performance(sample.0.10..5..0.01..1000),
  performance(sample.0.15..5..0.01..1000),
  performance(sample.0.20..5..0.01..1000),
  performance(sample.0.25..5..0.01..1000),
  performance(sample.0.30..5..0.01..1000),
  performance(sample.0.35..5..0.01..1000),
  performance(sample.0.40..5..0.01..1000),
  performance(sample.0.45..5..0.01..1000),
  performance(sample.0.50..5..0.01..1000),
  performance(sample.0.55..5..0.01..1000),
  performance(sample.0.60..5..0.01..1000),
  performance(sample.0.65..5..0.01..1000),
  performance(sample.0.70..5..0.01..1000),
  performance(sample.0.75..5..0.01..1000),
  performance(sample.0.80..5..0.01..1000),
  performance(sample.0.85..5..0.01..1000),
  performance(sample.0.90..5..0.01..1000),
  performance(sample.0.95..5..0.01..1000),
  performance(sample.1.00..5..0.01..1000),
  performance(sample.1.05..5..0.01..1000),
  performance(sample.1.10..5..0.01..1000),
  performance(sample.1.15..5..0.01..1000),
  performance(sample.1.20..5..0.01..1000),
  performance(sample.1.25..5..0.01..1000),
  performance(sample.1.30..5..0.01..1000),
  performance(sample.1.35..5..0.01..1000),
  performance(sample.1.40..5..0.01..1000),
  performance(sample.1.45..5..0.01..1000),
  performance(sample.1.50..5..0.01..1000)
)

points(0.05*c(1:30), 
     result[,1], 
     ylim = c(0,0.1),
     pch = 20)
points(0.05*c(1:30), result[,2], col = "red", pch = 20)
