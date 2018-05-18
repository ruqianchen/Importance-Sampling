library(FNN)
library(ggplot2)
V_err = 0.5
d_err = 5
err_dim = 0.01
N_rep = 100
n = 400 
m = 200
r = 1

U = runif(m, min = 0, max = 2*pi)



hist(
  U,
  col  = "purple",
  breaks = c(0:(ceiling(max(U))*2))*0.5,
  # prob = TRUE,
  # xaxt = 'n',
  xlab = "Component Index",
  main = "Histogram of Sampled Angles \n from the Uniform(0, 2pi) Distribution"
)


X0 = r* cos(U)
Y0 = r* sin(U)
circle = data.frame(X0, Y0)

plot(X0, 
     Y0, 
     main = "Sampled Points on the Unit Circle \n angle ~ Unif(0, 2pi) \n r = 1", 
     col = "purple", 
     pch = 20)
grid()


D = cbind(X0,Y0, matrix(rep(0, m*d_err), ncol=d_err))+
  cbind(matrix(rnorm(2*m, sd=0.1),ncol=2), 
        matrix(rnorm(d_err*m,sd=err_dim), ncol=d_err))


plot(data.frame(D),
     col = "purple",
     pch = ".",
     main = " \n  \n ")
title("Circle with Noise on Structural Dimensions \n and on Redudant Dimensions", cex.main=1)


V = rnorm(m, mean=exp(D[,1]+D[,2]), sd=V_err)


hist(
  V,
  col  = "purple",
  breaks = c((floor(min(V))*2):(ceiling(max(V))*2))*0.5,
  # prob = TRUE,
  # xaxt = 'n',
  xlab = "Component Index",
  main = "Histogram of V"
)


sp <- ggplot()  +
  ggtitle("Scatterplot of (X, Y) colored by V") + 
  geom_point(data = circle, 
             aes(x = X0,
                 y = Y0,
                 color = V)) + 
  scale_colour_gradient(low="purple", high="white")
sp 




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
# colors.grayscale <- paste("gray",  D1_fit$pred*100, sep = "")
color.tmp <- 101 - 100 * D1_fit$pred
colors.grayscale <-
  paste("gray", floor(100 -  3/4 * (100 * D1_fit$pred + 30)), sep = "")
plot(
  D1[,c(1,2)],
  col = colors.grayscale,
  bg = colors.grayscale,
  # col = gray.colors(100)[color.tmp],
  # bg = gray.colors(100)[color.tmp],
  main = "Samples for KNN Stage, colored by D1",
  pch = 20
)
points(D1[,c(1,2)], col = "purple")
grid()

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
