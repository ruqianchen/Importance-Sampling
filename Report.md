# Report 

## Intro

We designed a new importance sampling algorithm that is based on the naive sampling method. 


## Algorithm

We modify the naive sampling algorithm using KNN.


    # Setting up libraries and parameters
    
    library(FNN)
    library(ggplot2)
    V_err = 0.5
    d_err = 5
    err_dim = 0.01
    N_rep = 100
    n = 400 
    m = 200
    r = 1
    g_function = function(x){
      x>1.5
    }
    mu0 = 0.4092262



    # Sampling Angles
    
    U =  runif(m, min = 0, max = 2*pi)
    hist(
      U,
      col  = "purple",
      breaks = c(0:(ceiling(max(U))*2))*0.5,
      # prob = TRUE,
      # xaxt = 'n',
      xlab = "Component Index",
      main = "Histogram of Sampled Angles \n from the Uniform(0, 2pi) Distribution"
    )
![Angle sampled uniformly](https://d2mxuefqeaa7sj.cloudfront.net/s_048B61035C48B8887CB0369CD7CE976103AF89FC9C12F081E1F6FBC52F5989A8_1526668256815_image.png)

    X0 = r* cos(U)
    Y0 = r* sin(U)
    circle = data.frame(X0, Y0)
    
    plot(X0, 
         Y0, 
         main = "Sampled Points on the Unit Circle \n angle ~ Unif(0, 2pi) \n r = 1", 
         col = "purple", 
         pch = 20)
    grid()


![A circle resulting from the angles sampled above](https://d2mxuefqeaa7sj.cloudfront.net/s_048B61035C48B8887CB0369CD7CE976103AF89FC9C12F081E1F6FBC52F5989A8_1526668397896_image.png)



    # Adding Noise
    
    D = cbind(X0,Y0, matrix(rep(0, m*d_err), ncol=d_err))+
      cbind(matrix(rnorm(2*m, sd=0.1),ncol=2), 
            matrix(rnorm(d_err*m,sd=err_dim), ncol=d_err))
    
    plot(data.frame(D),
         col = "purple",
         pch = ".",
         main = " \n  \n ")
    title("Circle with Noise on Structural Dimensions \n and on Redudant Dimensions", cex.main=1)


![A seven dimension object, where the structural dimensions are the first two dimensions. The rest five dimensions are purely noise dimensions. Note the different axis scales in the first two dimensions versus the last five dimensions.](https://d2mxuefqeaa7sj.cloudfront.net/s_048B61035C48B8887CB0369CD7CE976103AF89FC9C12F081E1F6FBC52F5989A8_1526668422919_image.png)



    V = rnorm(m, mean=exp(D[,1]+D[,2]), sd=V_err)
    
    hist(
      V,
      col  = "purple",
      breaks = c((floor(min(V))*2):(ceiling(max(V))*2))*0.5,
      # prob = TRUE,
      # xaxt = 'n',
      # xlab = "Component Index",
      main = "Histogram of V"
    )


![](https://d2mxuefqeaa7sj.cloudfront.net/s_048B61035C48B8887CB0369CD7CE976103AF89FC9C12F081E1F6FBC52F5989A8_1526670591896_image.png)



    sp <- ggplot()  +
      ggtitle("Scatterplot of (X, Y) colored by V") + 
      geom_point(data = circle, 
                 aes(x = X0,
                     y = Y0,
                     color = V)) + 
      scale_colour_gradient(low="purple", high="white")
    sp 
    
![A plot of the training data without noise added yet. The lighter values correspond to regions where $$\exp(x_1+y_1)$$ are higher, e.g. the first quadrant.](https://d2mxuefqeaa7sj.cloudfront.net/s_048B61035C48B8887CB0369CD7CE976103AF89FC9C12F081E1F6FBC52F5989A8_1526668462742_image.png)




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


![A plot of the test data. The darker values correspond to regions where $$\exp(x_1+y_1)$$ are higher, e.g. the first quadrant.](https://d2mxuefqeaa7sj.cloudfront.net/s_048B61035C48B8887CB0369CD7CE976103AF89FC9C12F081E1F6FBC52F5989A8_1526668541305_image.png)



    # Importance
    p_max = max(sqrt(D1_fit$pred))
    
    idx_pass = runif(length(U))< (sqrt(D1_fit$pred)/p_max)
    
    D_pass = D1[idx_pass,]
    p_pass = (sqrt(D1_fit$pred)/p_max)[idx_pass]
    
    D_use = D_pass[1:n,]
    p_use = p_pass[1:n]
    V_use = rnorm(n, mean=exp(D_use[,1]+D_use[,2]), sd=V_err)
    
    mu_importance_sampling = mean(g_function(V_use)/p_use)*D1_fit_sum


    # Naive 
    U = runif(n, max=2*pi)
    X2 = r* cos(U)
    Y2 = r* sin(U)
    
    D2 = cbind(X2,Y2, matrix(rep(0, n*d_err), ncol=d_err))+
      cbind(matrix(rnorm(2*n, sd=0.1),ncol=2), 
            matrix(rnorm(d_err*n,sd=err_dim), ncol=d_err))
    
    V2 = rnorm(n, mean=exp(D2[,1]+D2[,2]), sd=V_err)
    
    mu = mean(g_function(V2))
    mu_naive_sampling = mu


## Performance

There are several parameters to vary in the algorithm. These include

- `V_err`
- `d_err`
- `err_dim`
- `n`
- `m`
- `N_rep`

We measure the performance using mean-squared error. The mean is obtained from


## Experiments

**Varying** $$err_{dim}$$

When we vary $$err_{dim}$$, we are changing the variance of coordinates in the error dimension. We notice that for small $$err_{dim}$$, the importance sampler outperforms the naive sampler, and the gap remains almost constant even when the error spans across fifty dimensions. 

For example, in the following three graphs, $$V_{err}$$ are chosen among 0.5, 1 and 2, while $$err_{dim}$$ is fixed at 0.01. We vary $$d_{err}$$ and see almost constant performance gap within each plot. This shows that when $$err_{dim}$$ is small, the importance sampler performs better than the naive sampler, even for high dimensions.




![N_rep = 1000 V_err = 0.5 d_err = 5, 10, 15, 20, …, 55 err_dim = 0.01](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1525986832047_image.png)



![N_rep = 1000 V_err = 1 d_err = 5, 10, 15, 20, …, 55 err_dim = 0.01](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1525990694525_image.png)



![N_rep = 1000 V_err = 2 d_err = 5, 10, 15, 20, …, 55 err_dim = 0.01](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1525991163527_image.png)


When $$err_{dim}$$ is even smaller, we see importance sampler outperforms the naive sampler.

![N_rep = 1000 V_err = 0.5 d_err = 5 err_dim = 0.0001, 0.01, 0.02, 0.03, …, 0.09, 0.10](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1525987128226_image.png)



We start to see more interesting behavior when the $$err_{dim}$$ is larger. For example, we can set $$err_{dim}$$ to be larger than $$1$$. The structural dimension is a circle of radius 1. The importance sampling estimator performs better than the naive sampler, until $$err_{dim}$$ gets larger (e.g. when its larger than the radius in the structural dimension). As we can see below, as the number of noise dimensions get higher, the performances of the two samplers cross over at smaller $$err_{dim}$$’s.


![N_rep = 1000 V_err = 0.5 d_err =  10 err_dim = varies r = 1](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526458778097_image.png)



![N_rep = 1000 V_err = 0.5 d_err =  20 err_dim = varies r = 1](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526460397658_image.png)



![N_rep = 1000 V_err = 0.5 d_err =  30 err_dim = varies r = 1](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526496413459_image.png)


 

![N_rep = 100 V_err = 0.5 d_err =  40 err_dim = varies r = 1](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526499352053_image.png)


 

![N_rep = 100 V_err = 0.5 d_err =  50 err_dim = varies r = 1](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526499324793_image.png)


 

![N_rep = 100 V_err = 0.5 d_err =  60 err_dim = varies r = 1](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526499450139_image.png)




**Varying** $$V_{err}$$
When we vary $$V_{err}$$, we are changing how noisy the $$Y|X$$ dependency is.  We expect that as the error increases, performances get worse for both estimators. For example,

![N_rep = 1000                 V_err = 0.05, 0.10, 0.15, …, 1.45, 1.50                 d_err = 5                 err_dim = 0.01](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1525997329453_image.png)


Here the gap between the two estimators gets closer as $$V_{err}$$ increases. However, the importance sampler still outperforms the naive sampler throughout the process. 

As we keep increasing $$V_{err}$$, we observe an unexpected trend of MSE, where the importance sampler’s performance increases and further outperforms the naive sampler. We are unsure of the cause of such improvement. 

![](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1525991277799_image.png)




**Varying** $$n$$

The variable $$n$$ controls how many samples in total we take for each sampler. It is also the number of times we evaluate $$Y$$. Hence $$n$$ measures our computational budget. We expect that as the computational budget gets larger, the sampler performance gets better.

This is indeed the case. For example, here we see that both samplers’ performance improve as $$n$$ gets larger.

![N_rep = 1000 V_err = 0.5 d_err = 5 err_dim = 0.01 n = 400, 600, 800, 1000, 1200, 1400, 2000, 3000 m = n/2](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526059984391_image.png)


**Varying** $$m$$ **while holding** $$n$$ **constant**
The variable $$m$$ is the number of times we sample at the first stage, while $$n-m$$ is the number of times we sample in the second stage. We see that in the experiment below, there is a crossover value of $$m$$, before which importance sampler outperforms the naive sampler.


![](https://d2mxuefqeaa7sj.cloudfront.net/s_3683266ED3E0611CD4F2BF11624F679A554CA6A5FD7B592307EEBAA43D0DB668_1526061069077_image.png)





