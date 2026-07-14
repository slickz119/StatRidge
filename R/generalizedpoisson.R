MCSRPois<-function(n, p, rho,beta_0){
  MSE<-rep(1,1000)
  MSE1<-rep(1,1000)
  MSE2<-rep(1,1000)
  MSE3<-rep(1,1000)
  MSE4<-rep(1,1000)
  
  for (r in 1:1000){
  B<-matrix(1:p,p,1)
    for (i in 1:p){B[i,1]<-sqrt(1/p)}
  w<-rnorm(200+n*p) #generate the 200 + n*p random numbers
  z<- matrix( w[seq(201,200+n*p)], n, p) #disregard the first 200 random numbers in w 
  x<-((1-rho^2)^(1/2))*z+(rho*z[,p]) #generating p predictors with multicollinearity
  mu<-exp(beta_0+x%*%B) 
    data <- data.frame(y = round(mu,0) , x) #creating data frame x and y
    # Define the design matrix X and the response vector y
    X <- as.matrix(data.frame(1, x))  # Adding intercept term
    y <- data$y 
    
    # Initialize beta
  beta1 <- matrix(0, nrow=ncol(X), ncol=1)
   #Compute mu.hat
  mu.hat<-exp(X%*%beta1)    
  
  # Computing k
  
  #to calculate alpha.hat and alpha.hat.max:
  xtx<- crossprod(x)
  Eigen<-eigen(xtx)
  Q<-Eigen$vectors
  Q<-t(Q)
  beta.ols=coef(lm(y~x-1)) # beta hat without intercept
  alpha.hat=Q%*%beta.ols
  alpha.hat.max<-max(alpha.hat^2)
  mu.hat1<-exp(x%*%beta.ols)  
  sigma.hat.sq= (sum((data$y-mu.hat)^2))/(n-p-1)
  
  # to compute K1
  
  k1 = sigma.hat.sq/alpha.hat.max
  
  #to compute MSE using k1
 
  #  learning rate alpha
  alpha <- 0.01

    
  # Gradient of the penalized log-likelihood function
  gradient <- function(X, y, beta1, k1) {
    mu <- exp(X %*% beta1)
    t(X) %*% (y - mu) - 2 * k1 * beta1
  }
  
  # Define the number of iterations
  iterations <- 1000
  
  # Gradient descent to update beta1
  for (i in 1:iterations) {
    grad <- gradient(X, y, beta1, k1)
    beta1 <- beta1 + alpha * grad
        }
  
  
  # Final coefficient
  beta1
  
    #compute for MSE
  MSE1[r]=sum(beta1-c(beta_0,B)^2)
}}
  