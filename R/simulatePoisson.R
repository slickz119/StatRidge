simulatePoisson <- function(n, p, rho, beta_0) {
  MSE1 <- rep(1, 1000)
  
  # Gradient of the penalized log-likelihood function
  gradient <- function(X, y, beta1, k1) {
    mu <- exp(X %*% beta1)
    mu[mu > 1e10] <- 1e10 # Limit mu to prevent overflow
    t(X) %*% (y - mu) - 2 * k1 * beta1
  }
  
  # Learning rate alpha
  alpha <- 0.01
  # Define the number of iterations
  iterations <- 1000
  
  for (r in 1:1000) {
    B <- matrix(1:p, p, 1)
    for (i in 1:p) {
      B[i, 1] <- sqrt(1 / p)
    }
    w <- rnorm(200 + n * p) # Generate the 200 + n*p random numbers
    z <- matrix(w[seq(201, 200 + n * p)], n, p) # Disregard the first 200 random numbers in w 
    x <- ((1 - rho^2)^(1/2)) * z + (rho * z[, p]) # Generating p predictors with multicollinearity
    mu <- exp(beta_0 + x %*% B) 
    mu[mu > 1e10] <- 1e10 # Limit mu to prevent overflow
    data <- data.frame(y = round(mu, 0), x) # Creating data frame x and y
    
    # Define the design matrix X and the response vector y
    X <- as.matrix(data.frame(1, x))  # Adding intercept term
    y <- data$y 
    
    # Initialize beta
    beta1 <- matrix(0, nrow = ncol(X), ncol = 1)
    # Compute mu.hat
    mu.hat <- exp(X %*% beta1)
    mu.hat[mu.hat > 1e10] <- 1e10 # Limit mu.hat to prevent overflow
    
    # Compute k
    xtx <- crossprod(x)
    Eigen <- eigen(xtx)
    Q <- Eigen$vectors
    Q <- t(Q)
    beta.ols <- coef(lm(y ~ x - 1)) # Beta hat without intercept
    alpha.hat <- Q %*% beta.ols
    alpha.hat.max <- max(alpha.hat^2)
    mu.hat1 <- exp(x %*% beta.ols)
    mu.hat1[mu.hat1 > 1e10] <- 1e10 # Limit mu.hat1 to prevent overflow
    
    sigma.hat.sq <- (sum((data$y - mu.hat)^2)) / (n - p - 1)
    if (is.nan(sigma.hat.sq) || sigma.hat.sq == 0) {
      print(paste("Iteration", r, ": sigma.hat.sq is NaN or zero"))
      next
    }
    
    # Compute k1
    k1 <- sigma.hat.sq / alpha.hat.max
    if (is.nan(k1) || k1 == 0) {
      print(paste("Iteration", r, ": k1 is NaN or zero"))
      next
    }
    
    # Gradient descent to update beta1
    for (i in 1:1000) {
      grad <- gradient(X, y, beta1, k1)
      if (any(is.nan(grad))) {
        print(paste("Iteration", r, ": gradient contains NaN at step", i))
        print("Gradient values:")
        print(grad)
        print("X values:")
        print(X)
        print("y values:")
        print(y)
        print("beta1 values:")
        print(beta1)
        next
      }
      beta1 <- beta1 + alpha * grad
    }
    
    if (any(is.nan(beta1))) {
      print(paste("Iteration", r, ": beta1 contains NaN"))
      print("beta1 values:")
      print(beta1)
      next
    }
    
    # Final coefficient
    BTrue <- matrix(c(beta_0, B), p + 1, 1)
    # Compute MSE
    MSE1[r] <- t(beta1 - BTrue) %*% (beta1 - BTrue)
  }
  
  if (any(is.nan(MSE1))) {
    print("MSE1 contains NaN values")
  }
  
  mean(MSE1, na.rm = TRUE)
}
