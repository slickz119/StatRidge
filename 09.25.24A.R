# Logistic Ridge Regression function using manual optimization
logistic_ridge_regression <- function(X, y, lambda) {
  # Loss function for logistic regression with ridge regularization
  loss_function <- function(beta) {
    # Compute the predicted probabilities using the logistic function
    p <- 1 / (1 + exp(-X %*% beta))
    
    # Clamp probabilities to avoid log(0) issues
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    
    # Calculate the negative log-likelihood with ridge penalty (L2 regularization)
    neg_log_likelihood <- -sum(y * log(p) + (1 - y) * log(1 - p))
    penalty <- lambda * sum(beta[-1]^2)  # Ridge regularization term, excluding the intercept
    return(neg_log_likelihood + penalty)
  }
  
  # Initialize beta parameters (initial guess)
  beta_initial <- rep(0, ncol(X))
  
  # Ensure that there are no non-finite values in beta_initial
  if (any(!is.finite(beta_initial))) {
    stop("Non-finite initial beta values")
  }
  
  # Call optim to minimize the loss function using BFGS
  result <- optim(beta_initial, loss_function, method = "BFGS", control = list(maxit = 1000))
  
  # Check for non-finite results in the optimization output
  if (any(!is.finite(result$par))) {
    stop("Optimization resulted in non-finite values")
  }
  
  return(result$par)  # Return the optimized beta parameters
}



# Monte Carlo Simulation for Logistic Ridge Regression
MCSRLr <- function(n, p, rho, beta_0) {
  num_replications <- 1000
  MSE0 <- numeric(num_replications)
  MSE1 <- numeric(num_replications)
  MSE2 <- numeric(num_replications)
  MSE3 <- numeric(num_replications)
  MSE4 <- numeric(num_replications)
  MSE5 <- numeric(num_replications)
  MSE6 <- numeric(num_replications)
  MSE7 <- numeric(num_replications)
  MSE8 <- numeric(num_replications)
  MSE9 <- numeric(num_replications)

  
  for (r in 1:num_replications) {
    # Generate data
    B <- rep(sqrt(1 / p), p)
    w <- rnorm(200 + n * p)  # Generate the 200 + n*p random numbers
    z <- matrix(w[201:(200 + n * p)], n, p)  # Disregard the first 200 random numbers in w 
    x <- ((1 - rho^2)^(1/2)) * z + (rho * z[, p])  # Generating p predictors with multicollinearity
    logits <- beta_0 + x %*% B
    probs <- 1 / (1 + exp(-logits))  # Logistic function to generate probabilities
    y <- rbinom(n, 1, probs)  # Binary response variable
    
    # Define the design matrix X and the response vector y
    X <- as.matrix(cbind(1, x))  # Adding intercept term
    
    # Compute lambda values
    # Compute lambda values
    xtx <- crossprod(x)
    Eigen <- eigen(xtx)
    Q<-Eigen$vectors
    Q<-t(Q)
    beta.ML <- logistic_ridge_regression(x,y,0)
    alpha.hat<-Q%*%beta.ML
    alpha.hat.max <- max(alpha.hat^2)
    mu.hat1 <- 1 / (1 + exp(-x %*% beta.ML))  # Predicted probabilities
    mu.hat1[mu.hat1 > 1e10] <- 1e10  # Limit mu.hat1 to prevent overflow
    
    
    sigma.hat.sq <- sum((y - mu.hat1)^2) / (n - p - 1)
    if (is.nan(sigma.hat.sq) || sigma.hat.sq == 0) {
      next
    }
    
    k0 <- 0
    k1 <- sigma.hat.sq / alpha.hat.max
    k2 <- 1 / alpha.hat.max
    k3 <- sigma.hat.sq / (prod(alpha.hat^2)^(1/p))
    k4 <- median(sqrt(sigma.hat.sq / alpha.hat^2))
    S <- rep(1, p) 
    for (i in 1:p) {
      S[i] <- Eigen$values[i] * sigma.hat.sq / ((n - p) * sigma.hat.sq + Eigen$values[i] * (alpha.hat[i])^2) 
    }
    
    k5 <- max(S)
    k6 <- max(1 / (sqrt(sigma.hat.sq / alpha.hat^2)))
    k7 <- (prod(1 / (sqrt(sigma.hat.sq / alpha.hat^2))))^(1/p)
    k8 <- median(1 / (sqrt(sigma.hat.sq / alpha.hat^2)))
    # After running the Monte Carlo simulation and obtaining MSE values
    MSE_values <- c(MSE0, MSE1, MSE2, MSE3, MSE4, MSE5, MSE6, MSE7, MSE8)  # Combine all MSE results
    
    # Calculate mean, median, and standard deviation of MSE
    mean_MSE <- mean(MSE_values, na.rm = TRUE)
    median_MSE <- median(MSE_values, na.rm = TRUE)
    std_MSE <- sd(MSE_values, na.rm = TRUE)
    
    # Small constant to ensure positive k
    c <- 0.01 
    
    # Compute k using the suggested formula
    k9 <- (mean_MSE + median_MSE) / (1 + std_MSE) + c
    
    
    # Compute coefficients using manual Logistic Ridge Regression
    beta_k0 <- logistic_ridge_regression(X, y, k0)
    beta_k1 <- logistic_ridge_regression(X, y, k1)
    beta_k2 <- logistic_ridge_regression(X, y, k2)
    beta_k3 <- logistic_ridge_regression(X, y, k3)
    beta_k4 <- logistic_ridge_regression(X, y, k4)
    beta_k5 <- logistic_ridge_regression(X, y, k5)
    beta_k6 <- logistic_ridge_regression(X, y, k6)
    beta_k7 <- logistic_ridge_regression(X, y, k7)
    beta_k8 <- logistic_ridge_regression(X, y, k8)
    beta_k9 <- logistic_ridge_regression(X, y, k9)
    
    # Adjusting BTrue for comparison
    BTrue <- c(beta_0, B)
    
    # Ensure BTrue matches the length of the coefficients
    if (length(BTrue) != length(beta_k0)) {
      BTrue <- c(BTrue, rep(0, length(beta_k0) - length(BTrue)))
    }
    
    # Compute MSE
    MSE0[r] <- sum((beta_k0 - BTrue)^2)
    MSE1[r] <- sum((beta_k1 - BTrue)^2)
    MSE2[r] <- sum((beta_k2 - BTrue)^2)
    MSE3[r] <- sum((beta_k3 - BTrue)^2)
    MSE4[r] <- sum((beta_k4 - BTrue)^2)
    MSE5[r] <- sum((beta_k5 - BTrue)^2)
    MSE6[r] <- sum((beta_k6 - BTrue)^2)
    MSE7[r] <- sum((beta_k7 - BTrue)^2)
    MSE8[r] <- sum((beta_k8 - BTrue)^2)
    MSE9[r] <- sum((beta_k9 - BTrue)^2)

  }
  
  # Summary statistics
  summary_stats <- data.frame(
    Lambda = c("k0", "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8","k9"),
    Min_MSE = c(min(MSE0, na.rm = TRUE), min(MSE1, na.rm = TRUE), min(MSE2, na.rm = TRUE), min(MSE3, na.rm = TRUE), min(MSE4, na.rm = TRUE), min(MSE5, na.rm = TRUE), min(MSE6, na.rm = TRUE), min(MSE7, na.rm = TRUE), min(MSE8, na.rm = TRUE),min(MSE9, na.rm = TRUE)),
    Max_MSE = c(max(MSE0, na.rm = TRUE), max(MSE1, na.rm = TRUE), max(MSE2, na.rm = TRUE), max(MSE3, na.rm = TRUE), max(MSE4, na.rm = TRUE), max(MSE5, na.rm = TRUE), max(MSE6, na.rm = TRUE), max(MSE7, na.rm = TRUE), max(MSE8, na.rm = TRUE),max(MSE9, na.rm = TRUE)),
    Mean_MSE = c(mean(MSE0, na.rm = TRUE), mean(MSE1, na.rm = TRUE), mean(MSE2, na.rm = TRUE), mean(MSE3, na.rm = TRUE), mean(MSE4, na.rm = TRUE), mean(MSE5, na.rm = TRUE), mean(MSE6, na.rm = TRUE), mean(MSE7, na.rm = TRUE), mean(MSE8, na.rm = TRUE),mean(MSE9, na.rm = TRUE)),
    Greater_Than_k0 = c(sum(MSE0 > MSE0), 
                        sum(MSE1 > MSE0), 
                        sum(MSE2 > MSE0),
                        sum(MSE3 > MSE0),
                        sum(MSE4 > MSE0),
                        sum(MSE5 > MSE0),
                        sum(MSE6 > MSE0),
                        sum(MSE7 > MSE0),
                        sum(MSE8 > MSE0),
                        sum(MSE9 > MSE0))
  )
  
  print(summary_stats)
}
