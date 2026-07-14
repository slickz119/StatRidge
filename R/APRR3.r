library(MASS)


generate_X <- function(n, p, rho){
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  X_raw <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- cbind(1, X_raw)
  return(X)
}


generate_Y <- function(X, beta, phi, contamination = 0){
  
  eta <- X %*% beta
  mu <- exp(eta)
  
  # Gamma noise for overdispersion
  v <- rgamma(length(mu), shape = 1/phi, scale = phi)
  mu_star <- mu * v
  
  y <- rpois(length(mu), mu_star)
  
  # Contamination (extreme counts)
  if(contamination > 0){
    n_cont <- floor(contamination * length(y))
    idx <- sample(1:length(y), n_cont)
    y[idx] <- y[idx] + rpois(n_cont, lambda = 10 * mean(mu))
  }
  
  return(y)
}


irls_poisson <- function(X, y, k = 0, weights_extra = NULL, max_iter = 50){
  
  X <- as.matrix(X)
  y <- as.vector(y)
  
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  
  for(iter in 1:max_iter){
    
    eta <- X %*% beta
    mu <- exp(eta)
    mu[mu > 1e6] <- 1e6
    mu[mu < 1e-6] <- 1e-6
    
    
    phi <- sum((y - mu)^2 / mu) / (n - p)
    
    # Ensure numeric vector (fixes your error)
    w_vec <- as.numeric(mu / phi)
    
    if(!is.null(weights_extra)){
      weights_extra <- as.numeric(weights_extra)
      
      if(length(weights_extra) != length(w_vec)){
        stop("weights_extra length mismatch")
      }
      
      w_vec <- w_vec * weights_extra
    }
    
    # Safe diagonal matrix
    W <- diag(w_vec)
    
    z <- eta + (y - mu)/(mu + 1e-6)
    
    A <- t(X) %*% W %*% X + k * diag(p)
    
    # Stabilization
    A <- A + 1e-6 * diag(p)
    
    beta_new <- tryCatch(
      solve(A, t(X) %*% W %*% z),
      error = function(e) MASS::ginv(A) %*% (t(X) %*% W %*% z)
    )
    
    if(max(abs(beta_new - beta)) < 1e-6) break
    beta <- beta_new
  }
  
  return(list(beta = as.vector(beta), mu = mu, phi = phi))
}

compute_robust_weights <- function(y, mu, phi, c = 1.345){
  
  r <- (y - mu) / sqrt(phi * mu + 1e-6)
  
  psi <- ifelse(abs(r) <= c, r, c * sign(r))
  omega <- psi / (r + 1e-6)
  
  return(omega)
}


compute_k_candidates <- function(X, y, beta_init){
  
  XtX <- t(X) %*% X
  eigvals <- eigen(XtX)$values
  
  sigma2 <- var(y)
  
  k1 <- sigma2 / max(eigvals)
  k2 <- 1 / max(eigvals)
  k3 <- sigma2 / (prod(eigvals)^(1/length(eigvals)))
  k4 <- median(sigma2 / eigvals)
  k5 <- max(eigvals / sum(eigvals))
  k6 <- max(1 / eigvals)
  k7 <- (prod(1 / eigvals))^(1/length(eigvals))
  k8 <- median(1 / eigvals)
  
  return(c(k1, k2, k3, k4, k5, k6, k7, k8))
}


compute_adaptive_k <- function(X, y, k_candidates){
  
  mse_vals <- numeric(length(k_candidates))
  
  for(i in 1:length(k_candidates)){
    fit <- irls_poisson(X, y, k = k_candidates[i])
    beta_hat <- fit$beta
    
    # Approximate MSE via fitted residuals
    mu_hat <- exp(X %*% beta_hat)
    mse_vals[i] <- mean((y - mu_hat)^2)
  }
  
  weights <- 1 / (mse_vals + 1e-6)
  k_adapt <- sum(k_candidates * weights) / sum(weights)
  
  return(k_adapt)
}

fit_ARPRR <- function(X, y){
  
  # Step 1: Initial fit
  init <- irls_poisson(X, y, k = 0)
  
  # Step 2: Robust weights
  omega <- compute_robust_weights(y, init$mu, init$phi)
  
  # Step 3: Candidate k
  k_candidates <- compute_k_candidates(X, y, init$beta)
  
  # Step 4: Adaptive k
  k_adapt <- compute_adaptive_k(X, y, k_candidates)
  k_adapt <- max(k_adapt, 1e-4)
  # Step 5: Final robust ridge fit
  final <- irls_poisson(X, y, k = k_adapt, weights_extra = omega)
  
  return(list(beta = final$beta, k = k_adapt))
}


compute_MSE <- function(beta_hat, beta_true){
  return(sum((beta_hat - beta_true)^2))
}


set.seed(123)

n <- 200
p <- 3
rho <- 0.99
phi <- 5
contam <- 0.20
R <- 200

beta_true <- c(0.5, 1, -1, 0.5)

# ✅ UPDATED results storage
results <- data.frame(
  PRR = numeric(R),
  ARPRR = numeric(R),
  k1 = numeric(R),
  k2 = numeric(R),
  k3 = numeric(R),
  k4 = numeric(R),
  k5 = numeric(R),
  k6 = numeric(R),
  k7 = numeric(R),
  k8 = numeric(R)
)

for(r in 1:R){
  
  # -------------------------------
  # STEP 1: Generate data
  # -------------------------------
  X <- generate_X(n, p, rho)
  X[, -1] <- scale(X[, -1])   # ✅ scaling HERE
  y <- generate_Y(X, beta_true, phi, contam)
  
  # -------------------------------
  # STEP 2: Compute candidate k’s
  # -------------------------------
  k_vals <- compute_k_candidates(X, y, beta_true)
  
  # -------------------------------
  # STEP 3: Evaluate k₁–k₈ estimators
  # -------------------------------
  for(i in 1:8){
    fit_k <- irls_poisson(X, y, k = k_vals[i])
    results[r, paste0("k", i)] <- compute_MSE(fit_k$beta, beta_true)
  }
  
  # -------------------------------
  # STEP 4: PRR (baseline)
  # -------------------------------
  prr <- irls_poisson(X, y, k = 0.0)
  results$PRR[r] <- compute_MSE(prr$beta, beta_true)
  
  # -------------------------------
  # STEP 5: ARPRR (your method)
  # -------------------------------
  arprr <- fit_ARPRR(X, y)
  results$ARPRR[r] <- compute_MSE(arprr$beta, beta_true)
}

# Summary
results1=colMeans(results)
results2=apply(results, 2, sd)
results3=colMeans(results) / colMeans(results)["ARPRR"]
results4= mean(results$ARPRR < results$PRR)

write.csv(results1, "results106_1.csv")
write.csv(results2, "results106_2.csv")
write.csv(results3, "results106_3.csv")
write.csv(results4, "results106_4.csv")

