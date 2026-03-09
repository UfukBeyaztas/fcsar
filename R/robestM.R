robestM <- function(y, X, W, maxIter = 100, tol = 1e-6, Type = "Andrews") {
  
  y <- as.numeric(y)
  X <- as.matrix(X)
  W <- as.matrix(W)
  
  n <- length(y)
  k <- ncol(X)
  I <- diag(n)
  
  kappa <- switch(Type,
                  "Andrews" = 0.829,
                  "Danish"  = 0.872)
  c_val <- switch(Type,
                  "Andrews" = 1.339,
                  "Danish"  = 2.767)
  
  psi <- switch(Type,
                "Andrews" = function(u) psi_Andrews(u, c_val),
                "Danish"  = function(u) psi_Danish(u, c_val))
  psi_prime <- switch(Type,
                      "Andrews" = function(u) psi_prime_Andrews(u, c_val),
                      "Danish"  = function(u) psi_prime_Danish(u, c_val))
  
  Wy <- W %*% y
  A <- cbind(Wy, X)
  fit_ols <- lm(y ~ A - 1)
  theta_ols <- coef(fit_ols)
  if (any(is.na(theta_ols))) {
    fit_ols <- lm(y ~ X - 1)
    theta_ols <- c(0, coef(fit_ols))
  }
  lambda <- theta_ols[1]
  gamma <- theta_ols[-1]
  e <- y - A %*% theta_ols
  sigma <- median(abs(e)) / 0.6745
  
  bias_lambda <- bias_acorr(W, X, lambda, sigma, gamma)
  lambda <- lambda - bias_lambda

  lambda <- max(min(lambda, 0.99), -0.99)
  e <- y - lambda * Wy - X %*% gamma
  sigma <- median(abs(e)) / 0.6745
  
  lambda_old <- lambda
  gamma_old <- gamma
  sigma_old <- sigma
  
  for (iter in 1:maxIter) {
    u <- e / sigma
    w <- weight_fn(u, Type, c_val)
    w <- as.vector(w)
    w[!is.finite(w)] <- 0
    Wmat <- diag(w)
    
    y_tilde <- y - lambda * Wy
    XtWX <- crossprod(X, Wmat %*% X)
    XtWy <- crossprod(X, Wmat %*% y_tilde)
    gamma_new <- tryCatch(
      solve(XtWX, XtWy),
      error = function(e) gamma
    )
    
    # ---- Inner Newton for lambda ----
    for (inner in 1:10) {
      e <- y - lambda * Wy - X %*% gamma_new
      u <- e / sigma
      psi_u <- psi(u)
      bias_vec <- bias_vector(W, lambda, sigma, kappa)
      S_lambda <- sum(psi_u * Wy) - sum(bias_vec)
      
      psi_prime_u <- psi_prime(u)
      dS_lambda <- sum(psi_prime_u * (-Wy / sigma) * Wy) - dbias_dv(W, lambda, sigma, kappa)
      
      step <- 0.5
      lambda_new <- lambda - step * S_lambda / dS_lambda
      lambda_new <- max(min(lambda_new, 0.99), -0.99)
      
      if (abs(lambda_new - lambda) < 1e-8) break
      lambda <- lambda_new
    }
    lambda <- lambda_new
    
    # ---- Update sigma ----
    e <- y - lambda * Wy - X %*% gamma_new
    u <- e / sigma
    w_new <- weight_fn(u, Type, c_val)
    w_new <- as.vector(w_new)
    w_new[!is.finite(w_new)] <- 0
    sigma_new <- sqrt(max(sum(w_new * e^2) / (n * kappa), 1e-6))
    
    # ---- Convergence check ----
    if (max(abs(gamma_new - gamma_old),
            abs(lambda - lambda_old),
            abs(sigma_new - sigma_old)) < tol) {
      gamma <- gamma_new
      sigma <- sigma_new
      break
    }
    
    # Prepare next iteration
    gamma <- gamma_new
    sigma <- sigma_new
    e <- y - lambda * Wy - X %*% gamma
    u <- e / sigma
    w <- weight_fn(u, Type, c_val)
    lambda_old <- lambda
    gamma_old <- gamma
    sigma_old <- sigma
  }
  
  return(list(
    lambda_hat = lambda,
    beta_hat = gamma,
    sigma = sigma,
    pweights = Wmat,
    iterations = iter
  ))
}
