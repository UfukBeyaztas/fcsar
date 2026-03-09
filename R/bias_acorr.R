bias_acorr <- function(W, X, lambda, sigma, beta) {
  
  n <- nrow(W)
  I_n <- diag(n)
  
  tr <- function(m) sum(diag(m))
  M <- I_n - X %*% solve(crossprod(X)) %*% t(X)
  
  S <- I_n - lambda * W
  S_inv <- solve(S)
  
  A <- t(W) %*% M %*% S
  B <- t(W) %*% M %*% W
  
  mu_y <- S_inv %*% X %*% beta
  Sigma_yy <- (sigma^2) * S_inv %*% t(S_inv)
  
  A_Sigma <- A %*% Sigma_yy
  B_Sigma <- B %*% Sigma_yy
  
  E_R <- tr(A_Sigma) + drop(t(mu_y) %*% A %*% mu_y)
  E_S <- tr(B_Sigma) + drop(t(mu_y) %*% B %*% mu_y)
  
  Cov_RS <- 2 * tr(A_Sigma %*% B_Sigma) + 4 * drop(t(mu_y) %*% A %*% Sigma_yy %*% B %*% mu_y)
  Var_S  <- 2 * tr(B_Sigma %*% B_Sigma) + 4 * drop(t(mu_y) %*% B %*% Sigma_yy %*% B %*% mu_y)
  
  bias_lambda <- (E_R / E_S) - (Cov_RS / (E_S^2)) + (Var_S * E_R) / (E_S^3)
  
  return(bias_lambda)
}
