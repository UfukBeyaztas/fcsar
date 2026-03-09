bias_vector <- function(W, lambda, sigma, kappa) {
  n <- nrow(W)
  S_inv <- solve(diag(n) - lambda * W)
  bias_i <- diag(W %*% S_inv) * sigma * kappa
  return(bias_i)
}
