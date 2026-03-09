dbias_dv <- function(W, lambda, sigma, kappa) {
  n <- nrow(W)
  I <- diag(n)
  S <- I - lambda * W
  S_inv <- solve(S)
  dS_inv <- S_inv %*% W %*% S_inv
  dbias_i <- diag(W %*% dS_inv) * sigma * kappa
  return(sum(dbias_i))
}
