psi_prime_Danish <- function(u, c = 2.767) {
  ifelse(abs(u) <= c, 1, exp(-u^2 / c) * (1 - 2 * u^2 / c))
}
