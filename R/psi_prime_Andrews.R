psi_prime_Andrews <- function(u, c = 1.339) {
  ifelse(abs(u) <= c * pi, cos(u / c), 0)
}
