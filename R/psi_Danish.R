psi_Danish <- function(u, c = 2.767) {
  ifelse(abs(u) <= c, u, u * exp(-u^2 / c))
}
