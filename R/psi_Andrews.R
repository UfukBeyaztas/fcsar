psi_Andrews <- function(u, c = 1.339) {
  ifelse(abs(u) <= c * pi, c * sin(u / c), 0)
}
