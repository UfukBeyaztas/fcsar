weight_fn <- function(u, Type, c_val) {
  if (Type == "Andrews") {
    psi_u <- psi_Andrews(u, c_val)
  } else if (Type == "Danish") {
    psi_u <- psi_Danish(u, c_val)
  } else {
    stop("Unknown Type")
  }
  w <- ifelse(u == 0, 1, psi_u / u)
  w[abs(u) > 1e8] <- 0
  return(w)
}
