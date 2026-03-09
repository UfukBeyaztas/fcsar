psi_fun <- function(z, c){
  z * pmin(1, c / abs(z))
}
