alpha_fun <- function(z, c){
  ifelse(z == 0, 1, psi_fun(z, c) * z^{-1})
}
