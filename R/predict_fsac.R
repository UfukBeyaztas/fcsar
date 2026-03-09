predict_fsac <- function(object, xnew, wnew){

  n <- dim(xnew)[1]
  rho <- object$rho
  b0 <- object$b0
  b <- object$b
  fpca <- object$fpca
  method.type <- fpca$method.type
  fsco.test <- getPCA_test(fpca, xnew)
  
  predicted.values <- ginv(diag(n) - rho*wnew) %*% as.matrix(rep(1,n)) * b0 +
    ginv(diag(n) - rho*wnew) %*% fsco.test %*% b
  
  return(predicted.values)
}
