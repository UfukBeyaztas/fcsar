fsac_pcaM <- function(y, x, nbasis = NULL, gpx = NULL, wei_mat,
                      Type = c('Andrews', 'Danish'),
                      maxIter = 100, tol = 1e-6){

  Type <- match.arg(Type)
  n <- length(c(y))
  px <- dim(x)[2]
  
  if(is.null(gpx))
    gpx <- seq(1/px, 1-1/px, len = px)
  
  if(is.null(nbasis))
    nbasis <- round(min(n/4, 40))
  
  fpca <- getPCA(x, nbasis, gpx, method.type = "robust")
  fcomp <- fpca$PCAcoef
  fsco <- cbind(1, fpca$PCAscore)
  evalbase <- fpca$evalbase
  
  sac_model <- robestM(y=y, X=fsco, W=wei_mat, maxIter=maxIter, tol=tol, Type=Type)
  
  b0 <- sac_model$beta_hat[1]
  b <- sac_model$beta_hat[-1]
  rho <- max(min(sac_model$lambda_hat, 1), -1)
  sig <- sac_model$sigma
  
  fits <- ginv(diag(n) - rho*wei_mat) %*% as.matrix(rep(1,n)) * b0 + 
    ginv(diag(n) - rho*wei_mat) %*% as.matrix(fsco[,-1]) %*% b
  
  b_hat_t <- evalbase %*% (fcomp$coefs %*% b)
  resids <- y - fits
  
  return(list(b = b, b0 = b0,
              bhat = b_hat_t,
              rho = rho,
              sig = sig,
              fitted.values = fits,
              residuals = resids,
              fpca = fpca,
              ncomp = fpca$ncomp))
}
