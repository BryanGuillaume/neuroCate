myFaPC <- function (Y, nZ, Ynormalisation = FALSE) {

  if (Ynormalisation)
  {
    nObs <- nrow(Y)
    diagVarY <- colMeans(Y^2)
    svdY <- svd(t(t(Y) / sqrt(diagVarY)), nu = nZ, nv = nZ) 
    
    Gamma <- svdY$v %*% diag(svdY$d[1:nZ], nZ, nZ) * (sqrt(diagVarY / nObs))
    Z     <- sqrt(nObs) * svdY$u
    Sigma <- colMeans((Y - tcrossprod(Z, Gamma))^2)
  }
  else 
  {
    nObs <- nrow(Y)
    svdY <- svd(Y, nu = nZ, nv = nZ) 
    
    Gamma <- svdY$v %*% (diag(svdY$d[1:nZ], nZ, nZ) / sqrt(nObs))
    Z     <- sqrt(nObs) * svdY$u
    Sigma <- colMeans((Y - tcrossprod(Z, Gamma))^2)
  }
  return(list(Gamma = Gamma, Z = Z, Sigma = Sigma))
}