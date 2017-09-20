myFaPXEM <- function (Y, nZ, tol = 1e-10, maxiter = 1000, estimationType = 1) 
{
  nVox <- ncol(Y)
  nObs <- nrow(Y)
  I <- diag(nZ)
  
  diagVarY <- colMeans(Y^2)
  
  convergence = FALSE
  
  # initialise BetaZ and sigma
  init <- cate::fa.pc(Y, nZ)
  # sigma = apply(Y,2, function(x) var(x))
  sigma <- init$Sigma
  tBetaZ <- init$Gamma
  invV <- I
  # sigma <- diagVarY
  # tBetaZ <- matrix(rnorm(nVox, nZ), nVox, nZ)
  eig <- eigen(I + crossprod((1/sqrt(sigma)) * tBetaZ), symmetric = TRUE)
  logDet <- sum(log(sigma)) + sum(log(eig$values))
  #   logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
  trSInvSigma <- sum(diagVarY / sigma) - sum((1/sqrt(eig$values) * t(eig$vectors) %*% t(Y %*% (tBetaZ / sigma )))^2) / nObs
  # cat( -logDet - trSInvSigma)

  sigmaRelDeltaMax = betaZRelDeltaMax = logLik <- numeric(maxiter)
  for (iter in 1:maxiter)
  {
    sigmaPrev <- sigma
    tBetaZPrev <- tBetaZ
    
    # E-step
    invSigma  <- 1 / sigma
    delta     <- invSigma * tBetaZ
    G         <- crossprod(Y, Y %*% delta) / nObs
    M2        <- invV + crossprod(tBetaZ, delta)
    tH        <- solve(M2, t(G))

    # M-step
    tmp <- tH %*% delta
    tmp <- solve(I +  t(tmp))
    
    invV <- tmp %*% M2

    tBetaZ <-  tcrossprod(G, tmp)
    sigma <- diagVarY - rowSums(tBetaZ * t(tH))
    
    # EZ  <- BetaZ %*% invMatrix %*% Y
    # EZZ <- I - tcrossprod(EZ, BetaZ)
    
    sigmaRelDeltaMax[iter] <- max(abs(((sigma-sigmaPrev)/sigmaPrev)))
    betaZRelDeltaMax[iter] <- max(abs(((tBetaZ-tBetaZPrev)/tBetaZPrev)))
    
    # log-likelihood
    # note: using Sylvester's determinant theorem
    cholV <- t(chol(solve(invV)))
    eig <- eigen(I + crossprod((1/sqrt(sigma)) * (tBetaZ %*% cholV)), symmetric = TRUE)
    logDet <- sum(log(sigma)) + sum(log(eig$values))
    #   logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
    trSInvSigma <- sum(diagVarY / sigma) - sum((1/sqrt(eig$values) * t(eig$vectors) %*% t(Y %*% ((tBetaZ %*% cholV) / sigma )))^2) / nObs
    logLik[iter] <- -logDet - trSInvSigma
    
    # stop if all the relative updates are less than the tolerance
    if (ifelse(iter > 1, logLik[iter] < logLik[iter-1], FALSE))
    {
      sigma <- sigmaPrev
      tBetaZ <- tBetaZPrev

      break
    }
    if (sigmaRelDeltaMax[iter] < tol & betaZRelDeltaMax[iter] < tol & ifelse(iter > 1, logLik[iter] - logLik[iter-1] < tol, FALSE))
    {
      convergence = TRUE
      break
    } 
  }
  tBetaZ <- tBetaZ %*% cholV
  invSigma  <- 1 / sigma
  delta <- invSigma * tBetaZ
  Z <-  Y %*% t(solve( I + crossprod(tBetaZ, delta), t(delta)))
  sigmaRelDeltaMax <- sigmaRelDeltaMax[1:iter]
  betaZRelDeltaMax <- betaZRelDeltaMax[1:iter]
  logLik <- logLik[1:iter] * (nObs/2)
  
  # if (restriction == TRUE)
  # {
  #   tBetaZ
  # }
  if (estimationType == 2)  # reestimate sigma based on the residuals 
  {
    sigma <- colMeans((Y - tcrossprod(Z, tBetaZ))^2)
  }
  if (estimationType == 3)  # reestimate betaZ and sigma based on the residuals 
  {
    # normalise first Z
    # tmp <- chol(crossprod(fit1$Z)/108)
    # Z <- Z %*% t(tmp)
    tBetaZ <- crossprod(Y,Z) / nObs
    sigma <- colMeans((Y - tcrossprod(Z, tBetaZ))^2)
  }
  return(list(Gamma = tBetaZ, Sigma = sigma, Z = Z, niter = iter,
              sigmaRelDeltaMax = sigmaRelDeltaMax, 
              betaZRelDeltaMax = betaZRelDeltaMax,
              logLik = logLik,
              convergence = convergence))
}