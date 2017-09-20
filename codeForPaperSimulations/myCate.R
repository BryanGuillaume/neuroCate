myCate <- function(Y, X, M, nZ = 0, faMeth = "PX-EM", estimationType = 1, Ynormalisation = FALSE, adjVar = TRUE, adjNullDistr = TRUE, boot = FALSE, typeBoot = "wbBeforeRotation", nB = 999, adjResBoot = TRUE, seedBoot = 0, resamplingMatrix = NULL) {
  
  nObs <- nrow(Y)
  nVox <- ncol(Y)
  nX   <- ncol(X)
  nM   <- ncol(M)
  XM   <- cbind(X,M)
  SigmaXM <- crossprod(XM)/nObs
  
  # check resampling matrix if supplied and needed
  if (!is.null(resamplingMatrix) & boot) {
    # is the number of columns correct?
    check = ncol(resamplingMatrix) == nB
    if(!check) stop("resamplingMatrix does not have the good number of columns")
    
    # is the number of rows correct?
    if (typeBoot == "wbBeforeRotation" | typeBoot == "resPermBeforeRotation"){
      check = check & (nrow(resamplingMatrix) == nObs)
    }
    if (typeBoot == "wbAfterRotation" | typeBoot == "resPermAfterRotation"){
      check = check & (nrow(resamplingMatrix) == (nObs - nX))
    }
    if(!check) stop("resamplingMatrix does not have the good number of rows")
    
    # is the numbers in the resampling matrix seems correct?
    if (typeBoot == "wbBeforeRotation" | typeBoot == "wbAfterRotation"){
      check = check & all(sort(unique(as.numeric(resamplingMatrix))) == c(-1,1))
      if(!check) stop("resamplingMatrix should contain -1 and 1")
    }
    if (typeBoot == "resPermBeforeRotation" ){
      check = check & all(apply(resamplingMatrix, 2, function(x) all(sort(unique(x)) == 1:nSubj)))
      if(!check) stop("resamplingMatrix should contain numbers from 1:nObs in each column")
    }    
    if (typeBoot == "resPermAfterRotation"){
      check = check & all(apply(resamplingMatrix, 2, function(x) all(sort(unique(x)) == 1:(nObs-nX))))
      if(!check) stop("resamplingMatrix should contain numbers from 1:(nObs-nX) in each column")
    }
    if(!check) stop("resamplingMatrix no well specified")
  }
  
  # run cate on the original data
  
  # compute the rotation matrix
  rotationMatrix <- t(qr.Q(qr(cbind(X,M)), complete = TRUE))
  rotatedMM <- rotationMatrix[(nX + 1): (nX + nM),] %*% M
  invRotatedMM <- solve(rotatedMM)

  # rotate data, keeping only the useful rotated data
  YR <- rotationMatrix[-(1:nX),] %*% Y
  rm(Y) # for memory
  
  #  factor analysis on the rotated data (removing the first row)
  if (nZ > 0){
    if (faMeth == "PX-EM")
    {
      faFit <- myFaPXEM(YR[-1,], nZ, estimationType = estimationType)
    } 
    
    if (faMeth == "ml")
    {
      
      faFit <- cate::factor.analysis(YR[-1,], r = nZ, method = faMeth)
    }

    if (faMeth =="pc")
    {
      faFit <- myFaPC(YR[-1,], nZ = nZ, Ynormalisation = Ynormalisation)
    }
 
    if (faMeth =="ppc")
    {
      faFit <- myFaPPC(YR[-1,], nZ = nZ)
    }
    
    
    # adjust variance estimates (cate does not account for the presence of missing factors)
    if (adjVar) faFit$Sigma <- faFit$Sigma * ((nObs - nX - nM) / (nObs - nX - nM - nZ))

    # TODO check the equation  below  should be sqrt(diag(invRotatedMM %*% t(invRotatedMM)))
    scaledSigma <- sqrt(diag(t(invRotatedMM) %*% invRotatedMM))
    rrFit <- cate::adjust.latent(corr.margin = t(invRotatedMM %*% YR[1:nM,]),
                                 n = nObs,
                                 X.cov = SigmaXM,
                                 Gamma = faFit$Gamma, 
                                 Sigma = faFit$Sigma,
                                 sigma.scales = scaledSigma,
                                 method = "rr")
    # compute t-score
    tScore <- (rrFit$beta / sqrt(rrFit$beta.cov.row)) * (sqrt(nObs) / sqrt(diag(rrFit$beta.cov.col)))
    
    # compute Z
    Z <- t(rotationMatrix) %*% rbind((rotationMatrix[1:(nX + nM),] %*% M)  %*% t(rrFit$alpha) , faFit$Z) 

  } else { # normal OLS regression
    faFit <- list()
    faFit$Sigma <- colMeans(YR[-1,]^2)
    tScore <- (YR[1,] / sqrt(faFit$Sigma)) * sign(rotatedMM)
    Z <- NULL
    rrFit <- list()
    rrFit$alpha <- NULL
  }
  #  if a resmapling method is used, save the first maxScore
  if (boot) {
    maxScore <- numeric(nB+1)
    maxScore[1] <- max(abs(tScore))
  }
  # compute p-values (2-tails t-test)
  if (adjNullDistr) {
    pParamUnc <- 2 * pt(abs(tScore), df = nObs - nX - nM - nZ, lower.tail = FALSE)
  } else {
    pParamUnc <- 2 * pnorm(abs(tScore),lower.tail = FALSE)
  }
  if (boot) {
    # set the random seed
    set.seed(seedBoot)
    # compute fitted data and adjusted residuals to resample
    if (nZ > 0){
      ZStar <- rbind((rotationMatrix[(nX+1):(nX + nM),] %*% M)  %*% t(rrFit$alpha) , faFit$Z) 
      fittedYR <- tcrossprod(ZStar, faFit$Gamma)
      adjResR <- YR - fittedYR
      if (adjResBoot) {
        adj <- 1 / sqrt(1 - hat(ZStar, intercept = FALSE))
        adjResR <- apply(adjResR, 2, function(x) adj * x) 
      }
    } else {
      fittedYR <- 0 * YR
      adjResR <- YR
    }
    rm(YR) # for memory
    # resample data and produce wild bootstrap scores
    pBootUnc = pBootFWE <- rep(1, nVox)
    cat("Bootstrap # 0")
    for (iB in 1:nB) {
      if(iB < 11) {
        cat("\b\b", iB)
      } else {
        if(iB < 101) {
          cat("\b\b\b", iB)
        } else {
          if(iB < 1001) {
            cat("\b\b\b\b", iB)
          } else { # assuming a number inZerior to 10000
            cat("\b\b\b\b\b", iB)
          }
        }
      }
      if (typeBoot == "wbAfterRotation"){ # not working well
        # generate realisations from Rademacher distribution
        if (is.null(resamplingMatrix)){
          f <- sample(c(-1,1), nObs - nX, replace = TRUE)
        } else {
          f <- resamplingMatrix[,iB]
        }
        YRB <- fittedYR  + apply(adjResR, 2, function(x) f * x)
      }
      if (typeBoot == "wbBeforeRotation"){ # not working well
        # generate realisations from Rademacher distribution
        if (is.null(resamplingMatrix)){
          f <- sample(c(-1,1), nObs, replace = TRUE)
        } else {
          f <- resamplingMatrix[,iB]
        }
        YRB <- fittedYR  + rotationMatrix[-(1:nX),] %*% apply(t(rotationMatrix[-(1:nX),]), 2, function(x) f * x) %*% adjResR
      }
      # if (typeBoot == "resBoot"){ # not working well
      #   f <- sample(1:(nObs - nX), nObs - nX, replace = TRUE)
      #   YRB <- fittedYR  + adjResR[f,]
      # }
      if (typeBoot == "resPermAfterRotation"){
        if (is.null(resamplingMatrix)){
          f <- sample(1:(nObs - nX), nObs - nX, replace = FALSE)
        } else {
          f <- resamplingMatrix[,iB]
        }
        YRB <- fittedYR  + adjResR[f,]
      } 
      if (typeBoot == "resPermBeforeRotation"){
        if (is.null(resamplingMatrix)){
          f <- sample(1:nObs, nObs, replace = FALSE)
        } else {
          f <- resamplingMatrix[,iB]
        }
        YRB <- fittedYR  + rotationMatrix[-(1:nX),] %*% t(rotationMatrix[-(1:nX),f]) %*% adjResR
      }       
      # factor analysis of this bootstrap rotated data (removing the first row)
      if (nZ > 0){
        if (faMeth == "PX-EM")
        {
          faFitboot <- myFaPXEM(YRB[-1,], nZ)
        } 
        if (faMeth == "ml")
        {
          faFitboot <- cate::factor.analysis(YRB[-1,], r = nZ, method = faMeth)
        }
        if (faMeth =="pc")
        {
          faFitboot <- myFaPC(YRB[-1,], nZ = nZ, Ynormalisation = Ynormalisation)
        }
        if (faMeth =="ppc")
        {
          faFitboot <- myFaPPC(YRB[-1,], nZ = nZ)
        }
        # adjust variance estimates (cate does not account for the presence of missing factors)
        if (adjVar) faFitboot$Sigma <- faFitboot$Sigma * ((nObs - nX - nM) / (nObs - nX - nM - nZ))
        
        # robust regression
        rrFitBoot <- cate::adjust.latent(corr.margin = t(invRotatedMM %*% YRB[1:nM,]),
                                         n = nObs,
                                         X.cov = SigmaXM,
                                         Gamma = faFitboot$Gamma, 
                                         Sigma = faFitboot$Sigma,
                                         sigma.scales = scaledSigma,
                                         method = "rr")
        
        # build the wild bootstrap score
        tScoreBoot <- (rrFitBoot$beta / sqrt(rrFitBoot$beta.cov.row)) * (sqrt(nObs) / sqrt(diag(rrFitBoot$beta.cov.col)))
        
      } else { # normal OLS regression
        faFitboot <-list()
        faFitboot$Sigma <- colMeans(YRB[-1,]^2)
        tScoreBoot <- (YRB[1,] / sqrt(faFitboot$Sigma)) * sign(rotatedMM)
      }
      maxScore[iB+1] <- max(abs(tScoreBoot))
      pBootUnc <- pBootUnc + 1 * (abs(tScore) <= abs(tScoreBoot))
      pBootFWE <- pBootFWE + 1 * (abs(tScore) <= maxScore[iB+1])
    } # end for (iB in 1:nB)
    pBootUnc <- pBootUnc / (nB + 1)
    pBootFWE <- pBootFWE / (nB + 1)
  } else {
    pBootUnc = pBootFWE = maxScore <- NULL
  }
  cat("\n")
  return(list(sigma = faFit$Sigma, tScore = tScore, Z = Z, betaZ = faFit$Gamma, alpha = rrFit$alpha, pParamUnc = pParamUnc, pBootUnc = pBootUnc, pBootFWE = pBootFWE, maxScore = maxScore))
}
