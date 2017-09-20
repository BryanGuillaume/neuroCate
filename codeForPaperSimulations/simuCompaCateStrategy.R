# Simulations for paper
# comparison between PCA & PX-EM (2)
# comparison between param, permAfter, permBefore, wbAfter, wbBefore (5)
# Adj. residual in the resampling?
# comparison between C = I, toeplitz_50, toeplitz_250, I_toeplitz_50, I_toeplitz_250 (5)
# comparison between nVox = 100 & 500 (2)
# comparison between nZ = 0, 2 (2)
# comparison between different noise type
# comparison between alpha = 0 0.1 0.3 (3)
Args = commandArgs(TRUE)
print(Args)
sessionInfo()

iSim <- as.numeric(Args[1])
nZ <- as.numeric(Args[2])
alpha <- as.numeric(Args[3])

nSubj = 114
index100Vox <- c(1:50, 251:300)

source('/home/users/nus/biebrlg/scripts/myCate.R')
source('/home/users/nus/biebrlg/scripts/myFaPC.R')
source('/home/users/nus/biebrlg/scripts/myFaPXEM.R')

load("~/data/covariates114Subj42173CpgNew.Rdata")
dfCovariates <- data.frame(int = rep(1, nSubj), 
                           ga = adjGA114Subj42173CpgNewCentered, 
                           female = gender114Subj42173CpG,
                           malay = as.numeric(ethnicity114Subj42173CpG == 2),
                           indian = as.numeric(ethnicity114Subj42173CpG == 3)
)
X <- as.matrix(dfCovariates)

# generate random variables (but the same accross simulations)
set.seed(0)
trueBetaZ <- matrix(rnorm(2*500), 2, 500)
trueBetaX <- matrix(rnorm(5*500), 5, 500)

tmp <- matrix(rnorm(nSubj*3), nSubj, 3)
tmp2 <- eigen(cov(tmp))
tmp <- tmp %*% (tmp2$vec %*% sqrt(diag(1/tmp2$val)) %*% t(tmp2$vec))
tmp <-apply(tmp, 2, function(x) x - mean(x))
Z <- tmp[,1:2]
M <- matrix(tmp[,3], nSubj, 1)

var <- rexp(500)
var[index100Vox] <- var[index100Vox] / mean(var[index100Vox])
var[-index100Vox] <- var[-index100Vox] / mean(var[-index100Vox])

effectX <- X %*% trueBetaX

if (iSim == 1 & nZ == 2 & alpha == 0)
{
  save(var, trueBetaZ, trueBetaX, Z,  M, file = "metaDataSimulationCompaCateStrategy.Rdata")
}


pUnc500 <-  array(NA, dim = c(16, 2, 5, 2, 500))
dimnames(pUnc500)[[1]] <- c("noCorr_Hom", "lowCorr_hom", "highCorr_hom", "lowHighCorr_hom", 
                            "expDistr_hom", "unifDistr_hom", "hetNorm1Distr_hom", "hetNorm2Distr_hom",
                            "noCorr_Het", "lowCorr_het", "highCorr_het", "lowHighCorr_het", 
                            "expDistr_het", "unifDistr_het", "hetNorm1Distr_het", "hetNorm2Distr_het")
dimnames(pUnc500)[[2]] <- c("PX-EM", "pc")
dimnames(pUnc500)[[3]] <- c("param", "wbBeforeRotation", "wbAfterRotation", "resPermBeforeRotation", "resPermAfterRotation")
dimnames(pUnc500)[[4]] <- c("adjBootRes", "nonAdjBootRes")

sigma500 <- pUnc500

sigma100 = pUnc100 <- pUnc500[,,,,1:100]
pFwe100 = pFwe500 <- pUnc500[,,,,1]

# generate the noises
set.seed(iSim)

noise <- array(0, dim = c(16, nSubj, 500))
dimnames(noise)[[1]] <- dimnames(pUnc500)[[1]]
noise["noCorr_Hom",,]  <- mvtnorm::rmvnorm(nSubj, sigma = diag(rep(1,500)))
noise["lowCorr_hom",,] <- mvtnorm::rmvnorm(nSubj, sigma = toeplitz(c(seq(1, 0.02, -0.02), rep(0, 500-50))))
noise["highCorr_hom",,] <- mvtnorm::rmvnorm(nSubj, sigma = toeplitz(c(seq(1, 0.004, -0.004), rep(0, 500-250))))
noise["lowHighCorr_hom",,] <- cbind(noise["lowCorr_hom",, 1:250], noise["highCorr_hom",,251:500])
noise["expDistr_hom",,] <- matrix(rexp(nSubj*500, rate = 1) - 1, nSubj, 500)
noise["unifDistr_hom",,] <- matrix(runif(nSubj * 500, -sqrt(3), sqrt(3)), nSubj, 500) 
noise["hetNorm1Distr_hom",,] <- matrix(rnorm(nSubj * 500, sd = sqrt(t(kronecker(matrix(rep(c(0.4, 1.6), nSubj/2), 1, nSubj), rep(1,500))))), nSubj, 500)
noise["hetNorm2Distr_hom",,] <- matrix(rnorm(nSubj * 500, sd = sqrt(t(kronecker(matrix(rep(c(1.6, 0.4), nSubj/2), 1, nSubj), rep(1,500))))), nSubj, 500)

noise["noCorr_Het",,]  <- mvtnorm::rmvnorm(nSubj, sigma = diag(var))
noise["lowCorr_het",,] <- mvtnorm::rmvnorm(nSubj, sigma = diag(sqrt(var)) %*% toeplitz(c(seq(1, 0.02, -0.02), rep(0, 500-50))) %*% diag(sqrt(var)))
noise["highCorr_het",,] <- mvtnorm::rmvnorm(nSubj, sigma = diag(sqrt(var)) %*% toeplitz(c(seq(1, 0.004, -0.004), rep(0, 500-250))) %*% diag(sqrt(var)))
noise["lowHighCorr_het",,] <- cbind(noise["lowCorr_het",, 1:250], noise["highCorr_het",, 251:500])
noise["expDistr_het",,] <- matrix(rexp(nSubj*500, rate = 1) - 1, nSubj, 500) * sqrt(t(kronecker(matrix(1,1,nSubj),var)))
noise["unifDistr_het",,] <- matrix(runif(nSubj * 500, -sqrt(3), sqrt(3)), nSubj, 500) * sqrt(t(kronecker(matrix(1,1,nSubj), var)))
noise["hetNorm1Distr_het",,] <- matrix(rnorm(nSubj * 500, sd = sqrt(t(kronecker(matrix(rep(c(0.4, 1.6), nSubj/2), 1, nSubj), var)))), nSubj, 500)
noise["hetNorm2Distr_het",,] <- matrix(rnorm(nSubj * 500, sd = sqrt(t(kronecker(matrix(rep(c(1.6, 0.4), nSubj/2), 1, nSubj), var)))), nSubj, 500)
  
a <- Sys.time()
for (iNoise in 1:16)
{
  cat(dimnames(noise)[[1]][iNoise], ":\n", sep="")
  Y <- effectX + (Z  + alpha * cbind(M, M)) %*% trueBetaZ + noise[iNoise,,]
  for (iFaMeth in c("PX-EM", "pc"))
  {
    
    for (iTypeBoot in c("wbBeforeRotation", "wbAfterRotation", "resPermBeforeRotation", "resPermAfterRotation"))
    {
      for (iAdjBootRes in 1:2)
      {
        fit <- myCate(Y, X, M,
                      nZ = nZ, 
                      faMeth = iFaMeth, 
                      boot = TRUE,
                      typeBoot = iTypeBoot, 
                      nB = 999, 
                      adjResBoot = c(TRUE,FALSE)[iAdjBootRes])
        
        if (iTypeBoot == "wbBeforeRotation")  # only needed once
        {
          pUnc500[iNoise,  iFaMeth, "param", 1, ] = 
            pUnc500[iNoise,  iFaMeth, "param", 2, ] <- fit$pParamUnc
          
          pFwe500[iNoise,  iFaMeth, "param", ] <- min(min(fit$pParamUnc) * 500, 1)
        }
        
        pUnc500[iNoise,  iFaMeth, iTypeBoot, iAdjBootRes, ] <- fit$pBootUnc
        pFwe500[iNoise,  iFaMeth, iTypeBoot, iAdjBootRes] <- min(fit$pBootFWE)
        sigma500[iNoise,  iFaMeth, iTypeBoot, iAdjBootRes, ] <- fit$sigma
        fit <- myCate(Y[, index100Vox], X, M,
                      nZ = nZ, 
                      faMeth = iFaMeth, 
                      boot = TRUE,
                      typeBoot = iTypeBoot, 
                      nB = 999, 
                      adjResBoot = c(TRUE,FALSE)[iAdjBootRes])
        
        if (iTypeBoot == "wbBeforeRotation")  # only needed once
        {
          pUnc100[iNoise, iFaMeth, "param", 1, ] = 
            pUnc100[iNoise, iFaMeth, "param", 2, ] <- fit$pParamUnc
          
          pFwe100[iNoise, iFaMeth, "param", ] <- min(min(fit$pParamUnc) * 100, 1)
        }
        
        pUnc100[iNoise, iFaMeth, iTypeBoot, iAdjBootRes, ] <- fit$pBootUnc
        pFwe100[iNoise, iFaMeth, iTypeBoot, iAdjBootRes] <- min(fit$pBootFWE)
        sigma100[iNoise, iFaMeth, iTypeBoot, iAdjBootRes, ] <- fit$sigma           
        
      }
    }
  }
}
cat(Sys.time() - a)

if(alpha == 0) alphaPrint = "000"
if(alpha < 0.1 & alpha > 0) alphaPrint = paste("00", alpha*100, sep ="")
if(alpha >= 0.1 & alpha <= 1) alphaPrint = paste("0", alpha*100, sep ="")

eval(parse(text=paste("save(pUnc500, pUnc100, pFwe500, pFwe100, sigma500, sigma100, file = \"/scratch/users/nus/biebrlg/simuCompaCateStrategy/simuCompaCateStrategy_",
                      nZ, "Z_", alphaPrint, "alpha_",iSim ,".Rdata\")",sep="")))

