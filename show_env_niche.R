###
# Correlation in the beta coefficients between species.
# i.e. correlations between species due to the environmental response
# Shows species with similar/dissimilar environmental niches
###

env_niches <- function (hM, start = 1, thin = 1, prob = 0.95) {
  # Combines a list of single or several MCMC chains into a single chain
  postList <- poolMcmcChains(hM$postList, start = start, thin = thin)
  
  # Extract betas, generate linear predictor and 
  X <- hM$X
  get.enviro.linpreds <- function(a) X %*% a$Beta
  enviro.linpreds <- lapply(postList,  get.enviro.linpreds)
  # Correlation and covariance matrices for each mcmc sample
  BetaCov <- lapply(enviro.linpreds, function(x) cov(x))
  BetaCor <- lapply(enviro.linpreds, function(x) cor(x))
  
  # Calculate means across MCMC samples
  mBetaCor <- apply(abind::abind(BetaCor, along = 3), c(1, 2), mean)
  mBetaCov <- apply(abind::abind(BetaCov, along = 3), c(1, 2), mean)
  
  # Calculate support/significance - HMSC way
  BetaCor2 <- lapply(BetaCor, function(a) return(a > 0))
  support1 <- apply(abind::abind(BetaCor2, along = 3), c(1, 2), mean)
  # Count the proportion of times the cor/cov is above or below 0
  signi <- (support1 > prob | support1 < (1-prob))*mBetaCor
  
  colnames(mBetaCor) <- colnames(signi)  <- hM$spNames
  rownames(mBetaCor) <- rownames(signi) <- hM$spNames
  
  list(cor = mBetaCor, signifiance = signi)
}
