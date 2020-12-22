###
# Correlation in the beta coefficients between species.
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
  sig_hmsc <- (support1 > prob | support1 < (1-prob))*mBetaCor
  
  ## calculate support/significance- boral way
  #cor_int <- apply(abind::abind(BetaCor, along = 3), c(1, 2), function(a) HPDinterval(as.mcmc(a), prob = prob))
  #cor_int <- aperm(cor_int, c(2, 3, 1))
  ## for boral, the number of times that the interval does not include 0
  #sig_boral <- (apply(cor_int, c(1, 2), function(z) (z[1] < 0 & z[2] < 0) | (z[1] > 0 & z[2] > 0)))*mBetaCor
  
  colnames(mBetaCor) <- colnames(sig_hmsc)  <- hM$spNames ## <-colnames(mBetaCov) <-colnames(sig_boral) <- hM$spNames
  rownames(mBetaCor) <- rownames(sig_hmsc) <- hM$spNames ## <-rownames(mBetaCov) <-rownames(sig_boral) <- hM$spNames
  
  list(cor = mBetaCor, signifiance = sig_hmsc)  #cov = mBetaCov, sig_boral=sig_boral)
}
