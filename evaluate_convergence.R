evaluate_convergence <- function(modelDir, samples_list, thin_list, nChains, panelsDir) {
  nst = length(thin_list)
  
  ma = NULL
  na = NULL
  for (Lst in 1:nst) {
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    
    
    filename = paste("models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),".Rdata",sep = "")
    filename = file.path(modelDir, filename)
    load(filename)
    nm = length(models)
    for(j in 1:nm){
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
      tmp = summary(psrf.beta)
      if(is.null(ma)){
        ma=psrf.beta[,1]
        na = paste0(as.character(thin),",",as.character(samples))
      } else {
        ma = cbind(ma,psrf.beta[,1])
        if(j==1){
          na = c(na,paste0(as.character(thin),",",as.character(samples)))
        } else {
          na = c(na,"")
        }
      }
    }
  }
  
  pdf(file=file.path(panelsDir, "MCMC_convergence.pdf"))
  par(mfrow=c(2,1))
  vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
  vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
  dev.off()
}


