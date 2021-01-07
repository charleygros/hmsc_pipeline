evaluate_convergence <- function(modelDir, samples_list, thin_list, nChains, panelsDir) {
  nst = length(thin_list)
  
  ma.beta = ma.V = NULL
  na = NULL
  for (Lst in 1:nst) {
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    
    #filenameout = paste("MCMCconvergence_thin_", as.character(thin),
    #                    "_samples_", as.character(samples),
    #                    "_chains_",as.character(nChains),
    #                    ".pdf",sep = "")
    #if (file.exists(filenameout)) {
    #  print(paste("File already exists -- Skipping:", filenameout))
    #
    #} else {
    filename = paste("models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),".Rdata",sep = "")
    filename = file.path(modelDir, filename)
    print(" ")
    print(paste("Examining: ", filename))
    load(filename)
    nm = length(models)
    for(j in 1:nm){
      print(modelnames[[j]])
      # We first extract the posterior distribution from the model object and 
      # convert it into a coda object (needed to examine its convergence).
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      
      # Explore the MCMC convergence
      # Effective size of the posterior sample.
      # beta-parameters (species niches) 
      ess.beta = effectiveSize(mpost$Beta)
      hist(ess.beta, xlab = expression("Effective sample size" ~ beta ~ ""))
      title(main="Histogram of ess.beta",
            sub=paste("Actual sample size: ", samples * nChains))
      # V-parameters (variation in species niches)
      ess.V = effectiveSize(mpost$V)
      hist(ess.V, xlab = "Effective sample size V")
      title(main="Histogram of ess.V",
            sub=paste("Actual sample size: ", samples * nChains))
      
      # Gelman diagnostics, i.e. the Potential scale reduction factors
      psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
      tmp = summary(psrf.beta)
      print("Gelman Diagnostic Beta parameters")
      print(tmp)
      psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
      tmp = summary(psrf.V)
      print("Gelman Diagnostic V parameters")
      print(tmp)
      if(is.null(ma.beta)){
        ma.beta=psrf.beta[,1]
        ma.V=psrf.V[,1]
        na = paste0(as.character(thin),",",as.character(samples))
      } else {
        ma.beta = cbind(ma.beta,psrf.beta[,1])
        ma.V = cbind(ma.V,psrf.V[,1])
        if(j==1){
          na = c(na,paste0(as.character(thin),",",as.character(samples)))
        } else {
          na = c(na,"")
        }
      }
    }
    
    #pdf(file=file.path(panelsDir, filenameout))
    #par(mfrow=c(2,1))
    vioplot(ma.beta,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma.beta)),main="psrf(beta)")
    vioplot(ma.beta,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
    vioplot(ma.V,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma.V)),main="psrf(V)")
    vioplot(ma.V,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(V)")
    #dev.off()
    #}
  }
}


