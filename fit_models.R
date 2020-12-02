fit_models <- function(modelDir, samples_list, thin_list, nChains) {
  load(file = file.path(modelDir, "unfitted_models"))
  
  for(Lst in 1:length(samples_list)){
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    nm = length(models)
    for (model in 1:nm) {
      print(paste0("model = ",modelnames[model]))
      m = models[[model]]
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains) 
      models[[model]] = m
    }
    filename = paste("models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")
    filename = file.path(modelDir, filename)
    save(models,modelnames,file=filename)
  }
}