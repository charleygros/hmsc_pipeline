fit_models <- function(modelDir, fname_unfitted, samples, thin, nChains) {
  # load unfitted model
  load(file = fname_unfitted)
  
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  # Output fname
  prefix = strsplit(fname_unfitted, "models")[[1]][2]
  if (is.na(prefix)) {
    prefix = ""
  }
  filename = paste("models", prefix,
                   "_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  filename = file.path(modelDir, filename)
  
  # Fit model
  if (file.exists(filename)) {
    print(paste("File already exists -- Skipping:", filename))
  } else {
    print(paste("Fitting:", filename))
    nm = length(models)
    for (model in 1:nm) {
      print(paste0("model = ",modelnames[model]))
      m = models[[model]]
      # Note re nParallel:  we use one  CPU for each of the nChains.
        # Setting nParallel > 1 has two consequences: tracing info
        # vanishes, and random sequences will change from nParallel=1 and
        # results are no longer reproducible compared to that choice.
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin), m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nChains) 
      models[[model]] = m
    }
    save(models, modelnames, file=filename)
  }
}