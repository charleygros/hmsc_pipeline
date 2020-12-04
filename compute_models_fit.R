compute_models_fit <- function(modelDir, samples_list, thin_list, nChains) {
  nst = length(thin_list)

  for (Lst in 1:nst) {
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    
    filename_out = paste("MF_thin_", as.character(thin),
                         "_samples_", as.character(samples),
                         "_chains_",as.character(nChains),
                         ".Rdata",sep = "")
    filename_out = file.path(modelDir, filename_out)
    if (file.exists(filename_out)) {
      print(paste("File already exists -- Skipping:", filename_out))
      
    } else {
      filename = paste("models_thin_", as.character(thin),
                       "_samples_", as.character(samples),
                       "_chains_",as.character(nChains),".Rdata",sep = "")
      filename = file.path(modelDir, filename)
      load(filename)
      nm = length(models)
      print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
      
      MF = list()
      MFCV = list()
      WAIC = list()
      
      for(model in 1:nm){
        print(paste0("model = ",as.character(model)))
        m = models[[model]]
        preds = computePredictedValues(m)
        MF[[model]] = evaluateModelFit(hM=m, predY=preds)
        partition = createPartition(m, nfolds = 2)
        preds = computePredictedValues(m,partition=partition) #nParallel = nChains
        MFCV[[model]] = evaluateModelFit(hM=m, predY=preds)
        WAIC[[model]] = computeWAIC(m)
      }
      
      save(MF,MFCV,WAIC,modelnames,file = filename_out)
    }
  }
}