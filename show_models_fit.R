show_models_fit <- function(modelDir, samples_list, thin_list, nChains, panelsDir, poissonModel = FALSE) {
  # TODO: Implement pseudoR2 for poissonModel = TRUE
  
  nst = length(thin_list)
  
  for (Lst in 1:nst) {
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    
    filenameout = paste("modelFit_thin_", as.character(thin),
                        "_samples_", as.character(samples),
                        "_chains_",as.character(nChains),
                        ".pdf",sep = "")
    filenameout = file.path(panelsDir, filenameout)
    if (file.exists(filenameout)) {
      print(paste("File already exists -- Skipping:", filenameout))
      
    } else {
      filename = paste("MF_thin_", as.character(thin),
                       "_samples_", as.character(samples),
                       "_chains_",as.character(nChains),
                       ".Rdata",sep = "")
      filename = file.path(modelDir, filename)
      load(filename)
    
      nm = length(MF)
      
      pdf(file = filenameout)
      for(j in 1:nm){
        cMF = MF[[j]]
        cMFCV = MFCV[[j]]
        if(!is.null(cMF$TjurR2)){
          plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
               xlab = "explanatory power",
               ylab = "predictive power",
               main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": Tjur R2"))
          abline(0,1)
          abline(v=0)
          abline(h=0)
        }
        if(!is.null(cMF$R2)){
          plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
               xlab = "explanatory power",
               ylab = "predictive power",
               main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": R2"))
          abline(0,1)
          abline(v=0)
          abline(h=0)
        }
        if(!is.null(cMF$AUC)){
          plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
               xlab = "explanatory power",
               ylab = "predictive power",
               main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": AUC"))
          abline(0,1)
          abline(v=0.5)
          abline(h=0.5)
        }
      }
      dev.off()
    }
  }
}

