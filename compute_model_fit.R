hmsc_CVpred_parallel <- function (hM, partition = NULL,  start = 1, thin = 1, 
                                  Yc = NULL, mcmcStep = 1, expected = TRUE, 
                                  initPar = NULL, nParallel = 1, 
                                  nChains = length(hM$postList), 
                                  updater = list(), verbose = hM$verbose, 
                                  alignPost = TRUE) 
{
  ###
  # Cross validation code to allow more efficient parallel processing at level of nfold.
  # Modified Hmsc package computePredictedValues function. Arguments are the same.
  ###
  
  train = (partition != k)
  val = (partition == k)
  dfPi = as.data.frame(matrix(NA, sum(train), hM$nr), 
                       stringsAsFactors = TRUE)
  colnames(dfPi) = hM$rLNames
  for (r in seq_len(hM$nr)) {
    dfPi[, r] = factor(hM$dfPi[train, r])
  }
  switch(class(hM$X)[1L], matrix = {
    XTrain = hM$X[train, , drop = FALSE]
    XVal = hM$X[val, , drop = FALSE]
  }, list = {
    XTrain = lapply(hM$X, function(a) a[train, , 
                                        drop = FALSE])
    XVal = lapply(hM$X, function(a) a[val, , drop = FALSE])
  })

  XRRRTrain = NULL
  XRRRVal = NULL

  hM1 = Hmsc(Y = hM$Y[train, , drop = FALSE], X = XTrain, 
             XRRR = XRRRTrain, ncRRR = hM$ncRRR, XSelect = hM$XSelect, 
             distr = hM$distr, studyDesign = dfPi, Tr = hM$Tr, 
             C = hM$C, ranLevels = hM$rL)
  setPriors(hM1, V0 = hM$V0, f0 = hM$f0, mGamma = hM$mGamma, 
            UGamma = hM$UGamma, aSigma = hM$aSigma, bSigma = hM$bSigma, 
            nu = hM$nu, a1 = hM$a1, b1 = hM$b1, a2 = hM$a2, 
            b2 = hM$b2, rhopw = hM$rhowp)
  hM1$YScalePar = hM$YScalePar
  hM1$YScaled = (hM1$Y - matrix(hM1$YScalePar[1, ], 
                                hM1$ny, hM1$ns, byrow = TRUE))/matrix(hM1$YScalePar[2, 
                                ], hM1$ny, hM1$ns, byrow = TRUE)
  hM1$XInterceptInd = hM$XInterceptInd
  hM1$XScalePar = hM$XScalePar
  switch(class(hM$X)[1L], matrix = {
    hM1$XScaled = (hM1$X - matrix(hM1$XScalePar[1, 
    ], hM1$ny, hM1$ncNRRR, byrow = TRUE))/matrix(hM1$XScalePar[2, 
    ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
  }, list = {
    hM1$XScaled = list()
    for (zz in seq_len(length(hM1$X))) {
      hM1$XScaled[[zz]] = (hM1$X[[zz]] - matrix(hM1$XScalePar[1, 
      ], hM1$ny, hM1$ncNRRR, byrow = TRUE))/matrix(hM1$XScalePar[2, 
      ], hM1$ny, hM1$ncNRRR, byrow = TRUE)
    }
  })
  if (hM1$ncRRR > 0) {
    hM1$XRRRScalePar = hM$XRRRScalePar
    hM1$XRRRScaled = (hM1$XRRR - matrix(hM1$XRRRScalePar[1, 
    ], hM1$ny, hM1$ncORRR, byrow = TRUE))/matrix(hM1$XRRRScalePar[2, 
    ], hM1$ny, hM1$ncORRR, byrow = TRUE)
  }
  hM1$TrInterceptInd = hM$TrInterceptInd
  hM1$TrScalePar = hM$TrScalePar
  hM1$TrScaled = (hM1$Tr - matrix(hM1$TrScalePar[1, 
  ], hM1$ns, hM1$nt, byrow = TRUE))/matrix(hM1$TrScalePar[2, 
  ], hM1$ns, hM1$nt, byrow = TRUE)
  hM1 = sampleMcmc(hM1, samples = hM$samples, thin = hM$thin, 
                   transient = hM$transient, adaptNf = hM$adaptNf, 
                   initPar = initPar, nChains = nChains, nParallel = nParallel, 
                   updater = updater, verbose = verbose, alignPost = alignPost)
  postList = poolMcmcChains(hM1$postList, start = start)
  dfPi = as.data.frame(matrix(NA, sum(val), hM$nr), 
                       stringsAsFactors = TRUE)
  colnames(dfPi) = hM$rLNames
  for (r in seq_len(hM$nr)) {
    dfPi[, r] = factor(hM$dfPi[val, r])
  }
  
  pred1 = predict(hM1, post = postList, X = XVal, 
                  XRRR = XRRRVal, studyDesign = dfPi, Yc = Yc[val, 
                                                              , drop = FALSE], mcmcStep = mcmcStep, expected = expected)
  pred1Array = abind(pred1, along = 3)
  
  return(pred1Array)
}


compute_model_fit <- function(modelDir, fname_model) {
  filename_out = paste("MF_", basename(fname_model), sep = "")
  filename_out = file.path(modelDir, filename_out)
  if (file.exists(filename_out)) {
    print(paste("File already exists -- Skipping:", filename_out))
    
  } else {
    load(fname_model)
    nm = length(models)

    MF = list()
    MFCV = list()
    WAIC = list()
    
    for(i in 1:nm){
      model_name_cur = modelnames[[i]]
      print(paste0("model = ",as.character(model_name_cur)))
      m = models[[i]]
      
      # Explanatory power
      preds = computePredictedValues(m)
      MF[[model_name_cur]] = evaluateModelFit(hM=m, predY=preds)
      
      # Predictive power
      partition = createPartition(m, nfolds = 5)
      preds = computePredictedValues(m,partition=partition)
      # TODO: CORRECT k variable w Nicole
      #preds = hmsc_CVpred_parallel(m, partition=partition)
      MFCV[[model_name_cur]] = evaluateModelFit(hM=m, predY=preds)
      WAIC[[model_name_cur]] = computeWAIC(m)
    }
    save(MF, MFCV, WAIC, modelnames, file = filename_out)
  }
}