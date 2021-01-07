data_prevalence_exploration <- function(Y) {
  print("Prevalence")
  print(colMeans(Y>0))
  
  hist(colMeans(Y>0), main="prevalence (proportion of cells occupied)", breaks=20)

  hist(as.matrix(log(Y[Y>0])),main="log abundance conditional on presence")
}

define_hurdle_model <- function(S, X, Y, XFormula, ModelDir, spatial=FALSE, 
                                predictor_selection=FALSE) {
  # predictor_selection: FALSE, spike_slab_jointly, spike_slab_separately, pca, rrr
  
  ofname = "unfitted_models"
  
  # Spatial model
  # Default params
  studyDesign = NULL
  ranLevels = NULL
  if (spatial == TRUE) {
    # Get cell IDs
    S$unique_cell = row.names(S)
    # Create study design DF
    studyDesign = data.frame(cell = as.factor(S$unique_cell))
    # Create spatial info DF
    xy = data.frame(x=S$Longitude,y=S$Latitude)
    rownames(xy) = studyDesign$cell
    print("Setting a spatial random level")
    # If we would define an unstructured random effect, we would use the units = ... argument.
    # Unstructured: with no reference to space or covariates.
    # As we wish to define a spatial random effect, we use instead the sData argument.
    rL.cell = HmscRandomLevel(sData = xy)
    ranLevels = {list("cell" = rL.cell)}
    ofname = paste0(ofname, "_spatial")
  }
  
  # Predictor selection
  # Default params
  XRRRData = NULL
  XRRRFormula = ~. - 1
  ncRRR = 2
  XSelect = NULL
  if (predictor_selection == FALSE) {}
  else if (startsWith(predictor_selection, "spike_slab")) {
    XSelect = list()
    qq = 0.1 # prior probability for a covariate to be included
    for (k in 2:ncol(X)) {
      covGroup = k
      if (endsWith(predictor_selection, "separately")) {
        spGroup = 1:ncol(Y)
      } else if (endsWith(predictor_selection, "jointly")) {
        spGroup = rep(1, ncol(Y))
      }
      q = rep(qq, max(spGroup))
      XSelect[[k-1]] = list(covGroup = covGroup, spGroup = spGroup, q = q)
    }
    if (endsWith(predictor_selection, "separately")) {
      ofname = paste0(ofname, "_spikeSlabSeparately")
    } else if (endsWith(predictor_selection, "jointly")) {
      ofname = paste0(ofname, "_spikeSlabJointly")
    }
  } else if (predictor_selection == "pca") {
    pc = princomp(X)
    X = data.frame(pc$scores[ ,1])
    XFormula = ~.
    ofname = paste0(ofname, "_pca")
  } else if (predictor_selection == "rrr") {
    XFormula = ~1
    XRRRData = X
    XRRRFormula = ~.-1
    ncRRR=1
    ofname = paste0(ofname, "_rrr")
  }

  # Presence Absence data
  # Even if the model definition HMSC accepts logical values (TRUE/FALSE) as the response matrix Y for presence-absence data,
  # some of the post-processing functions assume that Y has numerical 0/1 values. So it is safer to 
  # convert Y to numerical either as Y = 1*Y
  Ypa = 1*(Y>0)
  # Presence-absence model
  m1 = Hmsc(Y=Ypa, XData = X,  XFormula = XFormula,
            distr="probit",
            studyDesign=studyDesign,
            ranLevels=ranLevels,
            XSelect = XSelect,
            XRRRData = XRRRData, XRRRFormula = XRRRFormula, ncRRR = ncRRR)
  print(m1)
  print(head(m1$X))
  
  # Log-normal abundance data conditional on presence
  Yabu = Y
  Yabu[Yabu==0] = NA
  Yabu=log(Yabu)
  # Log-normal model for abundance conditional on presence
  m2 = Hmsc(Y=Yabu, YScale = TRUE,
            XData = X,  XFormula = XFormula,
            distr="normal",
            studyDesign=studyDesign,
            ranLevels=ranLevels,
            XSelect = XSelect,
            XRRRData = XRRRData, XRRRFormula = XRRRFormula, ncRRR = ncRRR)
  print(m2)
  print(head(m2$X))
  
  # Hurdle model
  models = list(m1,m2)
  modelnames = c("presence_absence","abundance_COP")
  print(paste("Constructed models:", modelnames))
  
  # Save unfitted models
  ofname = file.path(ModelDir, ofname)
  print(paste("Save Hurdle model in:", ofname))
  save(models, modelnames, file = ofname)
}