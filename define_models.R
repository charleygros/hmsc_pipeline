data_prevalence_exploration <- function(Y) {
  print("Prevalence")
  print(colMeans(Y>0))
  
  hist(colMeans(Y>0), main="prevalence (proportion of cells occupied)", breaks=20)

  hist(as.matrix(log(Y[Y>0])),main="log abundance conditional on presence")
}

define_hurdle_model <- function(S, X, Y, XFormula, ModelDir) {
  # Get cell IDs
  S$unique_cell = row.names(S)
  # Create study design DF
  studyDesign = data.frame(cell = as.factor(S$unique_cell))
  # Create spatial info DF
  xy = data.frame(x=S$Longitude,y=S$Latitude)
  rownames(xy) = studyDesign$cell
  print("Setting a spatial random level")
  rL.cell = HmscRandomLevel(sData = xy)
  
  # Presence Absence data
  Ypa = 1*(Y>0)
  # Presence-absence model
  m1 = Hmsc(Y=Ypa, XData = X,  XFormula = XFormula,
            distr="probit",
            studyDesign=studyDesign,
            ranLevels={list("cell" = rL.cell)})
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
            ranLevels={list("cell" = rL.cell)})
  print(head(m2$X))
  
  # Hurdle model
  models = list(m1,m2)
  modelnames = c("presence_absence","abundance_COP")
  print(paste("Constructed models:", modelnames))
  
  # Save unfitted models
  ofname = file.path(ModelDir, "unfitted_models")
  print(paste("Save Hurdle model in:", ofname))
  save(models, modelnames, file = ofname)
}