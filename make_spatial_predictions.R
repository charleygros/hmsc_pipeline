make_spatial_predictions <- function(S, X, model_path, fname_out) {
  #' Make predictions from a spatial grid.
  #'
  #' This function load a model and make some predictions given some predictors
  #' and spatial coordinates. The predictions are then saved under fname_out.
  #' 
  #' @param S Study design dataframe containing "Longitude" and "Latitude"
  #' @param X Predictors dataframe
  #' @param model_path Path towards the fitted models
  #' @param fname_out Predictions of each model will be saved under this 
  #' filename
  #' 
  #' @note TODO: Adapt to non spatial model
  
  # Construct the object xy.grid that have the coordinates
  # Get cell IDs
  grid$unique_cell = row.names(grid)
  # Create study design DF
  studyDesign = data.frame(cell = as.factor(grid$unique_cell))
  # Create spatial info DF
  xy.grid = data.frame(x=grid$Longitude, y=grid$Latitude)
  rownames(xy.grid) = studyDesign$cell
  
  # Construct the object XData.grid that have the environmental predictors
  XData.grid = X
  
  # Load model
  load(model_path)
  
  # Init results container
  results <- list()
  
  # Loop across models
  nm = length(models)
  for (model_idx in 1:nm) {
    print(paste0("Predictions of ", modelnames[model_idx]))
    m = models[[model_idx]]
    
    # Use the prepareGradient function to convert the environmental and spatial predictors into a format that can be used as input for the predict function.
    gradient = prepareGradient(m,
                               XDataNew = XData.grid,
                               sDataNew = list(cell=xy.grid))
    
    # Compute the posterior predictive distribution.
    predY = predict(m, Gradient=gradient, expected = TRUE, nParallel = nChains)
    results[[modelnames[model_idx]]] <- predY
  }
  
  # Save results
  print(paste0("Saving results in: ", fname_out
               ))
  save(results, file=fname_out)
}