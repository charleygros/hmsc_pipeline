# hmsc_pipeline
R functions to define, fit and evaluate HMSC models

## Available features

- Dataset sanity checks --> `read_data.R`
- Compute model fit, both explanatory and predictive powers. R2, TjurR2, RMSE, AUC are computed. A 5-fold cross validation is used for the predictive power. --> `compute_model_fit.R`
- Create a grid of environmental predictors by extracting data from a larger raster --> `create_grid.R`
- XX
