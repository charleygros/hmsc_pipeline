# hmsc_pipeline
R utility functions to define, fit and evaluate HMSC models. These functions build on top of [HMSC framework](https://github.com/hmsc-r/HMSC) to provide a full Species Distribution Modelling pipeline. See `demo.pdf` for an illustration of these tools.

## Requirements

- `Hmsc`
- `ggplot2`
- `raster`
- `vioplot`
- `docstring`

## Available features

WORK IN PROGRESS

1. Dataset sanity checks --> `read_data.R`
2. Create HMSC model instance. Different model are available, including spatial models or Hurdle modelling approach. Possibility to perform predictor selection (`spike_slab_jointly`, `spike_slab_separately`, `pca`, `rrr`). --> `define_model.R`
3. Fit the model by sampling the posterior with block-conditional Gibbs MCMC sampler. --> `fit_model.R`
4. Evaluage MCMC convergence. By looking at the effective size of the posterior sample: beta-parameters (species niches) and V-parameters (variation in species niches). By plotting the Gelman diagnostics, i.e. the Potential scale reduction factors. --> `evaluate_convergence.R`
5. Compute model fit, both explanatory and predictive powers. R2, TjurR2, RMSE, AUC are computed. A 5-fold cross validation is used for the predictive power. --> `compute_model_fit.R`
6. XX `make_predictions.R`
7. Create a grid of environmental predictors by extracting data from a larger raster --> `create_grid.R`
8. XX `make_spatial_predictions.R`
9. XX `show_env_niche.R`
10. XX `show_parameter_estimates.R`
11. XX `show_models_fit.R`

## Help

To get functions' description, you can use the following:

```r
?make_spatial_predictions
```

which provides the function description:
```r
Description
This function load a model and make some predictions given some predictors and spatial coordinates. The predictions are then saved under fname_out.

Usage
make_spatial_predictions(S, X, model_path, fname_out)

Arguments
S	: Study design dataframe containing "Longitude" and "Latitude"

X: Predictors dataframe

model_path: Path towards the fitted models

fname_out: Predictions of each model will be saved under this filename

Note:
TODO: Adapt to non spatial model
```

## References
See [HMSC Github documentation section](https://github.com/hmsc-r/HMSC#documentation).

## Acknowledgements
- N. Hill for `hmsc_CVpred_parallel` and `env_niches` functions.
- HMSC team for the fantastic 2020 HMSC online workshop.
