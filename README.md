# hmsc_pipeline
R utility functions to define, fit and evaluate HMSC models. These functions build on top of [HMSC framework](https://github.com/hmsc-r/HMSC) to provide a full Species Distribution Modelling pipeline.

## Available features

WORK IN PROGRESS

1. Dataset sanity checks --> `read_data.R`
2. Create HMSC model instance. Different model are available, including spatial models or Hurdle modelling approach. Possibility to perform predictor selection (`spike_slab_jointly`, `spike_slab_separately`, `pca`, `rrr`). --> `define_model.R`
3. Fit the model by sampling the posterior with block-conditional Gibbs MCMC sampler. --> `fit_model.R`
4. Evaluage MCMC convergence. By looking at the effective size of the posterior sample: beta-parameters (species niches) and V-parameters (variation in species niches). By plotting the Gelman diagnostics, i.e. the Potential scale reduction factors. --> `evaluate_convergence.R`
5. Compute model fit, both explanatory and predictive powers. R2, TjurR2, RMSE, AUC are computed. A 5-fold cross validation is used for the predictive power. --> `compute_model_fit.R`
6. Create a grid of environmental predictors by extracting data from a larger raster --> `create_grid.R`

## References
See [HMSC Github documentation section](https://github.com/hmsc-r/HMSC#documentation).

## Acknowledgements
- N. Hill for `hmsc_CVpred_parallel` and `env_niches` functions.
- HMSC team for the fantastic 2020 HMSC online workshop.
