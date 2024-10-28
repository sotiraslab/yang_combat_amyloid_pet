# Simulation experiments

## Simulating from LMER

See `lmer_simulation.R` for implementation

## Full simulation

To simulate clinical trial data, we will take the following approach:

1. **Fit ComBat model on healthy control data**
  - we fit ComBat to learn both the covariate coefficients and the batch-specific parameters
  - covariates to include: age, sex, APOE
    - here, we include age as a covariate, even though the treatment effect introduced will be a longitudinal rate of change. This is because the rate of change effect will be related to *time from treatment*, not age
  - TODO: fit on amyloid-positive data
2. **Generate subject data**
  - using the learned (EB adjusted) batch parameters of each tracer and the covariate parameters, we generate new data
  - for the treatment group, we further introduce a treatment effect that is dependent on *time from treatment* (or time from baseline)
3. **Use ComBat to remove batch effects**
  - learn ComBat model on synthesized data, apply on synthesized data
  - use CV scheme
4. **Test for group differences**
  - test for significant longitudinal parameter

### Notes on `neuroharmonize` package

1. Standardization step
 - `B_hat`: contains OLS parameters for each *covariate* and *batch*
   - shape: (n_covar + n_batch, n_feature)
   - $\hat{\beta}$ in the ComBat formula
 - `grand_mean`: weighted mean of OLS estimated batch shift across all batches
   - weighted by number of samples in each batch
   - $\hat{\alpha}$, or the mean across batches, in the ComBat formula
   - shape: (n_feature,)
 - `var_pooled`: pooled variance of data as compared to the grand mean
   - in other words, this is the mean squared error of the residuals of the OLS model
   - $\hat{\sigma}$ in the ComBat formula
   - used to scale the data in the standardization step (see section 3.1 of Johnson *et al.*)
   - shape: (n_feature, 1)
2. Empirical Bayes estimation and adjustment
 - functions `aprior` and `bprior` are used to estimate the prior distribution of batch scale effects, which is modeled by an inverse gamma distribution
   - method of moments is used to estimate the mean $\lambda$ (`aprior`) and variance $\theta$ (`bprior`) of this distribution with a closed form expression
 - `gamma_star`: EB adjusted batch shift parameter
 - `delta_star`: EB adjusted batch scale parameter
 - `gamma_hat`: mean standardized measurement for each batch
   - the empirical batch shift parameter
   - this is computed by taking the average within-batch standardized measurement $Z_{ijg}$
     - in the code, they compute this in a linear algebra way by solving a linear system: $Z = B \hat{gamma}$
     - $Z$ are the standardized measurements, and $B$ is the batch design matrix
     - use the pseudoinverse to obtain LS solution
     - essentially the same as computing the within-batch mean, but does it in one linear algebra step
 - `delta_hat`: variance of standardized measurement for each batch
   - the empirical batch scale parameter
 - `gamma_bar`: mean across all empirical gammas
   - set as the mean of the normal distribution of the prior on $\gamma$
 - `t2`: variance across all empirical gammas
   - $\bar{\tau}^2$ in the ComBat formula
   - set as the variance of the normal distribution of the prior on $\gamma$
 - `a_prior`: mean of inverse gamma distribution of the prior on $\delta$
   - $\bar{\lambda}$
 - `b_prior`: variance of inverse gamma distribution of the prior on $\delta$
   - $\bar{\theta}$ in the ComBat formula

### Experiment 1: generate longitudinal data, test for difference in rate-of-change using LRT

- simulate data using `n_subj` and `p_av45` parameters with `full_simulation.py`
- learn ComBat model, apply to simulated data
  - test out both ComBat with no covariates and ComBat with covariates
  - when using covariates, include treatment group as a covariate
- fit LMER to data in `full_simulation.R`
- test for significant group differences in rate-of-change using LRT
- repeat 1000 times

### Experiment 2: generate 2-timepoint data, test for difference in annual rate-of-change using t-test

- simulate data using `n_subj` and `p_av45` parameters with `full_simulation.py`
- learn ComBat model, apply to simulated data
  - test out both ComBat with no covariates and ComBat with covariates
  - when using covariates, include treatment group as a covariate
- compute annual rate-of-change between baseline and follow-up ($\frac{SUVR_2 - SUVR_1}{\Delta t}$)
- test for significant group differences in annual rate-of-change using t-test
