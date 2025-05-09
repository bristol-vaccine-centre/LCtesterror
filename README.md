
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LCtesterror

<!-- badges: start -->
<!-- badges: end -->

The goal of LCtesterror is to estimate true prevalence from test data
accounting for unknown test error.

## Installation

You can install the development version of LCtesterror from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("bristol-vaccine-centre/LCtesterror")
```

<!-- Need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. 
If adding plots need to commit and push resulting figure files, so they display on GitHub and CRAN.
-->

## Example

Example:

``` r
library(LCtesterror)

 # Example (no test dependence, delay, or time included)

 # Simulate data
 sim_data <- sim.test.data(sim_size = 5000)
 
 # Simulated input parameters
 sim_data$test_parameters
 # Simulated true prevalence = 21.4%
 
 # Simulated test data
 head(sim_data$test_results)

# Run LC model using simulated test results to correct prevalence for test error
 fit <- run.LC.model(sim_data$test_results, num_tests = 4,
                       data_ID = "sim", model_name = "basic_sim",
                      dependency_groups = NULL, covariates = NULL,
                      iter=2000, chains=4, warmup=1000)

 # View model outputs
 
 # Estimated prevalence
 fit$prev_mean
 # The model correctly infers the simulated true prevalence
 
 # Estimated test sensitivity and specificity 
 fit$median_sens
 fit$median_spec
 
 # Stan model traceplots
 fit$traceplots
 
```

    #> Using LCtesterror version: 0.0.0.9000
