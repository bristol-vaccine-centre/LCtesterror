
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
Knits using installed version of LCtesterror.
Example eval = FALSE because takes a while to run. Changed to TRUE for manual knit.
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
#> # A tibble: 4 × 12
#>   test_id  sens  spec p_performed sens_sim disease_prev_sim test_positivity_sim
#>   <chr>   <dbl> <dbl>       <dbl>    <dbl>            <dbl>               <dbl>
#> 1 test1    0.99  0.99         1       0.99            0.198               0.208
#> 2 test2    0.99  0.99         1       0.99            0.198               0.201
#> 3 test3    0.98  0.98         0.8     0.98            0.198               0.203
#> 4 test4    0.98  0.98         0.8     0.98            0.198               0.207
#> # ℹ 5 more variables: test_coverage_sim <dbl>, disease_prev_RG_est <dbl>,
#> #   overall_test_positivity_sim <dbl>, CI_lower <dbl>, CI_upper <dbl>
 # Simulated true prevalence = 20%
 
 # Simulated test data
 head(sim_data$test_results)
#> # A tibble: 6 × 4
#>   test1 test2 test3 test4
#>   <dbl> <dbl> <dbl> <dbl>
#> 1     1     1     1     1
#> 2     0     0     0     0
#> 3     0     0     0     0
#> 4     0     0     0     1
#> 5     0     0     0     0
#> 6     0     0     0     0

# Run LC model using simulated test results to correct prevalence for test error
 fit <- run.LC.model(sim_data$test_results, num_tests = 4,
                       data_ID = "sim", model_name = "basic_sim",
                      dependency_groups = NULL, covariates = NULL,
                      iter=2000, chains=4, warmup=1000, plot_chains = FALSE)

 
 # View model outputs
 # Estimated prevalence
 fit$prev_mean
#> [1] "Mean prevalence: 19.8% (18.7-20.9%)"
 # The model correctly infers the simulated true prevalence
 
 # Estimated test sensitivity and specificity 
 fit$median_sens
#>                   mean      2.5%     97.5%
#> Se_median[1] 0.9941189 0.9894686 0.9975216
#> Se_median[2] 0.9916728 0.9861221 0.9958106
#> Se_median[3] 0.9847242 0.9766545 0.9912738
#> Se_median[4] 0.9890827 0.9824261 0.9943707
 fit$median_spec
#>                   mean      2.5%     97.5%
#> Sp_median[1] 0.9909279 0.9883872 0.9931538
#> Sp_median[2] 0.9955418 0.9937000 0.9970597
#> Sp_median[3] 0.9880344 0.9846469 0.9909931
#> Sp_median[4] 0.9889894 0.9857982 0.9918422
 
 # Stan model traceplots
 fit$traceplots
#> $traceplot_lp_warmup
```

<img src="man/figures/README-example-1.png" width="100%" />

    #> 
    #> $traceplot_lp

<img src="man/figures/README-example-2.png" width="100%" />

    #> 
    #> $stan_dens_prev

<img src="man/figures/README-example-3.png" width="100%" />

    #> 
    #> $stan_dens_Se_mean

<img src="man/figures/README-example-4.png" width="100%" />

    #> 
    #> $stan_dens_Sp_mean

<img src="man/figures/README-example-5.png" width="100%" />

    #> 
    #> $stan_dens_Se_median

<img src="man/figures/README-example-6.png" width="100%" />

    #> 
    #> $stan_dens_Sp_median

<img src="man/figures/README-example-7.png" width="100%" />

    #> Using LCtesterror version: 0.0.0.9000
