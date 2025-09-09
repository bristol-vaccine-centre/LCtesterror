
#' @title Runs multiple latent class model simulations over time with different parameter values.
#' @description Runs multiple simulations of a bayesian LC stan model to infer true disease prevalence over time for a seasonal disease, test sensitivity, and specificity from
#' simulated test data with different parameters values (sensitivity, specificity, and probability of the test occurring).
#' Can simulate multiple test results but all with identical parameters.
#' Can run different parameter combinations in parallel using furrr::future_map. Need to first define future and number of cores to parallise (workers) using future::plan(future::multisession, workers= ).
#' This is separate to rstan parallel processing of chains which uses options(mc.cores = num_chains)
#' Parallel processing i.e workers >1 will not work if using devtools::load_all()
#'
#' @param num_tests A numeric value for the number of tests to simulate and model.
#' @param days Number of days to simulate. Default = 365.
#' @param spec_vec Values of specificity for all tests as a numeric vector for each simulation. Default = `c(1)`.
#' @param sens_vec Values of sensitivity for all tests as a numeric vector for each simulation. Default = `c(1)`.
#' @param p_performed_vec Probability that each test is performed as a numeric vector for each simulation. Default = `c(1)`.
#' @param sim_size Number of individuals to simulate in each dataset. Default = 1000.
#' @param iter The number of iterations for the stan model. Default = 1000.
#' @param chains The number of chains for the stan model. Default = 4.
#' @param warmup The number of warmup iterations for the stan model. Default = 500.
#' @param stan_arg Optional extra arguments to pass to the rstan::sampling function. Default = NULL.
#' @param data_ID Optional character identifier for labeling data outputs. Default = `"sims"`.
#' @param prior_spec Specification of specificity prior. A list of length equal to the value of num_tests with each element containing a vector of length two specifying the alpha and beta parameters for the Beta prior. Default = c(10, 1) for each test.
#' @param prior_sens Specification of sensitivity prior. A list of length equal to the value of num_tests with each element containing a vector of length two specifying the alpha and beta parameters for the Beta prior. Default = c(1, 1) for each test.
#' @param other_priors List of any other priors to be specified differently to the defaults. Given as character strings as written for stan. Defaults =
#' list(RE_prior= "normal(0,1), "bpos_prior= "gamma(1,1)", bneg_prior= "gamma(1,1)",
#' gaussian_prev_amplitude_prior= "beta(1,1)", gaussian_prev_baseline_prior= "beta(1,1)", mean_gaussian_prior= "uniform(0,52)", sigma_gaussian_prior = "normal(0,10)",
#' exp_prev_baseline_prior= "beta(1,1)", exp_growth_rate_prior= "gamma(1,5)")
#' @param seed Random seed for set.seed(). Default = 953.
#' @param Est_R_window Size, in number of days, of sliding window for custom R(t) estimation. Default = 14.
#' @param Est_R_n_samples Number of samples for uncertainty estimation in custom R(t). Default = 1000.
#' @param mean_gi Mean generation interval used for estimating R(t) from EpiEstim. Default = 1/gamma.
#' @param max_t Maximum time (days) used for calculating the gi_distribution for estimating R(t) from EpiEstim. Default = 5 * mean_gi
#' @param years Number of years to run the SIR model. Default = 50.
#' @param N Population size used in the SIR model. Default = 1.
#' @param init Initial state of the SIR model as a list. Default = init = list(init_S = 0.99, init_I = 0.01, init_R = 0).
#' @param SIR_params list of SIR parameters to simulate prevalence over time for a seasonal disease. Includes:
#' \describe{
#'   \item{beta_0}{Baseline transmission rate}
#'   \item{desired_R0}{Desired basic reproduction number (R0) - used to calculate beta_0 if beta_0 is null.}
#'   \item{beta_1}{Amplitude of seasonal forcing}
#'   \item{phi}{Phase shift of seasonal forcing}
#'   \item{gamma}{Recovery rate}
#'   \item{omega}{Waning immunity rate}
#' Default = beta0 = NULL, desired_R0 = 2.5, beta1 = 0.07, phi = 1.5, gamma = 0.03, omega = 0.001.
#' If beta0 is null it is calculated as desired_R0 * gamma.
#' }
#' @param exp_params list of exponential model parameters to simulate prevalence over time at the start of an epidemic. Includes:
#' \describe{
#'   \item{I0}{Initial disease prevalence}
#'   \item{beta_0}{Baseline transmission rate}
#'   \item{gamma}{Recovery rate}
#' Default = I0 = 0.0001, beta0 = 0.06, gamma = 0.04.
#' }
#' @param time_model If covariates includes "Time", which model should be used to infer changing prevalence over time - "gaussian" (for a seasonal peak) or "exponential". Default = "gaussian".
#' @return Stan model fit and various summary outputs.
#' A list containing:
#' \describe{
#'   \item{sim_inputs}{A data frame of all simulation input parameters per model run.
#'   Output from sim.test.data.time()$test_parameters}
#'    \item{sim_data}{A data frame of all simulated data per model run.
#'    Output from sim.test.data.time()$sim_data}
#'   \item{stan_results_df}{A data frame of posterior estimates and 95% credible intervals
#'     for sensitivity, specificity, and prevalence for each simulation for each week simulated.
#'     Sensitivity and specificity values are repeated for each week but are not time-varying parameters in this model.}
#'   \item{divergence_summary}{A data frame summarising any divergent transitions or low ESS in the Stan models.
#'   Output from check.divergent.transitions().}
#'   \item{R_estimates}{A data frame summarising of different R estimates and growth rate estimates combined from sim.test.data.time() function outputs}
#' }
#'
#' @details
#' For each combination of sensitivity, specificity, prevalence, and test performance probability,
#' the function:
#' \enumerate{
#'   \item Simulates test results for a population of individuals using the `sim.test.data()` function.
#'   \item Fits a Bayesian latent class model using `run.LC.model()`.
#'   \item Extracts posterior summaries for each test's sensitivity and specificity, as well as disease prevalence.
#'   \item Records model diagnostics, including checking for divergent transitions.
#' }
#'
#' Simulated test parameters and model estimates are stored in separate data frames, and intermediate
#' variables are cleared after each run to conserve memory.
#'
#' @importFrom rstan summary sampling rstan_options traceplot stan_dens
#' @importFrom dplyr filter group_by mutate rename row_number all_of starts_with n_distinct relocate summarise n select everything
#' @importFrom tidyr pivot_longer pivot_wider matches unite
#' @importFrom gt gt tab_header
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_density geom_point geom_ribbon labs theme_minimal theme_bw
#' @importFrom stats median na.omit setNames
#' @importFrom graphics pairs
#' @importFrom utils globalVariables
#' @importFrom purrr imap reduce
#' @importFrom future plan
#' @importFrom furrr future_map furrr_options
#' @name run.sims.LC.time
#' @examples
#' if (interactive()) {
#'  results <- run.sims.LC.time(num_tests = 4,
#'                        spec_vec = c(1),
#'                        sens_vec = c(1),
#'                        p_performed_vec= c(1),
#'                        sim_size = 100,
#'                        chains = 2,
#'                        data_ID = "sims_time")
#' }
#'
#'
#Function for running lots of LC runs using sim data with different parameters

utils::globalVariables(c("ess_threshold", "matches", "unite", "metric", "stat", "metric_stat",
                         "sim_sens", "sim_spec", "sim_prob", "everything"))

#' @export run.sims.LC.time
run.sims.LC.time <- function(num_tests, days = 365, spec_vec= c(1), sens_vec= c(1), p_performed_vec= c(1),
                                      sim_size=1000, iter=1000, chains=4, warmup=500, stan_arg=list(), data_ID = "sims",
                                      prior_spec = NULL, prior_sens = NULL,
                                      other_priors  = list(),
                                      seed = 953,
                                      Est_R_window = 14, Est_R_n_samples = 1000, #for my R method
                                      mean_gi = NULL, max_t = NULL, #For EpiEstim - Generation interval (mean and SD)
                                      years = 50, N = 1, init = list(init_S = 0.99,init_I = 0.01,init_R = 0),
                                      SIR_params =list(beta0 = NULL, desired_R0 = 2.5, beta1 = 0.07, phi = 1.5, gamma = 0.03, omega = 0.001), #for SIR
                                      exp_params = list(I0 = 0.0001, beta0 = 0.06, gamma = 0.04),
                                      time_model = "gaussian"
) {

  set.seed(seed)

  sim_stan_result_df_ALL <- data.frame()
  sim_inputs_df_ALL <- data.frame()
  sim_data_df_ALL <- data.frame()
  divergence_summary_df <- data.frame()
  R_est_df_ALL <- data.frame()

  param_grid <- expand.grid(
    sens = sens_vec,
    spec = spec_vec,
    p_performed = p_performed_vec,
    stringsAsFactors = FALSE
  )

  message(paste("Time_model", paste(time_model, collapse=", ")))

  # Wrapper function to run simulation/model for one parameter combo
  run.for.param.combo <- function(param_row) {

    s <- param_row$sens
    c <- param_row$spec
    r <- param_row$p_performed

    message("s = ", s, ", c = ", c, ", r = ", r)
    message(paste("param_row:", paste(param_row, collapse=", ")))

    # Construct test_params list
    test_params <- list()
    for (i in 1:num_tests) {
      testname <- paste0("test", i)
      test_params[[testname]] <- list(sens = s, spec = c, p_performed = r)
    }

    # Run model for each parameter combo:
    # Simulate data
  if (time_model == "exponential") {

    sim_data <- sim.test.data.exp(sim_size=sim_size, days=days, test_params=test_params,
                                   seed=seed, Est_R_window = Est_R_window, Est_R_n_samples = Est_R_n_samples,
                                   mean_gi = mean_gi, max_t = max_t,
                                   params =exp_params)

    result_name <- paste("sens_", s, "_spec_", c, "_p_", r, sep = "")

    stopifnot("day_of_year" %in% colnames(sim_data$test_results))

    #converts time in sim data from days to weeks
    test_results <- sim_data$test_results
    test_results <- test_results %>%
      dplyr::mutate(Time = ((day_of_year - 1) %/% 7) + 1 ) %>%
      dplyr::select(-day_of_year)

    stopifnot("Time" %in% colnames(test_results))

    # Run model
    result <- run.LC.model(data=test_results, num_tests=num_tests,
                           iter=iter, chains=chains, warmup=warmup, stan_arg=stan_arg,
                           data_ID = data_ID, prior_spec = prior_spec, prior_sens = prior_sens,
                           other_priors = other_priors,
                           covariates = c("Time"), time_model = time_model)


  } else {

    sim_data <- sim.test.data.time(sim_size=sim_size, days=days, test_params=test_params,
                                   seed=seed, Est_R_window = Est_R_window, Est_R_n_samples = Est_R_n_samples,
                                   mean_gi = mean_gi, max_t = max_t,
                                   years = years, N = N, init = init,
                                   params =SIR_params)

    result_name <- paste("sens_", s, "_spec_", c, "_p_", r, sep = "")

    stopifnot("day_of_year" %in% colnames(sim_data$test_results))

    #converts time in sim data from days to weeks
    test_results <- sim_data$test_results
    test_results <- test_results %>%
      dplyr::mutate(Time = ((day_of_year - 1) %/% 7) + 1 ) %>%
      dplyr::select(-day_of_year)

    stopifnot("Time" %in% colnames(test_results))

    # Run model
    result <- run.LC.model(data=test_results, num_tests=num_tests,
                           iter=iter, chains=chains, warmup=warmup, stan_arg=stan_arg,
                           data_ID = data_ID, prior_spec = prior_spec, prior_sens = prior_sens,
                           covariates = c("Time"))

  }

    stan_fit <- result$stan_fit
    fit_summary <- rstan::summary(stan_fit)
    stan_fit_summary <- fit_summary$summary
    fit_summary_df <- as.data.frame(stan_fit_summary)


    # Function to check model for divergent transitions
    divergence_check <- check.divergent.transitions(stan_fit, result_name)


    sim_stan_result_df <- data.frame()


    # Extract prevalence with CIs
    #Prev over time
    prev_time_rows <- grepl("^weekly_prevalence", rownames(fit_summary_df))
    stan_fit_summary_prev_time <- subset(fit_summary_df, prev_time_rows)
    stan_fit_summary_prev_time_df <- stan_fit_summary_prev_time %>%
      dplyr::select(mean, '2.5%', '97.5%') %>%
      dplyr::mutate(week = row_number()) %>%
      dplyr::rename(CI_min = '2.5%', CI_max = '97.5%', mean_prevalence = mean)

    n_weeks <- nrow(stan_fit_summary_prev_time_df)

    # Extract sens/spec medians
    sens_median_rows <- grepl("^Se_median", rownames(fit_summary_df))
    stan_fit_summary_sens_median <- subset(fit_summary_df, sens_median_rows)
    stan_fit_summary_sens_median_weekly <- stan_fit_summary_sens_median[rep(1, n_weeks), ]
    stan_fit_summary_sens_median_weekly$week <- seq_len(n_weeks)

    spec_median_rows <- grepl("^Sp_median", rownames(fit_summary_df))
    stan_fit_summary_spec_median <- subset(fit_summary_df, spec_median_rows)
    stan_fit_summary_spec_median_weekly <- stan_fit_summary_spec_median[rep(1, n_weeks), ]
    stan_fit_summary_spec_median_weekly$week <- seq_len(n_weeks)

    # Extract values from list name
    name_parts <- strsplit(result_name, "_")[[1]]

    sim_stan_result_df <- data.frame(
      sim_sens = as.numeric(name_parts[2]),
      sim_spec = as.numeric(name_parts[4]),
      sim_prob = as.numeric(name_parts[6])
    )

    #Repeat df rows by number of weeks and add week no. col
    sim_stan_result_df <- sim_stan_result_df[rep(1, n_weeks), ]
    sim_stan_result_df$week <- stan_fit_summary_prev_time_df$week
    # Add prev col
    sim_stan_result_df$stan_prev <- stan_fit_summary_prev_time_df$mean_prevalence
    sim_stan_result_df$stan_prev_CI_low <- stan_fit_summary_prev_time_df$CI_min
    sim_stan_result_df$stan_prev_CI_high <- stan_fit_summary_prev_time_df$CI_max

    # add columns for each test's sensitivity and specificity
    for (i in 1:num_tests) {
      test_name <- paste0("test", i)

      # Sensitivity (Median)
      sim_stan_result_df[[paste0("stan_sens_", test_name, "_median")]] <- stan_fit_summary_sens_median_weekly$mean[i]
      sim_stan_result_df[[paste0("stan_sens_", test_name, "_median_CI_low")]] <- stan_fit_summary_sens_median_weekly$`2.5%`[i]
      sim_stan_result_df[[paste0("stan_sens_", test_name, "_median_CI_high")]] <- stan_fit_summary_sens_median_weekly$`97.5%`[i]
      # Specificity (Median)
      sim_stan_result_df[[paste0("stan_spec_", test_name, "_median")]] <- stan_fit_summary_spec_median_weekly$mean[i]
      sim_stan_result_df[[paste0("stan_spec_", test_name, "_median_CI_low")]] <- stan_fit_summary_spec_median_weekly$`2.5%`[i]
      sim_stan_result_df[[paste0("stan_spec_", test_name, "_median_CI_high")]] <- stan_fit_summary_spec_median_weekly$`97.5%`[i]

    }


    #Final results
    #Stan
    #pivot longer results for each test
    sim_stan_result_df_long <- sim_stan_result_df %>%
      pivot_longer(
        cols = matches("^stan_(sens|spec)_test\\d+_.*"),
        names_to = c("metric", "test_id", "stat"),
        names_pattern = "stan_(sens|spec)_(test\\d+)_(.*)"
      ) %>%
      unite("metric_stat", metric, stat, sep = "_") %>%
      pivot_wider(
        names_from = metric_stat,
        values_from = value
      ) %>%
      mutate(test_id = as.factor(test_id)) %>%
      select(sim_sens, sim_spec, sim_prob, test_id, everything())
    sim_stan_result_df_long$model_id <- result_name
    #Sim
    sim_data$test_parameters$model_id <- result_name #model ID

    sim_data$sim_data$model_id <- result_name #model ID

    # Define the R-estimate data frames
    R_dfs <- list(
      custom_R_estimate              = sim_data$custom_R_estimate,
      custom_R_estimate_true         = sim_data$custom_R_estimate_true,
      custom_growth_rate_estimate    = sim_data$custom_growth_rate_estimate,
      custom_growth_rate_estimate_true = sim_data$custom_growth_rate_estimate_true
      #epi_R_estimate                 = sim_data$epi_R_estimate$R, #have different data structures
      #epi_R_estimate_true            = sim_data$epi_R_estimate_true$R
    )

    # Ensure all have the model_id and standardise column types
    R_dfs <- purrr::imap(R_dfs, ~ .x %>%
                           dplyr::mutate(model_id = result_name) %>%
                           dplyr::mutate(dplyr::across(tidyselect::where(is.factor), as.character))  # avoid merge issues
    )

    # Define join keys
    key_cols <- c("model_id", "day_of_year")

    # Rename columns to suffix them with estimate source
    R_dfs_renamed <- purrr::imap(
      R_dfs,
      function(df, suffix) {
        df %>%
          dplyr::rename_with(
            .fn = ~ ifelse(.x %in% key_cols, .x, paste0(.x, "_", suffix))
          )
      }
    )

    # Join all on model_id + day_of_year
    R_est_df <- purrr::reduce(R_dfs_renamed, dplyr::full_join, by = key_cols)

    out <- list(
      sim_stan_result_df_ALL_long = sim_stan_result_df_long,
      sim_inputs_df_ALL = sim_data$test_parameters,
      sim_data_df_ALL = sim_data$sim_data,
      divergence_summary_df = divergence_check,
      R_est_df_ALL = R_est_df
    )

    return(out)

  }

  # Run param combos in parallel:
  required_packages <- c("LCtesterror", "future", "furrr",
                         "rstan", "dplyr", "tidyr", "gt", "magrittr", "ggplot2",
                         "stats", "graphics", "utils", "purrr", "deSolve")
  results_list <- furrr::future_map(
    seq_len(nrow(param_grid)),
    function(f) run.for.param.combo(param_grid[f, ]),
    .options = furrr_options(seed = TRUE, packages = required_packages)
  )

  # Combine results across all param combos:
  sim_stan_result_df_ALL_long <- do.call(rbind, lapply(results_list, `[[`, "sim_stan_result_df_ALL_long"))
  sim_inputs_df_ALL <- do.call(rbind, lapply(results_list, `[[`, "sim_inputs_df_ALL"))
  sim_data_df_ALL <- do.call(rbind, lapply(results_list, `[[`, "sim_data_df_ALL"))
  divergence_summary_df <- do.call(rbind, lapply(results_list, `[[`, "divergence_summary_df"))
  R_est_df_ALL <- do.call(rbind, lapply(results_list, `[[`, "R_est_df_ALL"))


  #final list of outputs to return:
  sim_stan_model_output <- list(sim_inputs = NULL, stan_results_df = NULL, divergence_summary = NULL, R_estimates = NULL)

  sim_stan_model_output$sim_inputs <- sim_inputs_df_ALL
  sim_stan_model_output$sim_data <- sim_data_df_ALL
  sim_stan_model_output$stan_results_df <- sim_stan_result_df_ALL_long
  sim_stan_model_output$divergence_summary <- divergence_summary_df
  sim_stan_model_output$R_estimates <- R_est_df_ALL

  return(sim_stan_model_output)

}


