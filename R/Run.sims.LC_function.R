
#' @title Runs multiple latent class model simulations with different parameter values.
#' @description Runs multiple simulations of a bayesian LC stan model to infer true disease prevalence, test sensitivity and specificity from
#' simulated test data with different parameters values (prevalence, sensitivity, specificity, and probability of the test occurring).
#' Can simulate multiple test results but all with identical parameters.
#'
#' @param num_tests A numeric value for the number of tests to simulate and model.
#' @param prev_vec Values of disease prevalence to simulate as a numeric vector for each simulation. Default = `c(0.2)`.
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
#' @param set.seed Random seed for set.seed(). Default = 9876.
#' @return Stan model fit and various summary outputs:
#' @return A list containing:
#' \describe{
#'   \item{sim_inputs}{A data frame of all simulation input parameters per model run.
#'   Output from sim.test.data()$test_parameters}
#'   \item{stan_results_df}{A data frame of posterior estimates and 95% credible intervals
#'     for sensitivity, specificity, and prevalence for each simulation.}
#'   \item{divergence_summary}{A data frame summarising any divergent transitions or low ESS in the Stan models.
#'   Output from check.divergent.transitions().}
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
#' @name run.sims.LC
#' @export
#' @examples
#' if (interactive()) {
#'  results <- run.sims.LC(num_tests = 4,
#'                        prev_vec = c(0.1, 0.2),
#'                        spec_vec = c(1),
#'                        sens_vec = c(1),
#'                        p_performed_vec= c(1),
#'                        sim_size = 100,
#'                        chains = 2,
#'                        data_ID = "sims")
#' }
#'
#'
#Function for running lots of LC runs using sim data with different parameters

utils::globalVariables(c("ess_threshold", "matches", "unite", "metric", "stat", "metric_stat",
  "sim_sens", "sim_spec", "sim_prev", "sim_prob", "everything"))

run.sims.LC <- function(num_tests, prev_vec= c(0.2), spec_vec= c(1), sens_vec= c(1), p_performed_vec= c(1),
                        sim_size=1000, iter=1000, chains=4, warmup=500, stan_arg=list(), data_ID = "sims",
                        prior_spec = NULL, prior_sens = NULL, set.seed = 9876
                        ) {

  set.seed(set.seed)

  sim_stan_result_df_ALL <- data.frame()
  sim_inputs_df_ALL <- data.frame()
  divergence_summary_df <- data.frame()


   # Run model for each parameter combo:
  for (s in 1:length(sens_vec)) {
    for (c in 1:length(spec_vec)) {
      for (p in 1:length(prev_vec)) {
        for (r in 1:length(p_performed_vec)) {

          test_params <- list()

          for (i in 1:num_tests) {
            testname <- paste0("test", i)
            test_params[[testname]] <- list(sens = sens_vec[s], spec = spec_vec[c], p_performed = p_performed_vec[r])
          }


        # Simulate data
        sim_data <- sim.test.data(disease_prev=prev_vec[p], sim_size=sim_size, test_params=test_params)

        result_name <- paste("sens_", sens_vec[s], "_spec_", spec_vec[c], "_prev_", prev_vec[p], "_p_", p_performed_vec[r], sep = "")

        # Run model
        result <- run.LC.model(data=sim_data$test_results, num_tests=num_tests,
                               iter=iter, chains=chains, warmup=warmup, stan_arg=stan_arg,
                               data_ID = data_ID, prior_spec = prior_spec, prior_sens = prior_sens)

        stan_fit <- result$stan_fit
        fit_summary <- rstan::summary(stan_fit)
        stan_fit_summary <- fit_summary$summary
        fit_summary_df <- as.data.frame(stan_fit_summary)


        # Function to check model for divergent transitions
        divergence_check <- check.divergent.transitions(stan_fit, result_name)
        divergence_summary_df <- rbind(divergence_summary_df, divergence_check)


        sim_stan_result_df <- data.frame()

        # Extract prevalence with CIs
        prev_rows <- grepl("^prev$", rownames(fit_summary_df))
        stan_fit_summary_prev <- subset(fit_summary_df, prev_rows)
        stan_prev <- stan_fit_summary_prev$mean
        stan_prev_CI_low <- stan_fit_summary_prev$`2.5%`
        stan_prev_CI_high <- stan_fit_summary_prev$`97.5%`

        # Extract sens/spec medians
        sens_median_rows <- grepl("^Se_median", rownames(fit_summary_df))
        stan_fit_summary_sens_median <- subset(fit_summary_df, sens_median_rows)

        spec_median_rows <- grepl("^Sp_median", rownames(fit_summary_df))
        stan_fit_summary_spec_median <- subset(fit_summary_df, spec_median_rows)

        # Extract values from list name
        name_parts <- strsplit(result_name, "_")[[1]]

        sim_stan_result_df <- data.frame(
          sim_sens = as.numeric(name_parts[2]),
          sim_spec = as.numeric(name_parts[4]),
          sim_prev = as.numeric(name_parts[6]),
          sim_prob = as.numeric(name_parts[8])
        )

        # add columns for each test's sensitivity and specificity
        for (i in 1:num_tests) {
          test_name <- paste0("test", i)

          # Sensitivity (Median)
          sim_stan_result_df[[paste0("stan_sens_", test_name, "_median")]] <- stan_fit_summary_sens_median$mean[i]
          sim_stan_result_df[[paste0("stan_sens_", test_name, "_median_CI_low")]] <- stan_fit_summary_sens_median$`2.5%`[i]
          sim_stan_result_df[[paste0("stan_sens_", test_name, "_median_CI_high")]] <- stan_fit_summary_sens_median$`97.5%`[i]
          # Specificity (Median)
          sim_stan_result_df[[paste0("stan_spec_", test_name, "_median")]] <- stan_fit_summary_spec_median$mean[i]
          sim_stan_result_df[[paste0("stan_spec_", test_name, "_median_CI_low")]] <- stan_fit_summary_spec_median$`2.5%`[i]
          sim_stan_result_df[[paste0("stan_spec_", test_name, "_median_CI_high")]] <- stan_fit_summary_spec_median$`97.5%`[i]

        }

        # Add prev
        sim_stan_result_df$stan_prev <- stan_prev
        sim_stan_result_df$stan_prev_CI_low <- stan_prev_CI_low
        sim_stan_result_df$stan_prev_CI_high <- stan_prev_CI_high


        #Final results
        #Stan
        sim_stan_result_df_ALL <- rbind(sim_stan_result_df_ALL, sim_stan_result_df)
        #pivot longer results for each test
        sim_stan_result_df_ALL_long <- sim_stan_result_df_ALL %>%
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
          select(sim_sens, sim_spec, sim_prev, sim_prob, test_id, everything())
        sim_stan_result_df_ALL_long$model_id <- result_name
        #Sim
        sim_data$test_parameters$model_id <- result_name #model ID
        sim_inputs_df_ALL <- rbind(sim_inputs_df_ALL, sim_data$test_parameters)
        #Model
        #stan_model <- result$stan_fit #not saved b/c high memory cost


        #clear environment after each model run, but keeps neccessary bits
        required_vars <- c("num_tests", "sens_vec", "spec_vec", "prev_vec", "p_performed_vec",
                           "sim_size", "iter", "chains", "warmup", "s", "c", "p", "r",
                           "data_ID", "prior_spec", "prior_sens", "set.seed", "stan_arg",
                           "prior_spec", "prior_sens",
                           "sim_inputs_df_ALL", "sim_stan_result_df_ALL", "sim_stan_result_df_ALL_long", "divergence_summary_df")
        rm(list = setdiff(ls(), required_vars))
        gc()


        }
      }
    }
  }

  sim_stan_model_output <- list(sim_inputs = NULL, stan_results_df = NULL)

  sim_stan_model_output$sim_inputs <- sim_inputs_df_ALL
  sim_stan_model_output$stan_results_df <- sim_stan_result_df_ALL_long
  sim_stan_model_output$divergence_summary <- divergence_summary_df

  return(sim_stan_model_output)

}



