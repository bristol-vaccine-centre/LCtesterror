
#' @title Simulates test data based on true prevalence, sensitivity and specificity.
#' @description Simulates test data for individuals based on true prevalence, and test sensitivity and specificity.
#' Can simulate multiple test results for each individual with different test sens/spec parameters and different probabilities that each test is performed.
#' The simulated data does not take into account dependencies between tests.
#'
#' @param disease_prev True population disease prevalence to simulate (between 0-1). Default = 0.2.
#' @param sim_size Number of individuals to simulate test results for. Default = 1000.
#' @param test_params Test parameters used to simulate test results along with true prevalence. Given as a list of lists for each test containing sensitivity (sens =), specificity (spec =), and probability that the test is performed (p_performed = ).
#' Default = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1), test2 = list(sens = 0.99, spec = 0.99, p_performed = 1), test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8), test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),
#' @param seed for set.seed(). Default = 953.
#' @param delay Logical indicating whether to simulate effects of 'delay until testing' on sensitivity. Default = FALSE.
#' @param delay_distribution Optional specified distribution of delay effect. A function that takes an integer `n` and returns a vector of length `n` representing individual-specific delays. Default = sample(0:14, n, replace = TRUE).
#' @param delay_effect_fn Optional function to simulate effects of 'delay until testing' on sensitivity. Default = function(delay_day, sens) pmax(sens - 0.02 * delay_day, 0.5)
#' @return A table of test parameters (specified and simulated for comparison) and a table containing binary test results for each individual
#'   \describe{
#'   \item{test_parameters}{Test results table with row for each test (test_id) containing specified test parameters (sens; spec; p_performed; disease_prev), simulated test_positivity, test_coverage (based on p_performed), and the estimated true disease prevalence estimate (disease_prev_est: based on the Rogan-Gladen equation) }
#'   \item{test_results}{Simulated binary test results for each individual (N=sim_size), based on specified true prev, sens and spec parameters.
#'         If delay = TRUE, additional delay columns are included specifying the delay for each test.}
#' }
#' @export
#' @importFrom dplyr filter group_by summarise mutate n cross_join select case_when rowwise left_join
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom stats rbinom na.omit quantile
#' @importFrom utils globalVariables
#' @importFrom tidyselect starts_with
#' @name sim.test.data
#' @examples
#' if (interactive()) {
#' sim_data <- sim.test.data(disease_prev = 0.2, sim_size = 100,
#' test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1),
#' test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
#' test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8),
#' test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)))
#'
#' #head(sim_data$test_results)
#'
#' # Example with 'delay until test' simulated
#' sim_data <- sim.test.data(disease_prev = 0.2, sim_size = 100,
#' test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1),
#' test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
#' test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8),
#' test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),
#' delay = TRUE)
#' }
#'


utils::globalVariables(c("test_id", "p_performed", "true_prev", "sens", "spec",
                         "pat_id", "test_result", "any_positive",
                         "CI", "test_positivity", "delay_days", "adj_sens",
                         "test_positivity_sim", "sens_sim"))

default_delay_distribution <- function(n) sample(0:14, n, replace = TRUE)

default_delay_effect <- function(delay_day, sens) {
  pmax(sens - 0.02 * delay_day, 0.5)  # Ensures sens doesn't go below 0.5
}


sim.test.data <- function(disease_prev = 0.2, sim_size = 1000,
                          test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1), test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
                                             test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8), test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),
                          seed = 953,
                          delay = FALSE, delay_distribution = default_delay_distribution,
                          delay_effect_fn = default_delay_effect
                          ) {

  set.seed(seed)

  tests <- tibble::tibble(
    test_id = names(test_params),
    sens = sapply(test_params, function(x) x$sens),
    spec = sapply(test_params, function(x) x$spec),
    p_performed = sapply(test_params, function(x) x$p_performed)
  )

  synth_data <- tibble::tibble(
    pat_id = 1:sim_size,
    true_prev = sample(c(rep(1, sim_size * disease_prev), rep(0, sim_size * (1 - disease_prev))),
                       size = sim_size, replace = TRUE)
  ) %>%
    dplyr::cross_join(tests) %>%
    dplyr::mutate(delay_days = if (delay == TRUE) delay_distribution(n()) else 0
    ) %>%
    dplyr::group_by(test_id) %>%
    dplyr::mutate(adj_sens = if (delay == TRUE) delay_effect_fn(delay_days, sens) else sens
      ) %>%
    dplyr::mutate(
      test_result = ifelse(
        stats::rbinom(n(), 1, p_performed) == 0,
        NA, # set test results which were not performed to NA
        # otherwise get a result that is consistent with test sens and spec
        #i.e binary test result based on proportion of population that will test positive with these parameters
        true_prev * stats::rbinom(n(), 1, adj_sens) + (1 - true_prev) * (1 - stats::rbinom(n(), 1, spec))
        # test delay included here
      )
    )

  # data set with a test per column (1-pos, 0-neg, NA-not done)
  # Pivot wider
  select_vars <- c("pat_id", "true_prev", "test_id", "test_result")
  if (delay == TRUE) { select_vars <- c(select_vars, "delay_days") }
  wide_synth_data1 <- synth_data %>%
    dplyr::select(all_of(select_vars))

  if (delay == TRUE) {
    wide_synth_data1 <- wide_synth_data1 %>%
      tidyr::pivot_wider(
        names_from = test_id,
        values_from = c(test_result, delay_days),
        names_glue = "{.value}_{test_id}"
      ) %>%
      dplyr::rename_with(
        .cols = tidyselect::starts_with("test_result_"),
        .fn = ~ gsub("^test_result_", "", .x)
      ) %>%
      dplyr::rename_with(
        .cols = tidyselect::starts_with("delay_days_"),
        .fn = ~ gsub("^delay_days_", "delay_", .x)
      )
  } else {
    wide_synth_data1 <- wide_synth_data1 %>%
      tidyr::pivot_wider(
        names_from = test_id,
        values_from = test_result
      )
  }

test_result_cols <- setdiff(names(wide_synth_data1), c("pat_id", "true_prev", grep("^delay_", names(wide_synth_data1), value = TRUE)))
  wide_synth_data <- wide_synth_data1 %>%
    dplyr::filter(rowSums(!is.na(dplyr::across(all_of(test_result_cols)))) > 0) %>% #remove individuals with no test results
    dplyr::select(-pat_id, -true_prev)

  # Check simulation is doing the right thing using RG:
  rogan_gladen <- function(ap, sens, spec) {
    dplyr::case_when(
      ap <= 1 - spec ~ 0,
      sens <= ap ~ 1,
      TRUE ~ (ap + spec - 1) / (sens + spec - 1)
    )
  }

  # Calculate an individuals overall test result
  individual_positivity <- synth_data %>%
    dplyr::group_by(pat_id, sens, spec) %>%
    dplyr::summarise(any_positive = as.numeric(any(test_result == 1, na.rm = TRUE)), .groups = 'drop')
  # calculate overall test positivity for each parameter combination
  overall_test_positivity_summary <- individual_positivity %>%
    dplyr::group_by(sens, spec) %>%
    dplyr::summarise(overall_test_positivity_sim = mean(any_positive), .groups = 'drop')
  # Bootstrap to calculate confidence intervals (95%)
  bootstrap_CI <- function(data, n_bootstrap = 1000) {
    boot_results <- replicate(n_bootstrap, {
      sample_data <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
      mean(sample_data$any_positive)
    })
    c(stats::quantile(boot_results, 0.025), stats::quantile(boot_results, 0.975))
  }
  # Apply bootstrap CI to overall test positivity
  overall_test_positivity_summary <- overall_test_positivity_summary %>%
    dplyr::rowwise() %>%
    dplyr::mutate(CI = list(bootstrap_CI(individual_positivity[individual_positivity$sens == sens & individual_positivity$spec == spec, ]))) %>%
    dplyr::mutate(CI_lower = CI[1], CI_upper = CI[2]) %>%
    dplyr::select(-CI)
  # can reconstruct the input parameters from the data?
  reconstruct <- synth_data %>%
    # the input parameters
    dplyr::group_by(test_id, sens, spec, p_performed) %>%
    # equivalents from the data
    dplyr::summarise(
      sens_sim = mean(adj_sens),
      disease_prev_sim = mean(true_prev),
      test_positivity_sim = mean(test_result, na.rm = TRUE),
      test_coverage_sim = mean(!is.na(test_result)),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(disease_prev_RG_est = rogan_gladen(test_positivity_sim, sens_sim, spec)) %>% #estimates true disease prev using RG
    # add overall test positivity into df
    dplyr::left_join(overall_test_positivity_summary, by = c("sens", "spec"))

  sim_results <- list()
  sim_results$test_parameters <- reconstruct
  sim_results$test_results <- wide_synth_data

  return(sim_results)
}

