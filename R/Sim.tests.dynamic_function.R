
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
#' @param sequential_testing Logical. Are any tests only performed if the previous test / parallel set of tests was positive? Assumes the same sample is used for testing (i.e no time lag). FALSE means all tests are performed in parallel. Default = FALSE.
#' @param seq_order Must be specified if sequential_testing is TRUE. A list of test ID's (to match those in test_params) with a number specifying the order in which the test occurs if testing sequentially conditional on a positive test result.
#' Multiple tests can be performed in parallel at each testing stage and therefore will have the same seq_order number. E.g list(test1=1, test2=1, test3=2, test4=2) Default = NULL.
#' @param delay Logical indicating whether to simulate effects of 'delay until testing' on sensitivity. Default = FALSE.
#' @param delay_distribution Optional specified distribution of delay effect. A function that takes an integer `n` and returns a vector of length `n` representing individual-specific delays. Default = sample(0:14, n, replace = TRUE).
#' @param delay_effect_fn Optional function to simulate effects of 'delay until testing' on sensitivity. Default = function(delay_day, sens) pmax(sens - 0.02 * delay_day, 0.5)
#' @return A table of test parameters (specified and simulated for comparison) and a table containing binary test results for each individual
#'   \describe{
#'   \item{test_parameters}{Test results table with row for each test (test_id) containing specified test parameters (sens; spec; p_performed; disease_prev), simulated test_positivity, test_coverage (based on p_performed), and the estimated true disease prevalence estimate (disease_prev_est: based on the Rogan-Gladen equation) }
#'   \item{test_results}{Simulated binary test results for each individual (N=sim_size), based on specified true prev, sens and spec parameters.
#'         If delay = TRUE, additional delay columns are included specifying the delay for each test.}
#' }
#' @importFrom dplyr filter group_by ungroup summarise mutate n cross_join select case_when rowwise left_join arrange
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
                         "test_positivity_sim", "sens_sim", "test_result_updated"))

default_delay_distribution <- function(n) sample(0:14, n, replace = TRUE)

default_delay_effect <- function(delay_day, sens) {
  pmax(sens - 0.02 * delay_day, 0.5)  # Ensures sens doesn't go below 0.5
}

#' @export sim.test.data
sim.test.data <- function(disease_prev = 0.2, sim_size = 1000,
                          test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1), test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
                                             test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8), test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),
                          seed = 953,
                          sequential_testing = FALSE, seq_order = NULL,
                          delay = FALSE, delay_distribution = default_delay_distribution,
                          delay_effect_fn = default_delay_effect
                          ) {

  #** need to add overall_sens / spec and ppv/npv to output as in the time simu functions
  #(and think about how the delay effect is included)

  set.seed(seed)

  tests <- tibble::tibble(
    test_id = names(test_params),
    sens = sapply(test_params, function(x) x$sens),
    spec = sapply(test_params, function(x) x$spec),
    p_performed = sapply(test_params, function(x) x$p_performed),
    seq_order = 1 #default all tests done in parallel
  )

  max_seq <- max(tests$seq_order)


  #if sequential testing, update seq_order with parallel testing groups:
  if(sequential_testing) {

    #  Check seq_order correctly specified
    if (is.null(seq_order)) {
      stop("Error: 'seq_order' must be provided when sequential_testing = TRUE.")
    }
    if (!setequal(names(seq_order), names(test_params))) {
      stop("Error: Names in 'seq_order' must match test IDs in 'test_params'.")
    }

    seq_map <- tibble::tibble(
      test_id = names(seq_order),
      seq_order = unlist(seq_order)
    )

    tests <- tests %>%
      dplyr::select(-seq_order) %>%
      dplyr::left_join(seq_map, by=c("test_id"))

    #defining overall sens/spec of tests for sequenial or parallel testing
    overall_sens_parallel <- c()
    overall_spec_parallel <- c()

    for(i in 1:max_seq) {

      tests_parallel <- dplyr::filter(tests, seq_order == i)

      overall_sens_parallel[i] <- 1 - prod(1 - tests_parallel$sens)
      overall_spec_parallel[i] <- prod(tests_parallel$spec)

    }

    overall_sens <- prod(overall_sens_parallel)
    overall_spec <- 1 - prod(1 - overall_spec_parallel)

  } else {

    overall_sens <- 1 - prod(1 - tests$sens)
    overall_spec <- prod(tests$spec)

  }



  # Function to calculate PPV
  calculate_ppv <- function(sens, spec, prev) {
    ppv <- (sens * prev) / ((sens * prev) + ((1 - spec) * (1 - prev)))
    return(ppv)
  }

  # Function to calculate NPV
  calculate_npv <- function(sens, spec, prev) {
    npv <- (spec * (1 - prev)) / (((1 - sens) * prev) + (spec * (1 - prev)))
    return(npv)
  }


  #simulate data:
  synth_data <- tibble::tibble(
    pat_id = 1:sim_size,
    true_prev = sample(c(rep(1, sim_size * disease_prev), rep(0, sim_size * (1 - disease_prev))),
                       size = sim_size, replace = TRUE)
  ) %>%
    dplyr::cross_join(tests) %>%
    dplyr::mutate(delay_days = if (delay == TRUE) delay_distribution(n()) else 0
    ) %>%
    dplyr::group_by(test_id) %>%
    #adj_sens = sens after adjustment for delay effect (even if 0):
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


  #if sequential testing, then remove subsequent test results if no results in parallel testing group (seq_order) are positive:
  if(sequential_testing) {

    synth_data <- synth_data %>%
      dplyr::ungroup() %>%
      dplyr::arrange(pat_id, seq_order, test_id) %>%
      dplyr::group_by(pat_id) %>%
      # propagate NAs down seq_order
      dplyr::mutate(test_result_updated = {
        seq_orders <- unique(seq_order)
        result <- test_result
        for (i in seq_along(seq_orders)) {
          current_seq <- seq_orders[i]
          if (i == 1) next  # keep seq_order 1 as is
          prev_seq <- seq_orders[i - 1]
          # if all previous seq_order results are NA/0, set current to NA
          if (all(result[seq_order == prev_seq] %in% c(NA, 0))) { # result[seq_order == prev_seq] selects all test results from the previous step (group of parallel tests).
            result[seq_order == current_seq] <- NA_real_ #result[seq_order == current_seq] selects all rows in the current step
          }
        }
        result  }
      ) %>%
      ungroup() %>%
      dplyr::mutate(test_result = test_result_updated) %>%
      dplyr::select(-test_result_updated)
  }


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

  # Check simulation is doing the right thing using Rogan-Gladen:
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

