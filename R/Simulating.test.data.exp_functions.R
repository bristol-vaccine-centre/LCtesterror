
#' @title Generates disease prevalence assuming exponential growth
#' @description Generates disease prevalence assuming exponential growth to model true prevalence at the start of an epidemic.
#'
#' @param days Number of days to simulate disease growth.
#' @param I0 Initial disease prevalence. Default = 0.0001.
#' @param beta0 Base transmission rate. Default = 0.06.
#' @param gamma Recovery rate. Default = 0.04.
#' @return A dataframe of model results at each timepoint. Includes: Time = day; I_t = prevalence; R0_t = reproduction number .
#'
#' @examples
#' if (interactive()) {
#'
#'exp_mod_result <- run.exponential.model(days=365, I0 = 0.0001,  beta0 = 0.06, gamma = 0.04)
#'
#'# Plot the simulated cases over time
#'exp_mod_result %>%
#' ggplot() +
#'  geom_line(aes(x=time, y=I_t)) +
#'  theme_classic() +
#'  labs(x="Days", y="Prevalence")
#'
#'# Plot the reproduction number over time
#'exp_mod_result %>%
#'  ggplot() +
#'   geom_line(aes(x=time, y=R0_t)) +
#'   theme_classic() +
#'   labs(x="Days", y="Reproduction number (Rt)") +
#'   scale_y_continuous(labels = function(x) sprintf("%.2f", x))
#'
#' }
#' @name run.exponential.model
#' @export run.exponential.model
run.exponential.model <- function(days = 365, I0 = 0.0001, beta0 = 0.06, gamma = 0.04) {

  t <- 1:days
  I_t <- numeric(length(t))  # Prev at each time point
  R0_t <- numeric(length(t)) # Repro number at each time point

  # initial prev / R0
  I_t[1] <- I0
  R0_t[1] <- beta0 / gamma

  # Loop over each time point to calculate prevalence and R0
  for (i in 2:length(t)) {
    beta_t <- beta0
    # growth rate r(t) = beta(t) - gamma
    r_t <- beta_t - gamma
    # Calculate prevalence at time t[i] using exponential growth
    I_t[i] <- I_t[i - 1] * exp(r_t)
    # Cap prevalence at 1
    if (I_t[i] > 1) {
      I_t[i] <- 1
    }
    # Calculate the reproduction number R0(t)
    R0_t[i] <- beta_t / gamma
  }

  results <- data.frame(time = t, I_t = I_t, R0_t = R0_t)

  return(results)
}


#' @title Simulates test data over time assuming exponential growth
#' @description Simulates test data over time assuming exponential growth
#' Returned overall sens/spec combines individual test sens/spec parameters assuming multiple parallel testing where any test positive is assumed to be a disease positive individual.
#' Also estimates the reproduction number from simulated prevalence data using two different methods.
#' Method 1: As the true prevalence data and test positivity data were simulated using an exponential growth model, R was calculated as R = 1 + (growth rate/recovery rate).
#' The growth rate (r) was calculated from the slope of the test positivity data by fitting a linear model to the log of the test positivity data over a 14-day sliding window. The recovery rate = gamma parameter value used in the exponential growth model.
#' This will produce the warning message: "In summary.lm(lm_fit) : essentially perfect fit: summary may be unreliable" because the model perfectly fits the data - this is expected and can be ignored.
#' Method 1 produces the R outputs 'custom_R_estimate' and 'custom_R_estimate_true' and the growth rate outputs 'custom_growth_rate_estimate' and 'custom_growth_rate_estimate_true'.
#' Method 2: Uses the EpiEstim package with method = "non_parametric_si". The si distribution is based on the mean generation interval 1/gamma.
#' Method 2 produces the R outputs 'epi_R_estimate' and 'epi_R_estimate_true' as returned from EpiEstim::estimate_R().
#'
#' @param sim_size Number of individuals to simulate test results for at each timepoint. Default = 1000.
#' @param days Number of days to simulate. Default = 365.
#' @param test_params Test parameters used to simulate test results along with true prevalence. Given as a list of lists containing sensitivity (sens =), specificity (spec =), and probability that the test is performed (p_performed = ), for each test.
#' Default = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1), test2 = list(sens = 0.99, spec = 0.99, p_performed = 1), test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8), test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8))
#' @param seed For set.seed(). Default = 953.
#' @param Est_R_window Size, in number of days, of sliding window for custom R(t) estimation. Default = 14.
#' @param Est_R_n_samples Number of samples for uncertainty estimation in custom R(t). Default = 1000.
#' @param mean_gi Mean generation interval used for estimating R(t) from EpiEstim. Default = 1/gamma.
#' @param max_t Maximum time (days) used for calculating the gi_distribution for estimating R(t) from EpiEstim. Default = 5 * mean_gi
#' @param params Exponential model parameters as a list. Defaults = list(I0 = 0.0001, beta0 = 0.06, gamma = 0.004).
#' @return A list containing:
#'  \describe{
#'   \item{test_parameters}{Test results table with row for each test (test_id) containing specified overall test parameters summarised across time (sens; spec; p_performed; disease_prev), simulated test_positivity, test_coverage (based on p_performed), and the estimated true disease prevalence estimate (disease_prev_est: based on the Rogan-Gladen equation) }
#'   \item{test_results}{Simulated binary test results for each individual (N=sim_size) with a column for day of year (i.e time in days - can be >365), based on true prev, sens and spec parameters.}
#'   \item{sim_data}{Simulated test data, including true prevalence, true infection status, test results, and simulated PPV and NPV values for each day simulated}
#'   \item{incidence_data}{Daily test-derived and true case counts, as well as overall sensitivity and specificity (combined across all tests assuming parallel testing with any test positive assumed a disease positive).}
#'   \item{custom_R_estimate}{Daily mean estimated R(t) based on observed (test-positive) incidence from simulated test data with test error using custom method. Contains overall sens/spec (overall across tests)}
#'   \item{custom_R_estimate_true}{Daily mean estimated R(t) based on true incidence using custom method.Contains overall sens/spec (overall across tests)}
#'   \item{custom_growth_rate_estimate}{Daily mean estimated growth rate and its SE from observed incidence from simulated test data with test error.Contains overall sens/spec (overall across tests)}
#'   \item{custom_growth_rate_estimate_true}{Daily mean estimated growth rate and SE from true incidence.Contains overall sens/spec (overall across tests)}
#'   \item{epi_R_estimate}{Daily mean R(t) estimates from EpiEstim using observed incidence from simulated test data with test error. Output from Epiestim, but contains overall sens/spec columns (overall across tests)}
#'   \item{epi_R_estimate_true}{Daily mean R(t) estimates from EpiEstim using true incidence. Output from Epiestim, but contains overall sens/spec columns (overall across tests)}
#' }
#' @importFrom EpiEstim estimate_R make_config
#' @importFrom dplyr mutate group_by summarise ungroup select starts_with
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats coef lm rnorm sd pexp
#' @name sim.test.data.exp
#' @examples
#' if (interactive()) {
#' sim_data <- sim.test.data.exp(sim_size = 10, days=365,
#' test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1),
#' test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
#' test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8),
#' test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)))
#'
#' #head(sim_data$test_results)
#' }

utils::globalVariables(c("test_id", "p_performed", "true_prev", "sens", "spec",
                         "pat_id", "test_result", "any_positive",
                         "CI", "test_positivity", "test_positivity_sim",
                         "day_of_year", "true_disease", "ppv", "npv",
                         "test_result_overall"))
#' @export sim.test.data.exp
sim.test.data.exp <- function(sim_size = 1000, days = 365,
                               test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1), test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
                                                  test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8), test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),
                               seed=953,
                               Est_R_window = 14, Est_R_n_samples = 1000, #for my R method
                               mean_gi = NULL, max_t = NULL, #For EpiEstim - Generation interval (mean and SD)
                               params =list(I0 = 0.0001, beta0 = 0.06, gamma = 0.004)) #for exp model
{

  set.seed(seed)

  # Define test characteristics
  tests <- tibble::tibble(
    test_id = names(test_params),
    sens = sapply(test_params, function(x) x$sens),
    spec = sapply(test_params, function(x) x$spec),
    p_performed = sapply(test_params, function(x) x$p_performed)
  )

  overall_sens <- 1 - prod(1 - tests$sens)
  overall_spec <- prod(tests$spec)

  if (is.null(mean_gi)) mean_gi <- 1 / params$gamma
  if (is.null(max_t)) max_t <- 5 * mean_gi

  # EXP model output
  exp_output <- run.exponential.model(days = days, I0 = params$I0, beta0 = params$beta0, gamma = params$gamma)


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

  # Simulate data using exponential model prevalence
  synth_data <- expand.grid(
    pat_id = 1:sim_size,
    day_of_year = 1:days
  ) %>%
    dplyr::mutate(
      true_disease = sapply(day_of_year, function(day) {
        rbinom(1, 1, exp_output$I_t[day])
      }),
      true_prev = sapply(day_of_year, function(day) {
        exp_output$I_t[day]
      })
    ) %>%
    dplyr::cross_join(tests) %>%
    dplyr::group_by(test_id) %>%
    dplyr::mutate(
      test_result = ifelse(
        rbinom(n(), 1, p_performed) == 0,
        NA,
        true_disease * rbinom(n(), 1, sens) + (1 - true_disease) * (1 - rbinom(n(), 1, spec))
      ),
      ppv = calculate_ppv(sens, spec, true_prev),
      npv = calculate_npv(sens, spec, true_prev),
      sens = sens,
      spec = spec
    )

  # wide format
  wide_synth_data <- synth_data %>%
    dplyr::select(pat_id, day_of_year, true_prev, true_disease, test_id, test_result, ppv, npv, sens, spec) %>%
    tidyr::pivot_wider(names_from = test_id, values_from = c(test_result, ppv, npv, sens, spec))

  test_results <- wide_synth_data %>%
    dplyr::select(dplyr::starts_with("test_result"), day_of_year)

  # Check simulation is doing the right thing using Rogan-Gladen (for overall params):
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
      disease_prev_sim = mean(true_prev),
      test_positivity_sim = mean(test_result, na.rm = TRUE),
      test_coverage_sim = mean(!is.na(test_result)),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(disease_prev_RG_est = rogan_gladen(test_positivity_sim, sens, spec)) %>% #estimates true disease prev using RG
    # add overall test positivity into df
    dplyr::left_join(overall_test_positivity_summary, by = c("sens", "spec"))


  incidence_data <- synth_data %>%
    group_by(pat_id, day_of_year) %>%
    mutate(test_result_overall = as.numeric(case_when(any(test_result == 1, na.rm = TRUE) ~ 1, TRUE ~ 0))) %>%
    ungroup() %>%
    group_by(day_of_year) %>%
    summarise(case_count = sum(test_result_overall, na.rm = TRUE), case_count_true = sum(true_disease, na.rm=TRUE)) %>%
    mutate(overall_sens = overall_sens, overall_spec = overall_spec)


  # Calculation of R0(t) / Re(t) from daily observed cases (no adjustment for sens/spec)

  my.estimate.R <- function(incidence_data, gamma, day_of_year, window_size = Est_R_window , n_samples = Est_R_n_samples) {

    R_results <- list(R_t_df = NULL, growth_rate_df = NULL)

    R_t_df <- data.frame(day_of_year = day_of_year, R_t_mean = NA, R_t_sd = NA)
    growth_rate_df <- data.frame(day_of_year = day_of_year, growth_rate = NA, SE = NA)


    # Loops over each time point and calculates R_t using a sliding window
    for (i in seq(window_size, length(incidence_data))) {
      # subsets incidence data over window
      cases_window <- incidence_data[(i - window_size + 1):i]
      # calculates growth rate r as slope of log(cases) over the window
      log_cases <- log(cases_window + 1)  # (avoid log(0) by adding 1)
      lm_fit <- lm(log_cases ~ seq_along(log_cases))

      # Extract growth rate and se
      r <- stats::coef(lm_fit)[2]
      r_se <- summary(lm_fit)$coefficients[2, 2]  # SE of slope

      # R_t
      R_t <- 1 + r / gamma

      # Draw samples for R_t from normal distribution with mean R_t and SD based on SE
      R_t_samples <- rnorm(n_samples, mean = R_t, sd = r_se / gamma)

      # save mean and SD of R_t
      R_t_df$R_t_mean[i] <- mean(R_t_samples)
      R_t_df$R_t_sd[i] <- stats::sd(R_t_samples)
      # saves growth rate and SE
      growth_rate_df$growth_rate[i] <- r
      growth_rate_df$SE[i] <- r_se
    }


    R_results$R_t_df <- R_t_df
    R_results$growth_rate_df <- growth_rate_df

    return(R_results)
  }

  # Calculate R_t estimates using incidence_data for overall cases and true cases
  my_R_t_estimates <- my.estimate.R(incidence_data$case_count, params$gamma, incidence_data$day_of_year)
  my_R_t_estimates_true <- my.estimate.R(incidence_data$case_count_true, params$gamma, incidence_data$day_of_year)



  #Calculate using EpiEstim package

  # discrete approximation of exponential distribution
  gi_distribution <- sapply(0:max_t, function(t) {
    stats::pexp(t + 1, rate = 1 / mean_gi) - stats::pexp(t, rate = 1 / mean_gi)
  })
  # Ensure si_distr[1] = 0 and normalise the distribution
  gi_distribution[1] <- 0
  gi_distribution <- gi_distribution / sum(gi_distribution)


  #setting rolling time window (default = weekly sliding)
  #Every 2 weeks:
  t_start <- seq(2, length(incidence_data$case_count) - 13)
  t_end <- t_start + 13

  epi_R_estimate <- EpiEstim::estimate_R(incidence_data$case_count,
                                         method = "non_parametric_si",
                                         config = make_config(list(si_distr = gi_distribution, t_start = t_start,
                                                                   t_end = t_end)))

  epi_R_estimate$R <- epi_R_estimate$R %>%
    mutate(overall_sens = as.numeric(overall_sens), overall_spec = as.numeric(overall_spec))


  #true R
  epi_R_estimate_true <- EpiEstim::estimate_R(incidence_data$case_count_true,
                                              method = "non_parametric_si",
                                              config = make_config(list(si_distr = gi_distribution, t_start = t_start,
                                                                        t_end = t_end)))


  epi_R_estimate_true$R <- epi_R_estimate_true$R %>%
    mutate(overall_sens = as.numeric(overall_sens), overall_spec= as.numeric(overall_spec))


  #save results
  results <- list()

  results$test_results <- test_results
  results$test_parameters <- reconstruct
  results$sim_data <- wide_synth_data
  results$incidence_data <- incidence_data
  results$exp_params <- list(I0 = params$I0, beta0 = params$beta0, gamma=params$gamma,
                            mean_gi=mean_gi, max_t=max_t)
  results$custom_R_estimate <- my_R_t_estimates$R_t_df %>%
    mutate(overall_sens = as.numeric(overall_sens), overall_spec = as.numeric(overall_spec))
  results$custom_R_estimate_true <- my_R_t_estimates_true$R_t_df %>%
    mutate(overall_sens = as.numeric(overall_sens), overall_spec = as.numeric(overall_spec))

  results$custom_growth_rate_estimate <- my_R_t_estimates$growth_rate_df %>%
    mutate(overall_sens = as.numeric(overall_sens), overall_spec = as.numeric(overall_spec))
  results$custom_growth_rate_estimate_true <- my_R_t_estimates_true$growth_rate_df %>%
    mutate(overall_sens = as.numeric(overall_sens), overall_spec = as.numeric(overall_spec))

  results$epi_R_estimate <- epi_R_estimate
  results$epi_R_estimate_true <- epi_R_estimate_true

  return(results)
}













