
#These plotting functions do not currently account for different pre-test probabilities (sim_prob) values for the same simulated sens/spec combinations

#' @title Plots LC model simulation results.
#' @description Plots LC model simulation results generated from run.sims.LC() to evaluate model performance at different parameter values.
#' Compares LC model inferred prevalence, test sensitivities and test specificities against simulated parameters.
#' Will plot across any simulated any value of the probability of the test occurring (sim_prob),
#' so if multiple p_performed were simulated in run.sims.LC(), these will need to be pre-filtered to plot each sim_prob paramter value one at a time.
#'
#' @param sim_stan_results Named output from run.sims.LC().
#'
#' @return Plots of inferred parameters against simulated parameters. Returns 3 plots:
#' #' \describe{
#'   \item{1) Prevalence}{LC model inferred prevalence (mean across chains, mean across iterations)) vs simulated prevalence;}
#'   \item{2) Sensitivity}{LC model inferred test sensitivities (mean across chains, median across iterations) vs simulated test sensitivities;}
#'   \item{3) Specificity}{LC model inferred test specificities (mean across chains, median across iterations) vs simulated test specificities;}
#' }
#' @importFrom ggplot2 ggplot theme theme_classic theme_linedraw geom_abline geom_point geom_errorbar facet_grid labs coord_cartesian element_text element_rect element_blank unit label_both scale_color_manual scale_color_brewer
#' @importFrom dplyr rename left_join filter
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom utils globalVariables
#' @name sim.result.plot
#' @examples
#' if (interactive()) {
#' sim_results <- run.sims.LC(num_tests = 4,
#'                        prev_vec = c(0.1, 0.2),
#'                        spec_vec = c(1, 0.98),
#'                        sens_vec = c(1),
#'                        p_performed_vec= c(1),
#'                        sim_size = 100,
#'                        chains = 2,
#'                        data_ID = "sims")
#'
#' sim.result.plot(sim_results)
#' }
#'

utils::globalVariables(c("Simulated_prev", "Simulated_sens", "Simulated_spec", "stan_prev", "stan_prev_CI_low",
                         "stan_prev_CI_high", "sens_median", "Test_ID", "sens_median_CI_low", "sens_median_CI_high",
                         "spec_median", "spec_median_CI_low", "spec_median_CI_high", "overall_test_positivity_sim",
                         "Simulated_overall_test_positivity", "CI_lower", "CI_upper"))
#' @export
sim.result.plot <- function(sim_stan_results) {

  plots <- list()


  sim_stan_results_combined <- dplyr::left_join(sim_stan_results$stan_results_df, sim_stan_results$sim_inputs, by = c("model_id", "test_id"))


  sim_stan_results_combined <-  sim_stan_results_combined %>%
    dplyr::rename(Simulated_prev = sim_prev,
                  Simulated_sens = sim_sens,
                  Simulated_spec = sim_spec,
                  Test_ID = test_id,
                  Simulated_overall_test_positivity = overall_test_positivity_sim)

  theme <- theme_linedraw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              strip.background = element_rect(fill = "white"),
                              strip.text = element_text(size = 11, colour = "black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.key = element_rect(fill = "white", colour = NA),
                              legend.key.size = unit(1.2, "lines"),
                              legend.position = "bottom",
                              legend.text = element_text(size=10), legend.title = element_text(size=12),
                              axis.title = element_text(size = 14))


  plot_prev <- sim_stan_results_combined %>%
    dplyr::filter(Test_ID == "test1") %>% #all tests have the same results
    ggplot(aes(x = Simulated_prev)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(aes(y = stan_prev, colour="Inferred prevalence"), shape=21, fill="black", size=1) +
    geom_errorbar(aes(ymin = stan_prev_CI_low, ymax = stan_prev_CI_high, colour="Inferred prevalence"), width = 0, linewidth=1, alpha=0.8) +
    #overall test positivity (assuming any test positive is a positive):
    geom_point(aes(y = Simulated_overall_test_positivity, color = "Test positivity"), shape=21, fill="tomato2", size=0.8, alpha=0.8) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = "Test positivity"), width = 0, linewidth=0.8, alpha=0.6) +
    facet_grid(Simulated_spec ~ Simulated_sens, scales = "free",
               labeller = label_both) +
  scale_color_manual(name = "Estimate type",
      values = c("Inferred prevalence" = "black", "Test positivity" = "tomato2")) +
    labs(title = "Simuated vs inferred prevalence", x = "Simulated prevalence", y = "Inferred prevalence") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$prev <- plot_prev

  plot_sens_median <- sim_stan_results_combined %>%
    ggplot(aes(x = Simulated_sens, y = sens_median, colour=Test_ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_median_CI_low, ymax = sens_median_CI_high), width = 0, linewidth=1,alpha=0.4) +
    facet_grid(Simulated_spec ~ Simulated_prev, scales = "free",
               labeller = label_both) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "Simuated vs inferred sensitivity", x = "Simulated sensitivity", y = "Inferred sensitivity") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_median <-  plot_sens_median

  plot_spec_median <-  sim_stan_results_combined %>%
    ggplot(aes(x = Simulated_spec, y = spec_median, colour=Test_ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_median_CI_low, ymax = spec_median_CI_high), width = 0, linewidth=1, alpha=0.4) +
    facet_grid(Simulated_sens ~ Simulated_prev, scales = "free",
               labeller = label_both) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "Simuated vs inferred specificity", x = "Simulated specificity", y = "Inferred specificity") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$spec_median <- plot_spec_median


  return(plots)
}




#' @title Plots LC model simulation results for model with time.
#' @description Plots LC model simulation results generated from run.sims.LC() to evaluate model performance at different parameter values.
#' Compares LC model inferred prevalence, test sensitivities and test specificities against simulated parameters.
#' Will plot across any simulated any value of the probability of the test occurring (sim_prob),
#' so if multiple p_performed were simulated in run.sims.LC(), these will need to be pre-filtered to plot each sim_prob paramter value one at a time.
#'
#' @param sim_stan_results Named output from run.sims.LC.time().
#'
#' @return Plots of inferred parameters against simulated parameters. Returns 3 plots:
#' #' \describe{
#'   \item{1) Prevalence}{LC model inferred prevalence (mean across chains, mean across iterations)) vs simulated prevalence;}
#'   \item{2) Sensitivity}{LC model inferred test sensitivities (mean across chains, median across iterations) vs simulated test sensitivities;}
#'   \item{3) Specificity}{LC model inferred test specificities (mean across chains, median across iterations) vs simulated test specificities;}
#' }
#' @importFrom ggplot2 ggplot theme theme_classic theme_linedraw geom_abline geom_point geom_errorbar facet_grid labs coord_cartesian element_text element_rect element_blank unit label_both scale_color_manual scale_color_brewer
#' @importFrom dplyr rename left_join filter
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom utils globalVariables
#' @name sim.result.time.plot
#' @examples
#' if (interactive()) {
#' sim_results_time <- run.sims.LC.time(num_tests = 4,
#'                        days=50,
#'                        spec_vec = c(1, 0.98),
#'                        sens_vec = c(1),
#'                        p_performed_vec= c(1),
#'                        sim_size = 10,
#'                        chains = 2,
#'                        iter=500,
#'                        warmup=200,
#'                        data_ID = "sims_time")
#'
#' sim.result.time.plot(sim_results_time)
#' }
#'

utils::globalVariables(c("Simulated_prev", "Simulated_sens", "Simulated_spec", "stan_prev", "stan_prev_CI_low",
                         "stan_prev_CI_high", "sens_median", "Test_ID", "sens_median_CI_low", "sens_median_CI_high",
                         "spec_median", "spec_median_CI_low", "spec_median_CI_high", "overall_test_positivity_sim",
                         "Simulated_overall_test_positivity", "CI_lower", "CI_upper", "model_id", "overall_weekly_test_positivity"))
#' @export
sim.result.time.plot <- function(sim_stan_results) {

  plots <- list()


  sim_stan_results_combined <- dplyr::left_join(sim_stan_results$stan_results_df, sim_stan_results$sim_inputs, by = c("model_id", "test_id"))


  sim_stan_results_combined <-  sim_stan_results_combined %>%
    dplyr::rename(Simulated_sens = sim_sens,
                  Simulated_spec = sim_spec,
                  Test_ID = test_id,
                  Simulated_overall_test_positivity = overall_test_positivity_sim)

  theme <- theme_linedraw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                    strip.background = element_rect(fill = "white"),
                                    strip.text = element_text(size = 11, colour = "black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    legend.key = element_rect(fill = "white", colour = NA),
                                    legend.key.size = unit(1.2, "lines"),
                                    legend.position = "bottom",
                                    legend.text = element_text(size=10), legend.title = element_text(size=12),
                                    axis.title = element_text(size = 14))



  plot_sens_median <- sim_stan_results_combined %>%
    ggplot(aes(x = Simulated_sens, y = sens_median, colour=Test_ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_median_CI_low, ymax = sens_median_CI_high), width = 0, linewidth=1,alpha=0.4) +
    facet_grid(~Simulated_spec, scales = "free",
               labeller = label_both) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "Simuated vs inferred sensitivity", x = "Simulated sensitivity", y = "Inferred sensitivity") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_median <-  plot_sens_median

  plot_spec_median <-  sim_stan_results_combined %>%
    ggplot(aes(x = Simulated_spec, y = spec_median, colour=Test_ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_median_CI_low, ymax = spec_median_CI_high), width = 0, linewidth=1, alpha=0.4) +
    facet_grid(~Simulated_sens, scales = "free",
               labeller = label_both) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "Simuated vs inferred specificity", x = "Simulated specificity", y = "Inferred specificity") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$spec_median <- plot_spec_median

#prev plot over time:
  sim_data <- sim_stan_results$sim_data %>%
    rename(time = day_of_year) #assuming this is col name from sim.test.data.time()
  #assumes time is in days
  sim_data$week <- ((sim_data$time - 1) %/% 7) + 1
  #weekly data
  # Identify test result columns:
  test_result_cols <- sim_data %>%
    dplyr::select(dplyr::starts_with("test_result")) %>%
    names()
  # Combine test result columns into overall test result
  sim_data_week <- sim_data %>%
    dplyr::mutate(test_result_overall = as.factor(
      dplyr::if_else(rowSums(dplyr::across(dplyr::all_of(test_result_cols)) == "1", na.rm = TRUE) > 0, "1", "0"))) %>%
    dplyr::group_by(week, model_id) %>%
    dplyr::summarise(
      overall_weekly_test_positivity = mean(as.numeric(as.character(test_result_overall)), na.rm = TRUE),
      weekly_true_prev = mean(true_prev),
      .groups = "drop"
    )
  stan_data <- sim_stan_results$stan_results_df %>%
    filter(test_id == "test1")
  stan_sim_df_weekly <- dplyr::left_join(stan_data, sim_data_week, by=c('week', 'model_id'))

  #what if sim_prob is different for the same sens/spec?
  weekly_prev_plot <- stan_sim_df_weekly %>%
    dplyr::rename(Spec = "sim_spec", Sens = "sim_sens") %>%
    ggplot(aes(x = week)) +
    geom_point(aes(y = stan_prev, color = "Inferred prevalence")) +
    geom_ribbon(aes(ymin = stan_prev_CI_low, ymax = stan_prev_CI_high), alpha = 0.2) +  #inferred prev CI's
    geom_line(aes(y = weekly_true_prev, color = "True prevalence"), linewidth = 1, alpha=0.8) +
    geom_point(aes(y = overall_weekly_test_positivity, color = "Apparent prevalence"), alpha=0.8) +
    geom_smooth(aes(y = overall_weekly_test_positivity, color = "Apparent prevalence"), se = FALSE, method = "loess", alpha=0.8) +
    facet_grid(Sens ~ Spec , scales = "free",
               labeller = label_both) +
    scale_color_manual(values = c(
                          "Inferred prevalence" = "black",
                          "True prevalence" = "red",
                          "Apparent prevalence" = "blue")) +
    labs(y = "Prevalence", x = "Week") +
    theme

plots$prev_plot <- weekly_prev_plot
  return(plots)
}

