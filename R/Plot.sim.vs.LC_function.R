
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
#' @importFrom dplyr rename left_join filter pull
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom utils globalVariables
#' @importFrom RColorBrewer brewer.pal
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
                         "Simulated_overall_test_positivity", "CI_lower", "CI_upper",
                         "divergent_count", "low_ess", "Test_ID_divergent", "Prev", "Inferred_prev_divergent", "Sens", "Spec"))
#' @export
sim.result.plot <- function(sim_stan_results) {

  plots <- list()

  sim_stan_results_combined1 <- dplyr::left_join(sim_stan_results$stan_results_df, sim_stan_results$sim_inputs, by = c("model_id", "test_id"))
  sim_stan_results_combined2 <- left_join(sim_stan_results_combined1 , sim_stan_results$divergence_summary, by = c("model_id"))

  sim_stan_results_combined <-  sim_stan_results_combined2 %>%
    dplyr::rename(Prev = sim_prev,
                  Sens = sim_sens,
                  Spec = sim_spec,
                  Test_ID = test_id,
                  Simulated_overall_test_positivity = overall_test_positivity_sim) %>%
    mutate(Inferred_prev_divergent = factor(ifelse((divergent_count > 0 | low_ess == TRUE), "Not converged", "Inferred prevalence"))) %>%
    mutate(Test_ID_divergent = factor(ifelse((divergent_count > 0 | low_ess == TRUE), "Not converged", Test_ID)))

  #colours
  dark2_cols <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  # list of test_id levels
  test_ids <- sim_stan_results_combined %>%
    filter(Test_ID_divergent != "Not converged") %>%
    dplyr::pull(Test_ID_divergent) %>%
    unique()
  # named color vector
  named_colors <- setNames(dark2_cols[1:length(test_ids)], test_ids)
  named_colors["Not converged"] <- "grey85"

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
    ggplot(aes(x = Prev)) +
    geom_abline(intercept = 0, slope = 1, colour = "#D55E00", alpha=0.5) +
    geom_point(aes(y = stan_prev, colour=Inferred_prev_divergent, fill=Inferred_prev_divergent), shape=21, size=1) +
    geom_errorbar(aes(ymin = stan_prev_CI_low, ymax = stan_prev_CI_high, colour=Inferred_prev_divergent), width = 0, linewidth=0.5, alpha=0.5) +
    #overall test positivity (assuming any test positive is a positive):
    geom_point(aes(y = Simulated_overall_test_positivity, color = "Test positivity", fill="Test positivity"), shape=21, size=0.8, alpha=0.8) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, color = "Test positivity"), width = 0, linewidth=0.5, alpha=0.5) +
    facet_grid(Spec ~ Sens, scales = "free",
               labeller = label_both) +
    scale_color_manual(values = c("Inferred prevalence" = "black", "Not converged" = "grey85", "Test positivity" = "#009E73")) +
    scale_fill_manual(values = c("Inferred prevalence" = "black", "Not converged" = "grey85", "Test positivity" = "#009E73")) +
    labs(title = "Simuated vs inferred prevalence", x = "Simulated prevalence", y = "Inferred prevalence", color = NULL, fill = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme +
    guides(fill = "none", color = guide_legend() )

  plots$prev <- plot_prev

  plot_sens_median <- sim_stan_results_combined %>%
    ggplot(aes(x = Sens, y = sens_median, colour=Test_ID_divergent)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_median_CI_low, ymax = sens_median_CI_high), width = 0, linewidth=0.5,alpha=0.4) +
    facet_grid(Spec ~ Prev, scales = "free",
               labeller = label_both) +
    scale_color_manual(values = named_colors) +
    labs(title = "Simulated vs inferred sensitivity", x = "Simulated sensitivity", y = "Inferred sensitivity", color = NULL, fill = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_median <-  plot_sens_median

  plot_spec_median <-  sim_stan_results_combined %>%
    ggplot(aes(x = Spec, y = spec_median, colour=Test_ID_divergent)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_median_CI_low, ymax = spec_median_CI_high), width = 0, linewidth=0.5, alpha=0.4) +
    facet_grid(Sens ~ Prev, scales = "free",
               labeller = label_both) +
    scale_color_manual(values = named_colors) +
    labs(title = "Simulated vs inferred specificity", x = "Simulated specificity", y = "Inferred specificity", color = NULL, fill = NULL) +
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
#' @importFrom dplyr rename left_join filter pull
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom utils globalVariables
#' @importFrom RColorBrewer brewer.pal
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
                         "Simulated_overall_test_positivity", "CI_lower", "CI_upper", "model_id", "overall_weekly_test_positivity",
                         "divergent_count", "low_ess", "Test_ID_divergent", "Prev", "Inferred_prev_divergent", "Sens", "Spec"))

#' @export sim.result.time.plot
sim.result.time.plot <- function(sim_stan_results) {

  plots <- list()

  sim_stan_results_combined1 <- dplyr::left_join(sim_stan_results$stan_results_df, sim_stan_results$sim_inputs, by = c("model_id", "test_id"))
  sim_stan_results_combined2 <- left_join(sim_stan_results_combined1 , sim_stan_results$divergence_summary, by = c("model_id"))

  sim_stan_results_combined <-  sim_stan_results_combined2 %>%
    dplyr::rename(Sens = sim_sens,
                  Spec = sim_spec,
                  Test_ID = test_id,
                  Simulated_overall_test_positivity = overall_test_positivity_sim) %>%
    mutate(Inferred_prev_divergent = factor(ifelse((divergent_count > 0 | low_ess == TRUE), "Not converged", "Inferred prevalence"))) %>%
    mutate(Test_ID_divergent = factor(ifelse((divergent_count > 0 | low_ess == TRUE), "Not converged", Test_ID)))

  #colours
  dark2_cols <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
  # list of test_id levels
  test_ids <- sim_stan_results_combined %>%
    filter(Test_ID_divergent != "Not converged") %>%
    dplyr::pull(Test_ID_divergent) %>%
    unique()
  # named color vector
  named_colors <- setNames(dark2_cols[1:length(test_ids)], test_ids)
  named_colors["Not converged"] <- "grey85"

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
    ggplot(aes(x = Sens, y = sens_median, colour=Test_ID_divergent)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_median_CI_low, ymax = sens_median_CI_high), width = 0, linewidth=0.5,alpha=0.4) +
    facet_grid(~Spec, scales = "free",
               labeller = label_both) +
    scale_color_manual(values = named_colors) +
    labs(title = "Simulated vs inferred sensitivity", x = "Simulated sensitivity", y = "Inferred sensitivity", color = NULL, fill = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_median <-  plot_sens_median

  plot_spec_median <-  sim_stan_results_combined %>%
    ggplot(aes(x = Spec, y = spec_median, colour=Test_ID_divergent)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_median_CI_low, ymax = spec_median_CI_high), width = 0, linewidth=0.5, alpha=0.4) +
    facet_grid(~Sens, scales = "free",
               labeller = label_both) +
    scale_color_manual(values = named_colors) +
    labs(title = "Simulated vs inferred specificity", x = "Simulated specificity", y = "Inferred specificity", color = NULL, fill = NULL) +
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

  stan_sim_df_weekly1 <- dplyr::left_join(stan_data, sim_data_week, by=c('week', 'model_id'))
  stan_sim_df_weekly2 <- dplyr::left_join(stan_sim_df_weekly1, sim_stan_results$divergence_summary, by = c("model_id"))

  stan_sim_df_weekly <- stan_sim_df_weekly2 %>%
    mutate(Inferred_prev_divergent = factor(ifelse((divergent_count > 0 | low_ess == TRUE), "Not converged", "Inferred prevalence")))

  #what if sim_prob is different for the same sens/spec?
  weekly_prev_plot <- stan_sim_df_weekly %>%
    dplyr::rename(Spec = "sim_spec", Sens = "sim_sens") %>%
    ggplot(aes(x = week)) +
    geom_point(aes(y = stan_prev, colour=Inferred_prev_divergent), size=1) +
    geom_ribbon(aes(ymin = stan_prev_CI_low, ymax = stan_prev_CI_high, fill = Inferred_prev_divergent), alpha = 0.4, colour=NA) +  #inferred prev CI's
    geom_line(aes(y = weekly_true_prev, color = "True prevalence"),  alpha=0.6, linewidth = 1) +
    geom_point(aes(y = overall_weekly_test_positivity, color = "Test positivity"), size=1, alpha=0.6) +
    geom_line(aes(y = overall_weekly_test_positivity, color = "Test positivity"), stat="smooth", method="loess", alpha=0.6, linewidth = 1) + #using geom_line instead of geom_smooth so can use alpha on line colour
    facet_grid(Spec ~ Sens, scales = "free",
               labeller = label_both) +
    scale_color_manual(values = c("Inferred prevalence" = "black", "Not converged" = "grey", "True prevalence" = "#D55E00", "Test positivity" = "#009E73")) +
    scale_fill_manual(values = c("Inferred prevalence" = "black", "Not converged" = "grey", "Test positivity" = "#009E73")) +
    labs(y = "Prevalence", x = "Week", color = NULL, fill = NULL) +
    theme +
    guides(fill = "none", color = guide_legend() )

plots$prev_plot <- weekly_prev_plot

  return(plots)
}

