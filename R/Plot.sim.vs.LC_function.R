
#' @title Plots LC model simulation results.
#' @description Plots LC model simulation results generated from run.sims.LC() to evaluate model performance at different parameter values.
#' Compares LC model inferred prevalence, test sensitivities and test specificities against simulated parameters.
#' Will plot across any simulated any value of the probability of the test occurring (sim_prob),
#' so if multiple p_performed were simulated in run.sims.LC(), these will need to be pre-filtered to plot each sim_prob paramter value one at a time.
#'
#' @param sim_stan_results Dataframe output from run.sims.LC()$stan_results_df.
#'
#' @return Plots of inferred parameters against simulated parameters. Returns 3 plots:
#' #' \describe{
#'   \item{1) Prevalence}{LC model inferred prevalence (mean across chains, mean across iterations)) vs simulated prevalence;}
#'   \item{2) Sensitivity}{LC model inferred test sensitivities (mean across chains, median across iterations) vs simulated test sensitivities;}
#'   \item{3) Specificity}{LC model inferred test specificities (mean across chains, median across iterations) vs simulated test specificities;}
#' }
#' @importFrom ggplot2 ggplot theme theme_classic geom_abline geom_point geom_errorbar facet_grid labs coord_cartesian element_text element_rect element_blank unit label_both
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom utils globalVariables
#' @name sim.result.plot
#' @export
#' @examples
#' if (interactive()) {
#' sim_results <- run.sims.LC(num_tests = 4,
#'                        prev_vec = c(0.1, 0.2),
#'                        spec_vec = c(1),
#'                        sens_vec = c(1),
#'                        p_performed_vec= c(1),
#'                        sim_size = 100,
#'                        chains = 2,
#'                        data_ID = "sims")
#'
#' sim.result.plots(sim_results$stan_results_df)
#' }
#'

utils::globalVariables(c("Simulated_prev", "Simulated_sens", "Simulated_spec", "stan_prev", "stan_prev_CI_low",
                         "stan_prev_CI_high", "sens_median", "Test_ID", "sens_median_CI_low", "sens_median_CI_high",
                         "spec_median", "spec_median_CI_low", "spec_median_CI_high"))

sim.result.plot <- function(sim_stan_results) {

  plots <- list()

  sim_stan_results <- sim_stan_results %>%
    dplyr::rename(Simulated_prev = sim_prev,
                  Simulated_sens = sim_sens,
                  Simulated_spec = sim_spec,
                  Test_ID = test_id)

  theme <- theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              strip.background = element_rect(fill = "white"),
                              strip.text = element_text(size = 11),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.key = element_rect(fill = "white", colour = NA),
                              legend.key.size = unit(1.2, "lines"),
                              legend.position = "bottom",
                              legend.text = element_text(size=10), legend.title = element_text(size=12),
                              axis.title = element_text(size = 14))


  plot_prev <- sim_stan_results %>%
    ggplot(aes(x = Simulated_prev, y = stan_prev)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(shape=21, fill="black", colour="black", size=1) +
    geom_errorbar(aes(ymin = stan_prev_CI_low, ymax = stan_prev_CI_high), width = 0, colour="black", linewidth=1, alpha=0.6) +
    facet_grid(Simulated_spec ~ Simulated_sens, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated vs inferred prevalence", x = "Simulated prevalence", y = "Inferred prevalence") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$prev <- plot_prev

  plot_sens_median <- sim_stan_results %>%
    ggplot(aes(x = Simulated_sens, y = sens_median, colour=Test_ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_median_CI_low, ymax = sens_median_CI_high), width = 0, linewidth=1,alpha=0.5) +
    facet_grid(Simulated_spec ~ Simulated_prev, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated vs inferred sensitivity", x = "Simulated sensitivity", y = "Inferred sensitivity") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_median <-  plot_sens_median

  plot_spec_median <-  sim_stan_results %>%
    ggplot(aes(x = Simulated_spec, y = spec_median, colour=Test_ID)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_median_CI_low, ymax = spec_median_CI_high), width = 0, linewidth=1, alpha=0.5) +
    facet_grid(Simulated_sens ~ Simulated_prev, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated vs inferred specificity", x = "Simulated specificity", y = "Inferred specificity") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$spec_median <- plot_spec_median


  return(plots)
}


