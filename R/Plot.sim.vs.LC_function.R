


#input will be sim test reconstruct? and






#plot individual test results sens/spec, mean and median

#Function to plot results
#medians
sim.result.plot.individual.tests <- function(sim_result_df) {

  plots <- list()

  theme <- theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              strip.background = element_rect(fill = "white"),
                              strip.text = element_text(size = 11),
                              panel.grid.major = element_blank(), # Remove gridlines
                              panel.grid.minor = element_blank(),
                              legend.key = element_rect(fill = "white", colour = NA),  # White background for legend symbols
                              legend.key.size = unit(1.2, "lines"),
                              legend.position = "bottom",
                              legend.text = element_text(size=10), legend.title = element_text(size=12),
                              axis.title = element_text(size = 14))

  plot_spec_median <-  sim_result_df %>%
    ggplot(aes(x = sim_spec, y = spec_median, colour=test_number)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_median_CI_low, ymax = spec_median_CI_high), width = 0, linewidth=1, alpha=0.5) +
    facet_grid(sim_sens ~ sim_prev, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated Vs LC model median specificity", x = "Simulated spec", y = "LC output spec") +
    # xlim(c(0.8, 1)) + ylim(c(0.8, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$spec_median <- plot_spec_median

  plot_sens_median <- sim_result_df %>%
    ggplot(aes(x = sim_sens, y = sens_median, colour=test_number)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_median_CI_low, ymax = sens_median_CI_high), width = 0, linewidth=1,alpha=0.5) +
    facet_grid(sim_spec ~ sim_prev, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated Vs LC model median sensitivity", x = "Simulated sens", y = "LC output sens") +
    # xlim(c(0.2, 1)) + ylim(c(0.2, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_median <-  plot_sens_median

  plot_spec_mean <-  sim_result_df %>%
    ggplot(aes(x = sim_spec, y = spec_mean, colour=test_number)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = spec_mean_CI_low, ymax = spec_mean_CI_high), width = 0, linewidth=1, alpha=0.5) +
    facet_grid(sim_sens ~ sim_prev, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated Vs LC model mean specificity", x = "Simulated spec", y = "LC output spec") +
    # xlim(c(0.8, 1)) + ylim(c(0.8, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$spec_mean <- plot_spec_mean

  plot_sens_mean <- sim_result_df %>%
    ggplot(aes(x = sim_sens, y = sens_mean, colour=test_number)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(size=1, alpha=0.8) +
    geom_errorbar(aes(ymin = sens_mean_CI_low, ymax = sens_mean_CI_high), width = 0, linewidth=1,alpha=0.5) +
    facet_grid(sim_spec ~ sim_prev, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated Vs LC model mean sensitivity", x = "Simulated sens", y = "LC output sens") +
    # xlim(c(0.2, 1)) + ylim(c(0.2, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$sens_mean <- plot_sens_mean

  plot_prev <- sim_result_df %>%
    ggplot(aes(x = sim_prev, y = stan_prev)) +
    geom_abline(intercept = 0, slope = 1, colour = "grey") +
    geom_point(shape=21, fill="pink3", colour="pink3", size=1) +
    geom_errorbar(aes(ymin = stan_prev_CI_low, ymax = stan_prev_CI_high), width = 0, colour="pink3", linewidth=1, alpha=0.6) +
    facet_grid(sim_spec ~ sim_sens, scales = "free",
               labeller = label_both) +
    labs(title = "Simuated Vs LC model prevalence", x = "Simulated prev", y = "LC output prev") +
    # xlim(c(0, 0.2)) + ylim(c(0, 0.2)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme

  plots$prev <- plot_prev

  return(plots)
}


sim.result.plot.individual.tests(stan_sim_results_combined  )
