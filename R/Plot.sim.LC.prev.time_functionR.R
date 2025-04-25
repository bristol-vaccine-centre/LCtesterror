
#' @title Plots inferred weekly prevalence from run.LC.model with covariates = c("Time")
#' @description Plots inferred weekly prevalence from run.LC.model with covariates = c("Time"). If sim.test.data.time() was used to simulate true prevalence and test results,
#' this function output can be provided to the sim_time_data argument to compare inferred prev with true prevalence.
#' @param run.LC.model_output Returned run.LC.model() model object.
#' @param sim_time_data Returned output from sim.test.data.time() if simulated data were used. Allows plotting of true prevalence over time alongside inferred. Default = NULL.
#' @param time_resolution Time resolution of test data used in run.LC.model(). Either "days" or "weeks". Default = "days".
#' @return Stan model fit and various summary outputs:
#'   \describe{
#'   \item{weekly_inferred_prev_plot}{Plot of weekly inferred prevalence from the LC model which adjusts test data for test error.}
#'   \item{weekly_prev_plot}{Plot of weekly inferred prevalence from the LC model which adjusts test data for test error, with plotted test positivity (apparent prevalence), and true prevalence if simulated}
#'   \item{weekly_data}{Data aggregated to weekly}
#' }
#' @export
#' @importFrom dplyr group_by mutate rename case_when left_join if_else across all_of summarise
#' @importFrom ggplot2 ggplot aes geom_point geom_ribbon geom_line geom_smooth labs theme_bw scale_color_manual scale_fill_manual guides guide_legend
#' @name prev.time.plots
#'
#' @examples
#' if (interactive()) {
#' # Example with simulated data
#' # (In this example the model will not converge
#' # and CI's will be very wide because of small amount of simulated data)
#' sim_data <- sim.test.data.time(sim_size=10)
#' fit <- run.LC.model(sim_data$test_results, num_tests = 4, covariates = c("time"))
#' plots <- prev.time.plots(fit, sim_time_data = sim_data)
#' plots$weekly_prev_plot
#'  }

utils::globalVariables(c(
  "week", "mean", "CI_min", "CI_max", "mean_prevalence",
  "weekly_test_positivity", "weekly_true_prev",
  "test_result_overall", "true_prev"
))


prev.time.plots <- function(run.LC.model_output, sim_time_data = NULL, time_resolution = "days") {

  stan_fit_summary_df <- run.LC.model_output$stan_fit_summary_df

#Prev over time
prev_time_rows <- grepl("^weekly_prevalence", rownames(stan_fit_summary_df))
stan_fit_summary_prev_time <- subset(stan_fit_summary_df, prev_time_rows)
stan_fit_summary_prev_time_df <- stan_fit_summary_prev_time %>%
  dplyr::select(mean, '2.5%', '97.5%') %>%
  dplyr::mutate(week = row_number()) %>%
  dplyr::rename(CI_min = '2.5%', CI_max = '97.5%', mean_prevalence = mean)

weekly_inferred_prev_plot <- stan_fit_summary_prev_time_df %>%
  ggplot(aes(y=mean_prevalence, x=week)) +
  geom_point() +
  geom_ribbon(aes(ymin = CI_min, ymax = CI_max), alpha=0.3) +
  theme_bw()

#if no data is provided
if(is.null(sim_time_data)) {

sim_data <- run.LC.model_output$test_data_input
names(sim_data)[ncol(sim_data)] <- "time"

if (time_resolution == "days") {
sim_data$week <- ((sim_data$time - 1) %/% 7) + 1
 } else if (time_resolution == "weeks") {
   sim_data$week <- sim_data$time

 }

#weekly data
# Identify test result columns: everything except 'time', 'week'
non_test_cols <- c("time", "week")
test_result_cols <- setdiff(names(sim_data), non_test_cols)
# Combine test result columns into overall test result
sim_data_week <- sim_data %>%
  dplyr::mutate(test_result_overall = as.factor(
    dplyr::if_else(rowSums(dplyr::across(dplyr::all_of(test_result_cols)) == "1", na.rm = TRUE) > 0, "1", "0"))) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(
    weekly_test_positivity = mean(as.numeric(as.character(test_result_overall)), na.rm=TRUE),
    .groups = "drop"
  )

stan_sim_df_weekly <- dplyr::left_join(stan_fit_summary_prev_time_df, sim_data_week, by=c('week'))

weekly_prev_plot <- stan_sim_df_weekly %>%
  ggplot(aes(x = week)) +
  geom_point(aes(y = mean_prevalence, color = "Inferred prevalence")) +
  geom_ribbon(aes(ymin = CI_min, ymax = CI_max), alpha = 0.2) +  #inferred prev CI's
  geom_point(aes(y = weekly_test_positivity, color = "Apparent prevalence"), alpha=0.8) +
  geom_smooth(aes(y = weekly_test_positivity, color = "Apparent prevalence"), se = FALSE, method = "loess", alpha=0.8) +
  scale_color_manual( name = "Legend",
                      values = c(
                        "Inferred prevalence" = "black",
                        "Apparent prevalence" = "blue")) +
  labs(y = "Prevalence", x = "Week") +
  theme_bw()

}

#if provided data is from sim.test.data.time()
if(!is.null(sim_time_data)) {

  sim_data <- sim_time_data$sim_data

  sim_data <- sim_data %>%
    rename(time = day_of_year) #assuming this is col name from sim.test.data.time()

  if (time_resolution == "days") {
    sim_data$week <- ((sim_data$time - 1) %/% 7) + 1
  } else if (time_resolution == "weeks") {
    sim_data$week <- sim_data$time #assuming column naming is the same

  }

  #weekly data
  # Identify test result columns:
  test_result_cols <- sim_data %>%
             dplyr::select(dplyr::starts_with("test_result")) %>%
              names()
  # Combine test result columns into overall test result
  sim_data_week <- sim_data %>%
    dplyr::mutate(test_result_overall = as.factor(
      dplyr::if_else(rowSums(dplyr::across(dplyr::all_of(test_result_cols)) == "1", na.rm = TRUE) > 0, "1", "0"))) %>%
    dplyr::group_by(week) %>%
    dplyr::summarise(
      weekly_test_positivity = mean(as.numeric(as.character(test_result_overall)), na.rm = TRUE),
      weekly_true_prev = mean(true_prev),
      .groups = "drop"
    )

  stan_sim_df_weekly <- dplyr::left_join(stan_fit_summary_prev_time_df, sim_data_week, by=c('week'))

  weekly_prev_plot <- stan_sim_df_weekly %>%
    ggplot(aes(x = week)) +
    geom_point(aes(y = mean_prevalence, color = "Inferred prevalence")) +
    geom_ribbon(aes(ymin = CI_min, ymax = CI_max), alpha = 0.2) +  #inferred prev CI's
    geom_line(aes(y = weekly_true_prev, color = "True prevalence"), linewidth = 1, alpha=0.8) +
    geom_point(aes(y = weekly_test_positivity, color = "Apparent prevalence"), alpha=0.8) +
    geom_smooth(aes(y = weekly_test_positivity, color = "Apparent prevalence"), se = FALSE, method = "loess", alpha=0.8) +
    scale_color_manual( name = "Legend",
      values = c(
        "Inferred prevalence" = "black",
        "True prevalence" = "red",
        "Apparent prevalence" = "blue")) +
    labs(y = "Prevalence", x = "Week") +
    theme_bw()

}



results <- list()

results$weekly_inferred_prev_plot <- weekly_inferred_prev_plot
results$weekly_prev_plot <- weekly_prev_plot
results$weekly_data <- stan_sim_df_weekly

return(results)

 }
