#' @title Identifies divergent transitions and low bulk ESS in a fitted stan model.
#' @description Counts the number of divergent transitions in a fitted stan model and specifies whether the bulk Effective Sample Size (ESS) is too low
#' for prev/sens/spec parameters.The function focuses only on these top-level population parameters. Individual-level parameters (such as random effects) are ignored
#' when checking for low ESS. Low ESS is based on a threshold of 100 x number of chains.
#'
#' @param fit Stan model fit object.
#' @param model_name Optional model name identifier. Default = NA.
#' @return Dataframe with the specified model name, the number of divergent transitions, and logical indicating if bulk ESS is too low for prev/sens/spec parameters.
#' @importFrom rstan get_sampler_params summary
#' @importFrom dplyr filter group_by summarise mutate n cross_join select case_when rowwise left_join
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom stats rbinom na.omit quantile
#' @importFrom stringr str_detect
#' @importFrom utils globalVariables
#' @name check.divergent.transitions
#' @examples
#' if (interactive()) {
#' # Simulate data
#' sim_data <- sim.test.data(sim_size = 100)
#'
#' # Run LC model using simulated test results
#' fit <- run.LC.model(sim_data$test_results, num_tests = 4,
#'                       data_ID = "sim", model_name = "basic_sim",
#'                       dependency_groups = NULL, covariates = NULL,
#'                       iter=1000, chains=4, warmup=500)
#'
#' check.divergent.transitions(fit$stan_fit, model_name = "Test_model")
#' }
#'

# Function to check model for divergent transitions
#' @export
check.divergent.transitions <- function(fit, model_name = NA_character_) {

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- do.call(rbind, sampler_params)[, "divergent__"]
  num_divergent <- sum(divergent == 1)

  summary_df <- rstan::summary(fit)$summary
  param_names <- rownames(summary_df)
  num_chains <- length(sampler_params)
  ess_threshold <- 100 * num_chains

  target_params <- c("prev", "Se_median", "Sp_median")
  # Find matching parameters
  matches <- param_names[grepl(paste0("^(", paste(target_params, collapse = "|"), ")"), param_names)]

  if (length(matches) == 0) {
    warning("No target parameters found in model output.")
    return(data.frame(model_id = model_name, divergent_count = num_divergent, low_ESS = NA))
  }

  # Check bulk ESS for matching parameters
  ess_bulk <- summary_df[matches, "n_eff"]
  any_low_ess <- any(ess_bulk < ess_threshold, na.rm = TRUE)

  msg <- paste0("Model: ", if (is.na(model_name)) "unnamed" else model_name,
                " -Divergent transitions: ", num_divergent,
                ", Prev/sens/spec parameters with low ESS: ", any_low_ess)
  message(msg)

  if (num_divergent > 0) {
    warning(paste("Model", if (!is.na(model_name)) paste0("(", model_name, ")"),
                  "has", num_divergent, "divergent transitions."))
  }

  if (any_low_ess == TRUE) {
    warning(paste("Model", model_name, "has prev/sens/spec parameters with low effective sample size (ESS < ", ess_threshold, ")."))
  }


  return(data.frame(model_id = model_name, divergent_count = num_divergent, low_ess = any_low_ess))
}



