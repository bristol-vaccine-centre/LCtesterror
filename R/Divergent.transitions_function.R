#' @title Counts the number of divergent transitions in a fitted stan model.
#' @description Counts the number of divergent transitions in a fitted stan model.
#'
#' @param fit Stan model fit object.
#' @param model_name Optional model name identifier.
#' @return Specified model name with the number of divergent transitions
#' @export
#' @importFrom dplyr filter group_by summarise mutate n cross_join select case_when rowwise left_join
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom stats rbinom na.omit quantile
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
check.divergent.transitions <- function(fit, model_name = NULL) {

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent <- do.call(rbind, sampler_params)[, "divergent__"]
  num_divergent <- sum(divergent == 1)

  msg <- paste("Divergent transitions:", num_divergent)
  if (!is.null(model_name)) {
    msg <- paste("Model:", model_name, "-", msg)
  }
  message(msg)

  if (num_divergent > 0) {
    warning(paste("Model", if (!is.null(model_name)) paste0("(", model_name, ")"),
                  "has", num_divergent, "divergent transitions."))
  }

  return(data.frame(model = model_name, divergent_count = num_divergent))
}
