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
