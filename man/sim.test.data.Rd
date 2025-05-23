% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sim.tests.dynamic_function.R
\name{sim.test.data}
\alias{sim.test.data}
\title{Simulates test data based on true prevalence, sensitivity and specificity.}
\arguments{
\item{disease_prev}{True population disease prevalence to simulate (between 0-1). Default = 0.2.}

\item{sim_size}{Number of individuals to simulate test results for. Default = 1000.}

\item{test_params}{Test parameters used to simulate test results along with true prevalence. Given as a list of lists for each test containing sensitivity (sens =), specificity (spec =), and probability that the test is performed (p_performed = ).
Default = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1), test2 = list(sens = 0.99, spec = 0.99, p_performed = 1), test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8), test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),}

\item{seed}{for set.seed(). Default = 953.}

\item{delay}{Logical indicating whether to simulate effects of 'delay until testing' on sensitivity. Default = FALSE.}

\item{delay_distribution}{Optional specified distribution of delay effect. A function that takes an integer \code{n} and returns a vector of length \code{n} representing individual-specific delays. Default = sample(0:14, n, replace = TRUE).}

\item{delay_effect_fn}{Optional function to simulate effects of 'delay until testing' on sensitivity. Default = function(delay_day, sens) pmax(sens - 0.02 * delay_day, 0.5)}
}
\value{
A table of test parameters (specified and simulated for comparison) and a table containing binary test results for each individual
\describe{
\item{test_parameters}{Test results table with row for each test (test_id) containing specified test parameters (sens; spec; p_performed; disease_prev), simulated test_positivity, test_coverage (based on p_performed), and the estimated true disease prevalence estimate (disease_prev_est: based on the Rogan-Gladen equation) }
\item{test_results}{Simulated binary test results for each individual (N=sim_size), based on specified true prev, sens and spec parameters.
If delay = TRUE, additional delay columns are included specifying the delay for each test.}
}
}
\description{
Simulates test data for individuals based on true prevalence, and test sensitivity and specificity.
Can simulate multiple test results for each individual with different test sens/spec parameters and different probabilities that each test is performed.
The simulated data does not take into account dependencies between tests.
}
\examples{
if (interactive()) {
sim_data <- sim.test.data(disease_prev = 0.2, sim_size = 100,
test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1),
test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8),
test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)))

#head(sim_data$test_results)

# Example with 'delay until test' simulated
sim_data <- sim.test.data(disease_prev = 0.2, sim_size = 100,
test_params = list(test1 = list(sens = 0.99, spec = 0.99, p_performed = 1),
test2 = list(sens = 0.99, spec = 0.99, p_performed = 1),
test3 = list(sens = 0.98, spec = 0.98, p_performed = 0.8),
test4 = list(sens = 0.98, spec = 0.98, p_performed = 0.8)),
delay = TRUE)
}

}
