test_that("LC model runs with sim data (no delay or dependence)", {

  test_params <- list(test1 = list(sens = 0.95, spec = 0.98, p_performed = 1), test2 = list(sens = 0.90, spec = 0.97, p_performed = 0.8),
                      test3 = list(sens = 0.95, spec = 0.98, p_performed = 1), test4 = list(sens = 0.90, spec = 0.97, p_performed = 0.8))


  num_tests=4
  iter=500
  chains=2
  warmup=200
  data_ID = "sim"

  sim_results <- sim.test.data(disease_prev = 0.2, sim_size = 1000, test_params = test_params)

  test <- run.LC.model(data=sim_results$test_results, num_tests=num_tests,
                       iter=iter, chains=chains, warmup=warmup)

  expect_false(is.null(test$stan_fit))

})

test_that("LC model runs with sim data with model dependence (no delay)", {
  #dependence doesn't exist in data but is included in model just to test function
  test_params <- list(test1 = list(sens = 0.95, spec = 0.98, p_performed = 1), test2 = list(sens = 0.90, spec = 0.97, p_performed = 0.8),
                      test3 = list(sens = 0.95, spec = 0.98, p_performed = 1), test4 = list(sens = 0.90, spec = 0.97, p_performed = 0.8))


  num_tests=4
  iter=500
  chains=2
  warmup=200
  data_ID = "sim_dependence"

  sim_results <- sim.test.data(disease_prev = 0.2, sim_size = 1000, test_params = test_params)

  test <- run.LC.model(data=sim_results$test_results, num_tests=num_tests,
                       iter=iter, chains=chains, warmup=warmup,
                       dependency_groups = list(c(1, 2), c(3, 4)))

  expect_false(is.null(test$stan_fit))

})

test_that("LC model runs with sim data with delay (no dependence)", {

  test_params <- list(test1 = list(sens = 0.95, spec = 0.98, p_performed = 1), test2 = list(sens = 0.90, spec = 0.97, p_performed = 0.8),
                      test3 = list(sens = 0.95, spec = 0.98, p_performed = 1), test4 = list(sens = 0.90, spec = 0.97, p_performed = 0.8))


  num_tests=4
  iter=500
  chains=2
  warmup=200
  data_ID = "sim_delay"

  sim_results <- sim.test.data(disease_prev = 0.2, sim_size = 1000, test_params = test_params, delay = TRUE)

  test <- run.LC.model(data=sim_results$test_results, num_tests=num_tests,
                       iter=iter, chains=chains, warmup=warmup,
                       covariates = c("delay"))

  expect_false(is.null(test$stan_fit))

})

test_that("LC model simulation runs", {

  num_tests=4
  iter=500
  chains=2
  warmup=200
  prev_vec = c(0.1, 0.2)
  spec_vec = c(1, 0.98)
  sens_vec = c(1)
  p_performed_vec= c(1)

  results <- run.sims.LC(num_tests = num_tests,
                         prev_vec = prev_vec,
                         spec_vec = spec_vec,
                         sens_vec = sens_vec,
                         p_performed_vec= p_performed_vec,
                         sim_size = 100, chains = chains, iter = iter, warmup=warmup, data_ID = "sims")

  expect_false(is.null(results$sim_inputs))
  expect_false(is.null(results$stan_results_df))
  expect_false(is.null(results$divergence_summary))

  # Calculate expected number of parameter combinations
  n_param_combos <- length(prev_vec) * length(spec_vec) * length(sens_vec) * length(p_performed_vec)
  # Expect number of rows in sim_inputs to equal num_tests * param_combos
  expect_equal(nrow(results$sim_inputs), num_tests * n_param_combos)
  # Expect one row per param combo in stan_results_df
  expect_equal(nrow(results$stan_results_df), n_param_combos)
})
