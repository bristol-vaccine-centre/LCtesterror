test_that("LC model runs with generate sim data (no delay or dependence)", {

  test_params <- list(test1 = list(sens = 0.95, spec = 0.98, p_performed = 1), test2 = list(sens = 0.90, spec = 0.97, p_performed = 0.8),
                      test3 = list(sens = 0.95, spec = 0.98, p_performed = 1), test4 = list(sens = 0.90, spec = 0.97, p_performed = 0.8))

  test_params

  num_tests=4
  iter=500
  chains=2
  warmup=200
  data_ID = "sim"

  sim_results <- sim.test.data(disease_prev = 0.2, sim_size = 1000, test_params = test_params)
  sim_results

  test <- run.LC.model(data=sim_results$test_results, num_tests=num_tests,
                       iter=iter, chains=chains, warmup=warmup)
})
