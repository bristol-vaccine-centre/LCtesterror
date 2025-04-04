
#Function for running lots of LC runs using sim data with different parameters


run.sims.LC <- function(prev_vec= c(0.5), spec_vec= c(1), sens_vec= c(1), p_performed_vec= c(1), sim_size=1000) {

  for (s in 1:length(sens_vec)) {
    for (c in 1:length(spec_vec)) {
      for (p in 1:length(prev_vec)) {

        list(test1 = list(sens = 0.95, spec = 0.98, p_performed = 1), test2 = list(sens = 0.90, spec = 0.97, p_performed = 0.8))

        sim_data <- sim.test.data(prev_vec[p], sim_size=sim_size)

        result_name <- paste("sens_", sens_vec[s], "_spec_", spec_vec[c], "_prev_", prev_vec[p], sep = "")

        result <- sim.stan.model(rsv_prev = prev_vec[p], sim_size = 10000,
                                 test1_sens = sens_vec[s] , test2_sens = sens_vec[s], test3_sens = sens_vec[s], test4_sens = sens_vec[s],
                                 test1_spec = spec_vec[c], test2_spec = spec_vec[c], test3_spec = spec_vec[c], test4_spec = spec_vec[c],
                                 test1_prob = 1, test2_prob = 1, test3_prob = 1, test4_prob = 1,
                                 stan_model_code = stan_model_code, iter = iter, chain = chain, warmup = warmup)

      }
    }
  }

}


list(test1 = list(sens = 0.95, spec = 0.98, p_performed = 1), test2 = list(sens = 0.90, spec = 0.97, p_performed = 0.8))
