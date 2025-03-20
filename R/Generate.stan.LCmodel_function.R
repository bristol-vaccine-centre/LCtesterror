
#' Creates stan model code for different versions of the LC model that can be run in R using rstan.
#'
#' @param num_tests Number of tests in data.
#' @param include_time Logical indicating whether time should be included in the model. Default = FALSE.
#' @param include_delay Logical indicating whether delay until testing should be included in the model. Default = FALSE.
#' @param dependency_groups A list of vectors specifying dependencies between tests. All test numbers that are not independent of each other are specified in a single vector. Default = NULL.
#' @return Stan LC model code.
#' @export
#' 
generate.stan.model <- function(num_tests, include_time = FALSE, include_delay = FALSE, dependency_groups = list()) {
  
  # Identify valid dependency groups (ignoring single-test groups)
  valid_groups <- Filter(function(g) length(g) > 1, dependency_groups)  
  num_valid_groups <- length(valid_groups)
  G <- num_valid_groups
  print(paste("valid_groups:", paste(valid_groups, collapse=", ")))
  print(paste("num_valid_groups:", paste(num_valid_groups, collapse=", ")))
  print(paste("G:", paste(G, collapse=", ")))
  #b_pos/b_neg index vector
  b_index <- rep(0, num_tests)  # Default 0 (no bpos)
  group_counter <- 1
  for (g in valid_groups) {  # Iterate over valid dependency groups
    for (t in g) {
      b_index[t] <- group_counter  # Assign bpos group number
    }
    group_counter <- group_counter + 1
  }
  print(paste("b_index:", paste(b_index, collapse=", ")))
  # Determine if dependency terms (bpos/bneg) should be added
  all_tests_separate <- length(dependency_groups) == num_tests && all(sapply(dependency_groups, length) == 1) #are all tests in seperate groups? If so add dependency will = FALSE
  add_dependency <- !(num_valid_groups == 0 || 
                        (length(dependency_groups) == 1 && length(dependency_groups[[1]]) == num_tests) ||
                        all_tests_separate)
  print(paste("add_dependency:", paste(add_dependency, collapse=", ")))
  # Create a vector: 1 if test has a valid dependency group, 0 otherwise
  dependency_group_existence <- rep(0, num_tests)
  for (group in valid_groups) {
    dependency_group_existence[group] <- 1
  }
  print(paste(" dependency_group_existence:", paste( dependency_group_existence, collapse=", ")))
  
  T <- as.numeric(num_tests)
  
  # Start building Stan code
  stan_code <- "
  data {
    int<lower=1> T; // Number of different tests
    int<lower=1> N; // Number of individuals
    int<lower=1> O; // Number of observations
    array[O] int<lower=1, upper=T> tt; // Test index for each observation
    array[O] int<lower=0, upper=1> y; // Test results (0/1)
    int s[N]; // Number of tests per individual
    "
  if (add_dependency) {
    stan_code <- paste0(stan_code,
                        "array[T] int test_to_group;  //Dependency group mapping for each test - think i can remove
    int<lower=1> G; //Number of dependency groups
                        ")
  }
  if (include_time) {
    stan_code <- paste0(stan_code, "vector[N] time_points;\n")
  }
  
  if (include_delay) {
    stan_code <- paste0(stan_code, "matrix[N, T] delay_days;\n")
  }
  
  stan_code <- paste0(stan_code, 
                      "} \n transformed data { \n")
  
  if (include_delay) {
    stan_code <- paste0(stan_code, "matrix[N, T] delay_unscaled = delay_days;\n")
  }
  
  stan_code <- paste0(stan_code, 
                      "} \n parameters { \n vector[N] RE; \n")
  
  if (add_dependency) {
    stan_code <- paste0(stan_code, 
                        "vector<lower=0,upper=4>[G] bpos;
     vector<lower=0,upper=4>[G] bneg;
     ")
  }
  
  for (t in 1:T){
    stan_code <- paste0(stan_code,
                        "real<lower=0,upper=1> a",t,"1;
                     ") 
  }
  
  for (t in 1:T){
    if(dependency_group_existence[t] == 1) {
      stan_code <- paste0(stan_code,
                          "real<lower=1-inv_logit(logit(a",t,"1)-bpos[",b_index[t],"]*2),upper=1> a",t,"2;
                    ")
    } else {
      stan_code <- paste0(stan_code, 
                          "real<lower=1-inv_logit(logit(a",t,"1)*2),upper=1> a",t,"2;
                    ") 
    }
  } 
  
  if (include_time) {
    stan_code <- paste0(stan_code, "
    real<lower=0,upper=1> prev_amplitude;
    real<lower=0,upper=1> prev_baseline;
    real<lower=0> sigma_gaussian;
    real<lower=0, upper=52> mean_gaussian;
    ")
  } else {
    stan_code <- paste0(stan_code, "
    real<lower=0,upper=1> prev;
      ")
  }
  
  if (include_delay) {
    stan_code <- paste0(stan_code, "real<upper=0> delay_pos;\n")
  }
  
  stan_code <- paste0(stan_code, "
  }

  transformed parameters {
    simplex[2] theta;
    vector[N] prob[T, 2];
    
 ")
  
  if (include_time) {
    stan_code <- paste0(stan_code, "
  for (n in 1:N) {
    real seasonal_effect = prev_amplitude * exp(-0.5 * square((time_points[n] - mean_gaussian) / sigma_gaussian));
    real prevalence = prev_baseline + seasonal_effect;
    
    prevalence = fmin(fmax(prevalence, 0.0), 1.0);
    
    theta[1] = 1 - prevalence;
    theta[2] = prevalence;
   }
  ") 
  } else {
    stan_code <- paste0(stan_code, "
    theta[1] = 1 - prev;
    theta[2] = prev;
    
     ")
  }
  
  if (include_delay) {
    stan_code <- paste0(stan_code, "matrix[N, T] delay = delay_unscaled / 10;\n")
  }
  
  for (t in 1:T){
    if(dependency_group_existence[t] == 1) {
      stan_code <- paste0(stan_code,
                          "prob[",t,", 1] = inv_logit(logit(1 - a",t,"1) + bneg[",b_index[t],"] * RE);
                              ")
      if (include_delay) {
        stan_code <- paste0(stan_code,
                            "prob[",t,", 2] = inv_logit(logit(a",t,"2) + delay_pos * delay[,",t,"] + bpos[",b_index[t],"] * RE);
                                ")
      } else {
        stan_code <- paste0(stan_code, 
                            "prob[",t,", 2] = inv_logit(logit(a",t,"2) + bpos[",b_index[t],"] * RE);
                     ")
      } 
    } 
    
    else {
      stan_code <- paste0(stan_code, 
                          "prob[",t,", 1] = inv_logit(logit(1 - a",t,"1) + RE);
                    ") 
      if (include_delay) {
        stan_code <- paste0(stan_code, 
                            "prob[",t,", 2] = inv_logit(logit(a",t,"2) + delay_pos * delay[,",t,"] + RE);
          ")
      } else {
        stan_code <- paste0(stan_code, 
                            "prob[",t,", 2] = inv_logit(logit(a",t,"2) + RE);
                    ") 
      } 
    }
  }
  
  stan_code <- paste0(stan_code, "      
     }
  
  model {
    real ps[2];
    int pos;
    vector[N] log_lik; 
    ")
  
  for (t in 1:T){ 
    stan_code <- paste0(stan_code, 
                        "a",t,"1 ~ beta(10, 1);
                           a",t,"2 ~ beta(1, 1);
                          ")
  }
  
  if (add_dependency) {
  stan_code <- paste0(stan_code, 
                      "bpos ~ gamma(1,1);
                             bneg ~ gamma(1,1);
                          ")
  }
  
  if (include_time) {
    stan_code <- paste0(stan_code, "
    prev_amplitude ~ beta(1, 1);
    prev_baseline ~ beta(1, 1);
    mean_gaussian ~ uniform(0, 52);
    sigma_gaussian ~ normal(0, 10);
    ")
  } else {  stan_code <- paste0(stan_code, "
  prev ~ beta(1, 1);
  ")
  }
  
  if (include_delay) {
    stan_code <- paste0(stan_code, "delay_pos ~ normal(-0.5, 10);\n
    ")
  }
  
  stan_code <- paste0(stan_code, "
    RE ~ normal(0, 1);

    pos = 1;
    
  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) + binomial_lpmf(segment(y, pos, s[n]) | 1, prob[segment(tt, pos, s[n]),k,n]);
    }
  target += log_sum_exp(ps); 
  pos = pos + s[n];
  log_lik[n] = log_sum_exp(ps);
  
   }
 }
  generated quantities {
    for (i in 1:1) {
      print(\"Generating quantities1\");
    }

    real Se_mean[T];
    real Sp_mean[T];
    real Se_median[T];
    real Sp_median[T];

    for (i in 1:1) {
      print(\"Generating quantities2\");
    }

    for (t in 1:T) {
      print(\"Sp/Se loop - t =\", t);
      Se_mean[t] = mean(prob[t,2,]); // Change to median if needed
      Sp_mean[t] = mean(1 - prob[t,1,]);

      // Median values
      Se_median[t] = quantile(prob[t,2,], 0.5);
      Sp_median[t] = quantile(1 - prob[t,1,], 0.5);
    }
        
  ")
  
  if (include_time) {
    stan_code <- paste0(stan_code, "
  real weekly_prevalence[52];
  
  // Calculate weekly prevalence using Gaussian seasonality:
  for (w in 1:52) {
    real seasonal_effect = prev_amplitude * exp(-0.5 * square((w - mean_gaussian) / sigma_gaussian));
    real prevalence = prev_baseline + seasonal_effect;
    weekly_prevalence[w] = fmin(fmax(prevalence, 0.0), 1.0);
  }
  ")
  }
  
  stan_code <- paste0(stan_code, "
  }
  ")
  
  return(stan_code)
}

