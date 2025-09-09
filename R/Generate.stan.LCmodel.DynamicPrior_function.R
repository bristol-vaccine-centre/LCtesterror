
#' @title Generates stan model code for different model versions.
#' @description Generates stan model code for different LC model versions.
#' Options to specify the number of tests, to include 'time' or 'delay' until testing as covariates (but not both), and to specify dependency relationships between tests.
#' The prior distributions for specificity and sensitivity can be specified as beta distributions.
#' Function is automatically used in run.LC.model() function
#'
#' @param num_tests Number of tests in data.
#' @param include_time Logical indicating whether time should be included in the model. Default = FALSE.
#' @param include_delay Logical indicating whether delay until testing should be included in the model. Default = FALSE.
#' @param dependency_groups A list of vectors specifying dependencies between tests. All test numbers that are not independent of each other are specified in a single vector. Test numbers refer to the order they are specified in the data columns.
#' Example for scenario where tests 1+2 are dependent on each other and where tests 3+4 are dependent on each other: "list(c(1, 2), c(3, 4))". Default = NULL.
#' @param time_model If covariates includes "Time", which model should be used to infer changing prevalence over time - "gaussian" (for a seasonal peak) or "exponential". Default = "gaussian".
#' @param priors List of any other priors, excluding sens/spec/delay priors (which are specified as stan data), to be specified differently to the defaults. Given as character strings as written for stan. Defaults =
#' list(prev_prior= "beta(1,1)", RE_prior= "normal(0,1), "bpos_prior= "gamma(1,1)", bneg_prior= "gamma(1,1)",
#' gaussian_prev_amplitude_prior= "beta(1,1)", gaussian_prev_baseline_prior= "beta(1,1)", mean_gaussian_prior= "uniform(0,52)", sigma_gaussian_prior = "normal(0,10)",
#' exp_prev_baseline_prior= "beta(1,1)", exp_growth_rate_prior= "gamma(1,5)")
#' @return Stan LC model code.
#' @importFrom magrittr %>%
#' @importFrom gridExtra grid.arrange
#' @name generate.stan.model
#' @examples
#' if (interactive()) {
#' # Basic model with 4 tests:
#' generate.stan.model(num_tests = 4)
#'
#' # Model with 4 tests that have a dependence structure
#' # (in this example tests 1+2 are correlated and tests 3+4 are correlated):
#' generate.stan.model(num_tests = 4, dependency_groups = list(c(1, 2), c(3, 4)))
#'
#' # Model with 4 tests and delay or time modelled as a covariate:
#' generate.stan.model(num_tests = 4, include_delay = TRUE)
#' generate.stan.model(num_tests = 4, include_time = TRUE)
#'
#' # Model with 4 tests that have a dependence structure
#' # (in this example tests 1 is independent and tests 3+4 are correlated)
#' # and have either delay or time covariates:
#' generate.stan.model(num_tests = 4, dependency_groups = list(c(1), c(2, 3, 4)), include_delay = TRUE)
#' generate.stan.model(num_tests = 4, dependency_groups = list(c(1), c(2, 3, 4)), include_time = TRUE)
#' }
#'

#' @export generate.stan.model
generate.stan.model <- function(num_tests, include_time = FALSE, include_delay = FALSE, dependency_groups = list(),
                                time_model = "gaussian", priors = list()) {

  # Identify valid dependency groups (ignoring single-test groups)
  valid_groups <- Filter(function(g) length(g) > 1, dependency_groups)
  num_valid_groups <- length(valid_groups)
  G <- num_valid_groups
  time_exp <- time_model == "exponential"
  message(paste("valid_groups:", paste(valid_groups, collapse=", ")))
  message(paste("num_valid_groups:", paste(num_valid_groups, collapse=", ")))
  message(paste("G:", paste(G, collapse=", ")))
  #b_pos/b_neg index vector
  b_index <- rep(0, num_tests)  # Default 0 (no bpos)
  group_counter <- 1
  for (g in valid_groups) {  # Iterate over valid dependency groups
    for (t in g) {
      b_index[t] <- group_counter  # Assign bpos group number
    }
    group_counter <- group_counter + 1
  }
  message(paste("b_index:", paste(b_index, collapse=", ")))
  # Determine if dependency terms (bpos/bneg) should be added
  all_tests_separate <- length(dependency_groups) == num_tests && all(sapply(dependency_groups, length) == 1) #are all tests in seperate groups? If so add dependency will = FALSE
  add_dependency <- !(num_valid_groups == 0 ||
                        (length(dependency_groups) == 1 && length(dependency_groups[[1]]) == num_tests) ||
                        all_tests_separate)
  message(paste("add_dependency:", paste(add_dependency, collapse=", ")))
  # Create a vector: 1 if test has a valid dependency group, 0 otherwise
  dependency_group_existence <- rep(0, num_tests)
  for (group in valid_groups) {
    dependency_group_existence[group] <- 1
  }
  message(paste(" dependency_group_existence:", paste( dependency_group_existence, collapse=", ")))

  T <- as.numeric(num_tests)

  # priors
  if (is.null(priors[["bpos_prior"]])) priors[["bpos_prior"]] <- "gamma(1,1)"
  if (is.null(priors[["bneg_prior"]])) priors[["bneg_prior"]] <- "gamma(1,1)"
  if (is.null(priors[["exp_prev_baseline_prior"]])) priors[["exp_prev_baseline_prior"]] <- "beta(1, 1)"
  if (is.null(priors[["exp_growth_rate_prior"]])) priors[["exp_growth_rate_prior"]] <- "gamma(1, 5)"
  if (is.null(priors[["gaussian_prev_amplitude_prior"]])) priors[["gaussian_prev_amplitude_prior"]] <- "beta(1, 1)"
  if (is.null(priors[["gaussian_prev_baseline_prior"]])) priors[["gaussian_prev_baseline_prior"]] <- "beta(1, 1)"
  if (is.null(priors[["mean_gaussian_prior"]])) priors[["mean_gaussian_prior"]] <- "uniform(0, 52)"
  if (is.null(priors[["sigma_gaussian_prior"]])) priors[["sigma_gaussian_prior"]] <- "normal(0, 10)"
  if (is.null(priors[["prev_prior"]])) priors[["prev_prior"]] <- "beta(1, 1)"
  if (is.null(priors[["RE_prior"]])) priors[["RE_prior"]] <- "normal(0, 1)"



  # Start building Stan code
  stan_code <- "
  data {
    int<lower=1> T; // Number of different tests
    int<lower=1> N; // Number of individuals
    int<lower=1> O; // Number of observations
    array[O] int<lower=1, upper=T> tt; // Test index for each observation
    array[O] int<lower=0, upper=1> y; // Test results (0/1)
    int s[N]; // Number of tests per individual

  // Dynamic prior parameters for sens/spec
  real<lower=0> alpha_spec[T];
  real<lower=0> beta_spec[T];
  real<lower=0> alpha_sens[T];
  real<lower=0> beta_sens[T];
    "
  if (add_dependency) {
    stan_code <- paste0(stan_code,
                        "array[T] int test_to_group;  //Dependency group mapping for each test - think i can remove
    int<lower=1> G; //Number of dependency groups
                        ")
  }
  if (include_time) {
    stan_code <- paste0(stan_code, "vector[N] time_points;\n
                        int<lower=1> max_time;
                        ")
  }

  if (include_delay) {
    stan_code <- paste0(stan_code, "matrix[N, T] delay_days;\n
                        real mu_delay;
                        real sigma_delay;
                        ")
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
    if (time_exp) {  stan_code <- paste0(stan_code, "
         real<lower=0,upper=1> prev_baseline;
         real<lower=0> growth_rate;
      ")
    } else {  stan_code <- paste0(stan_code, "
    real<lower=0,upper=1> prev_amplitude;
    real<lower=0,upper=1> prev_baseline;
    real<lower=0> sigma_gaussian;
    real<lower=0, upper=52> mean_gaussian;
    ")
  }
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
    vector[N] prob[T, 2];

 ")

  if (include_time) {
    if (time_exp) {  stan_code <- paste0(stan_code, "

    simplex[2] theta[N];

     for (n in 1:N) {
      real prevalence = prev_baseline * exp(growth_rate * time_points[n]);
      prevalence = fmin(fmax(prevalence, 0.0), 1.0);

    theta[n][1] = 1 - prevalence;
    theta[n][2] = prevalence;
   }
    ")
    } else {
    stan_code <- paste0(stan_code, "

     simplex[2] theta[N];

  for (n in 1:N) {
    real seasonal_effect = prev_amplitude * exp(-0.5 * square((time_points[n] - mean_gaussian) / sigma_gaussian));
    real prevalence = prev_baseline + seasonal_effect;

    prevalence = fmin(fmax(prevalence, 0.0), 1.0);

    theta[n][1] = 1 - prevalence;
    theta[n][2] = prevalence;
   }
  ")
  }
    } else {
    stan_code <- paste0(stan_code, "
    simplex[2] theta;
    theta[1] = 1 - prev;
    theta[2] = prev;

     ")
  }


  if (include_delay) {
    stan_code <- paste0(stan_code, "matrix[N, T] delay;
  for (t in 1:T) {
    vector[N] col = delay_unscaled[, t];
    delay[, t] = (col - mean(col)) / sd(col);
  }
  ")
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

 #priors
  for (t in 1:T){
    stan_code <- paste0(stan_code, "
  target += beta_lpdf(a",t,"1 | alpha_spec[",t,"], beta_spec[",t,"]);
  target += beta_lpdf(a",t,"2 | alpha_sens[",t,"], beta_sens[",t,"]);
                          ")
  }

  if (add_dependency) {
  stan_code <- paste0(stan_code,
                      "bpos ~ ",priors[["bpos_prior"]],";
                             bneg ~ ",priors[["bneg_prior"]],";
                          ")
  }

  if (include_time) {
    if(time_exp) { stan_code <- paste0(stan_code, "
     prev_baseline ~ ",priors[["exp_prev_baseline_prior"]],";
     growth_rate ~ ",priors[["exp_growth_rate_prior"]],";
    ")
  } else {
    stan_code <- paste0(stan_code, "
    prev_amplitude ~ ",priors[["gaussian_prev_amplitude_prior"]],";
    prev_baseline ~ ",priors[["gaussian_prev_baseline_prior"]],";
    mean_gaussian ~ ",priors[["mean_gaussian_prior"]],";
    sigma_gaussian ~ ",priors[["sigma_gaussian_prior"]],";
    ")
   }
    } else {  stan_code <- paste0(stan_code, "
  prev ~ ",priors[["prev_prior"]],";
  ")
   }


  if (include_delay) {
    stan_code <- paste0(stan_code, "
      target += normal_lpdf(delay_pos | mu_delay, sigma_delay);
    ")
  }

  stan_code <- paste0(stan_code, "
    RE ~ ",priors[["RE_prior"]],";

    pos = 1;

  ")

  if (include_time) {
    stan_code <- paste0(stan_code, "

    for(n in 1:N){
      for(k in 1:2){

      ps[k] = log(theta[n][k]) + binomial_lpmf(segment(y, pos, s[n]) | 1, prob[segment(tt, pos, s[n]),k,n]);

    }

  target += log_sum_exp(ps); //calculates log of sum of exponentials of ps values for each latent class.
  pos = pos + s[n];

  }
     ")
  } else {
  stan_code <- paste0(stan_code, "

  for(n in 1:N){
    for(k in 1:2){
      ps[k] = log(theta[k]) + binomial_lpmf(segment(y, pos, s[n]) | 1, prob[segment(tt, pos, s[n]),k,n]);
    }
  target += log_sum_exp(ps); //calculates log of sum of exponentials of ps values for each latent class.
  pos = pos + s[n];
  log_lik[n] = log_sum_exp(ps);

  }

    ")
  }

 stan_code <- paste0(stan_code, "
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
      Se_mean[t] = mean(prob[t,2,]);
      Sp_mean[t] = mean(1 - prob[t,1,]);

      // Median values
      Se_median[t] = quantile(prob[t,2,], 0.5);
      Sp_median[t] = quantile(1 - prob[t,1,], 0.5);
    }

  ")

  if (include_time) {
    if(time_exp) { stan_code <- paste0(stan_code, "
      real weekly_prevalence[max_time];
        // Calculate weekly prevalence using exponential growth:
   for (w in 1:max_time) {
      real prevalence = prev_baseline * exp(growth_rate * w);
      weekly_prevalence[w] = fmin(fmax(prevalence, 0.0), 1.0);
   }
      ")
    } else {
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
  }

  stan_code <- paste0(stan_code, "
  }
  ")

  return(stan_code)
 }

