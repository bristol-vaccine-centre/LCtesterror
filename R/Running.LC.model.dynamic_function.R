
#' @title Runs model to infer disease prevalence from test data accounting for unknown test error.
#' @description Runs bayesian LC stan model to infer true disease prevalence from test data accounting for unknown test error.
#' Options to specify the number of tests, to include 'time' or 'delay' until testing as covariates, and to specify dependency relationships between tests.
#' Model compilation (which occurs after the code is printed) can take a few minutes.
#'
#' @param data A dataframe of binary test results with a column for each test and row for each individual. If delay until test covariates are to be included in the model, a column with delay data for each test should be included after all of the test result data columns, in the same test order as the test result data.
#' If time is included, time data should be included as the last column.
#' @param num_tests A numeric value for the number of tests.
#' @param test_names_defined A character vector of identifiers / names for each test, given in the same order as the data. Default = NULL.
#' @param data_ID An optional name identifier for the dataset used. Default = NULL.
#' @param dependency_groups A list of vectors specifying dependencies between tests. All test numbers that are not independent of each other are specified in a single vector. Test numbers refer to the order they are specified in the data columns.
#' Example for scenario where tests 1+2 are dependent on each other and where tests 3+4 are dependent on each other: "list(c(1, 2), c(3, 4))". Default = NULL.
#' @param prior_spec Specification of specificity prior. A list of length equal to the value of num_tests with each element containing a vector of length two specifying the alpha and beta parameters for the Beta prior. Default = c(10, 1) for each test.
#' @param prior_sens Specification of sensitivity prior. A list of length equal to the value of num_tests with each element containing a vector of length two specifying the alpha and beta parameters for the Beta prior. Default = c(1, 1) for each test.
#' @param prior_delay Specification of delay prior, if included in covariates. A vector of length two specifying the mu and sigma parameters for the Normal prior. Default = c(-0.5, 2).
#' @param iter The number of iterations for the stan model. Default = 1000.
#' @param chains The number of chains for the stan model. Default = 4.
#' @param warmup The number of warmup iterations for the stan model. Default = 500.
#' @param stan_arg Optional extra arguments to pass to the rstan::sampling function. Default = NULL.
#' @param covariates Optional vector of covariates to include in the model. Either c("Time") or c("Delay"). If specified, relevant individual-level data must be provided in data. Default = NULL.
#' @param prior_list Optional list of prior distributions-specified to generate in R- to plot with posteriors.
#' @param n_samples If prior_list provided, number of times to samples from the specified prior distributions. Default = 1000.
#' @param model_name Optional name of model. Default = NULL.
#' @param plot_chains For plotting parameter posteriors using rstan::stan_dens, should the density for each chain be plotted? FALSE = for each parameter the draws from all chains are combined. Default = TRUE.
#' @return Stan model fit and various summary outputs:
#'   \describe{
#'   \item{stan_fit}{Fitted stan LC model as returned from rstan::sampling.}
#'   \item{data_inputs}{List of input data used for model fitting.}
#'   \item{test_data_input}{Test data used for model fitting in format given to function.}
#'   \item{stan_fit_summary_df}{A data frame summarising model parameters from rstan::summary.}
#'   \item{prev_mean}{The estimated mean prevalence with 95% credible intervals (in models without time).}
#'   \item{prev_time}{In models with time: A data frame of inferred prevalence estimates over time with 95% credible intervals.}
#'   \item{prev_time_plot}{In models with time: Plotted inferred prevalence over time.}
#'   \item{median_sens}{Median inferred sensitivity for each test with 95% credible intervals.}
#'   \item{median_spec}{Median inferred specificity for each test with 95% credible intervals.}
#'   \item{median_sens_spec_table}{A formatted table summarising each test's inferred sensitivity and specificity with 95% credible intervals .}
#'   \item{traceplots}{Plots traceplots of log posterior, and posterior density plots for prevalence, sensitivity and specifcity parameters.}
#'   \item{pairs.plot_function}{A function to plot a pairs plot for prev, sens and spec parameters.}
#'   \item{mean_Rhat}{The mean R-hat value across parameters.}
#'   \item{mean_ESS}{The mean effective sample size (ESS) across parameters.}
#'   \item{priors.posteriors.plots_function}{Plots prior vs posterior plots if prior_list was provided.}
#'   \item{stan_code}{The generated Stan model code.}
#'   \item{include_time}{Logical; whether time was included in the model.}
#'   \item{include_delay}{Logical; whether delay was included in the model.}
#' }
#' @importFrom rstan summary sampling rstan_options traceplot stan_dens
#' @importFrom dplyr filter group_by mutate rename row_number all_of starts_with n_distinct relocate summarise n select
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom gt gt tab_header
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_density geom_point geom_ribbon labs theme_minimal theme_bw
#' @importFrom stats median na.omit setNames
#' @importFrom graphics pairs
#' @importFrom utils globalVariables
#' @name run.LC.model
#' @examples
#' if (interactive()) {
#' # Example with no test dependence, delay, or time included
#' # Simulate data
#' sim_data <- sim.test.data(sim_size = 100)
#'
#' # Run LC model using simulated test results
#' fit <- run.LC.model(sim_data$test_results, num_tests = 4,
#'                       data_ID = "sim", model_name = "basic_sim",
#'                       dependency_groups = NULL, covariates = NULL,
#'                       iter=1000, chains=4, warmup=500)
#'
#' # View estimated prevalence, sens and spec, and traceplots
#' fit$prev_mean
#' fit$median_sens_spec_table
#' fit$traceplots
#'
#'
#' # Example with test dependence
#' # Simulate data
#' sim_data <- sim.test.data(sim_size = 100)
#'
#' # Run LC model using simulated test results
#' fit <- run.LC.model(sim_data$test_results, num_tests = 4,
#'                       data_ID = "sim", model_name = "dependence_sim",
#'                       dependency_groups = list(c(1, 2), c(3, 4)), covariates = NULL,
#'                       iter=1000, chains=4, warmup=500)
#'
#'
#' # Example with 'delay until test' data
#' # Simulate data
#' sim_data <- sim.test.data(sim_size = 100, delay=TRUE)
#'
#' # Run LC model using simulated test results
#' fit <- run.LC.model(sim_data$test_results, num_tests = 4,
#'                       data_ID = "sim_delay", model_name = "delay_sim",
#'                       dependency_groups = NULL, covariates = c("delay"),
#'                       iter=1000, chains=4, warmup=500)
#'
#'
#' # Example with changing prevalence over time (for a disease with a seasonal peak)
#' # Simulate data
#' sim_data <- sim.test.data.time(sim_size = 10)
#'
#' # Run LC model using simulated test results
#' fit <- run.LC.model(sim_data$test_results, num_tests = 4,
#'                       data_ID = "sim_time", model_name = "time_sim",
#'                       dependency_groups = NULL, covariates = c("time"),
#'                       iter=1000, chains=4, warmup=500)
#'
#' # Plot prev over time
#' plots <- plot.prev.time(fit, sim_time_data = sim_data)
#' plots$weekly_prev_plot
#' }
#'

utils::globalVariables(c("delay", "Time", "time", "week", "CI_min", "CI_max", "mean_prevalence",
                         ".", "2.5%", "97.5%", "Sensitivity", "Specificity",
                         "week", "CI_min", "CI_max", "value", "type", "n"))

#' @export
run.LC.model <- function(data, num_tests, test_names_defined=NULL, data_ID = NULL,
                         dependency_groups = list(),
                         prior_spec = NULL, prior_sens = NULL,
                         prior_delay = NULL,
                         iter=1000, chains=4, warmup=500, stan_arg=list(),
                         covariates = NULL, prior_list = NULL, n_samples=iter,
                         model_name=NULL, plot_chains=TRUE) {

  # Use tryCatch to handle errors inside the function
  result <- tryCatch({

    Individuals <- c(1:(nrow(data)))
    data_long <- cbind(Individuals, data)
    test_names <- as.character(seq_len(num_tests))
    message(paste("test_names:", paste(test_names, collapse=", ")))

    if (!is.null(covariates)) {
      col_names <- c("N", test_names, tolower(covariates))  # Include covariate names if present
    } else {
      col_names <- c("N", test_names)
    }

    if ("delay" %in% tolower(covariates)) {
      delay_names <- list()
      for(i in 1:num_tests) {
        delay_name <- paste0("delay_", i)
        delay_names[i] <- delay_name
      }
      col_names <- c("N", test_names, delay_names)
    }

    message(paste("col_names_data:", paste(colnames(data_long), collapse=", ")))
    message(paste("col_names_defined:", paste(col_names, collapse=", ")))
    colnames(data_long) <- col_names
    message(paste("col_names:", paste(colnames(data_long), collapse=", ")))

    numeric_cols <- grep("^[0-9]+$", colnames(data_long), value = TRUE)
    message(paste("numeric_cols:", paste(numeric_cols, collapse=", ")))

    test_result_data_long <- data_long %>%
      pivot_longer(cols = all_of(numeric_cols), names_to = "T", values_to = "y")
    #message("test_result_data_long:")
    #message(test_result_data_long)

    if ("delay" %in% tolower(covariates)) {

      test_result_data_long <- test_result_data_long %>%
        dplyr::select(-starts_with("delay"))
      #message("test_result_data_long (after removing delay columns):")
      #message(test_result_data_long)

      delay_data_long <- data_long %>%
        pivot_longer(cols = all_of(starts_with("delay")), names_to = "T_delay", values_to = "delay") %>%
        dplyr::select(delay)
      #message("delay_data_long:")
      #message(delay_data_long)

      test_data_long <- cbind(test_result_data_long, delay_data_long)
      #message("test_data_long (after merging delay data):")
      #message(test_data_long)

    } else { test_data_long <- test_result_data_long
    #message("test_data_long (no delay data):")
    #message(test_data_long)
    }


    test_data_long <- na.omit(test_data_long)
    #message("test_data_long (after removing NAs):")
    #message(test_data_long)

    tt <- as.numeric(test_data_long$T)
    nn <- as.numeric(test_data_long$N)
    y <- as.numeric(as.character(test_data_long$y))
    T <- as.numeric(length(seq_len(num_tests)))
    O <- nrow(test_data_long)
    N <- n_distinct(test_data_long$N)
    # Matrix of all test combinations
    test_names <- paste0("test", seq_len(num_tests))
    test_list <- setNames(replicate(num_tests, c(0, 1), simplify = FALSE), test_names)
    comb_list <- expand.grid(test_list) # Generate all test combinations
    comb_list <- comb_list[rowSums(comb_list) > 0, ] # Remove rows where all tests are zero
    comb_list$test_count <- rowSums(comb_list) # Create test count column
    # Order by number of tests, and then by priority of test columns
    order_cols <- c("test_count", rev(test_names)) # Prioritise test1, test2, etc.
    comb_list <- comb_list[do.call(order, comb_list[order_cols]), ]
    comb_list <- comb_list[, test_names] # Remove test_count column
    # Convert to matrix
    comb_matrix <- as.matrix(comb_list)
    rownames(comb_matrix) <- NULL
    colnames(comb_matrix) <- NULL
    # Get number of unique test combinations
    C <- as.numeric(nrow(comb_matrix))
    comb_matrix_df <- as.data.frame(comb_matrix)
    rownames(comb_matrix_df) <- NULL
    colnames(comb_matrix_df) <- NULL

    #add covariates if given
    time_points <- NULL
    d <- NULL

    # only check for "time" and "delay" if covariates is not NULL
    if (!is.null(covariates)) {

      # Time - add if specified as a covariate
      if ("time" %in% tolower(covariates)) {
        if (!"time" %in% colnames(test_data_long)) {
          stop("Error: Time column is missing from dataset but time is in covariates.")
        }
        time_points_data <- test_data_long %>%
          group_by(N) %>%
          summarise(time_points = mean(time, na.rm = TRUE))
        time_points <- as.numeric(time_points_data$time_points)
      }

      # Delay - add if specified as a covariate
      if ("delay" %in% tolower(covariates)) {
        if (!"delay" %in% colnames(test_data_long)) {
          stop("Error: `delay` column is missing from dataset but `delay` is in covariates.")
        }
        d_NA <- test_data_long %>%
          dplyr::select(-y) %>%
          tidyr::pivot_wider(names_from = T, values_from = delay, names_prefix = "T_delay") %>%
          dplyr::relocate(N) %>%  # Ensure "N" column stays first
          dplyr::select(N, order(as.numeric(gsub("T_delay", "", colnames(.)[-1])))+1) %>%
          dplyr::select(-N) %>%
          as.matrix()
        d <- apply(d_NA, 2, function(col) replace(col, is.na(col), round(median(col, na.rm = TRUE))))
        message(d)
      }
    }

    #for ragged model
    s <- test_data_long %>% group_by(N) %>% summarise(group = n())
    s <- as.numeric(s$group)

    # Assign test depdency groups - think i can remove?
    test_to_group <- rep(0, num_tests)  # Default to 0 if no group assigned
    for (g in seq_along(dependency_groups)) {
      test_to_group[dependency_groups[[g]]] <- g  # Assign group number
    }


    # Identify valid dependency groups (ignoring single-test groups)
    valid_groups <- Filter(function(g) length(g) > 1, dependency_groups)
    num_valid_groups <- length(valid_groups)
    G <- num_valid_groups


    # Sens/spec priors
    if(is.null(prior_spec)) {
      prior_spec <- lapply(1:T, function(x) c(10, 1))  #Default Beta(10,1)
    }

    if (is.null(prior_sens)) {
      prior_sens <- lapply(1:T, function(x) c(1, 1))  #Default Beta(1,1)
    }

    # Ensure prior_spec and prior_sens are correctly formatted
    if (length(prior_spec) != T || length(prior_sens) != T) {
      stop("prior_spec and prior_sens must be lists of length equal to num_tests")
    }

    # Extract prior alpha and beta parameters
    alpha_spec <- sapply(prior_spec, function(p) p[1])
    beta_spec <- sapply(prior_spec, function(p) p[2])
    alpha_sens <- sapply(prior_sens, function(p) p[1])
    beta_sens <- sapply(prior_sens, function(p) p[2])


    #delay prior
    if (is.null(prior_delay)) {
      mu_delay <- as.numeric(-0.5)
      sigma_delay <- as.numeric(2)
    } else {
      mu_delay <- as.numeric(prior_delay[1])
      sigma_delay <- as.numeric(prior_delay[2])
    }

    #Data for model
    stan_data <- list(T=T, N=N, O=O, tt=tt, y=y, s=s, C=C,
                      comb_matrix=comb_matrix_df,
                      test_to_group = test_to_group, G=G,
                      alpha_spec = alpha_spec, beta_spec = beta_spec,
                      alpha_sens = alpha_sens, beta_sens = beta_sens
                      )
    if (!is.null(time_points)) { #add time points if time data given
      stan_data$time_points <- time_points
    }
    if (!is.null(d)) { #add delay data if delay data is given
      stan_data$delay_days <- d
      stan_data$mu_delay <- mu_delay
      stan_data$sigma_delay <- sigma_delay
    }

    #dynamically generate stan code
    include_time <- FALSE  # Default to FALSE
    include_delay <- FALSE  # Default to FALSE

    if (!is.null(covariates)) {
      if ("time" %in% tolower(covariates)) {
        include_time <- TRUE
      }
      if ("delay" %in% tolower(covariates)) {
        include_delay <- TRUE
      }
    }

    stan_code <- generate.stan.model(num_tests, include_time, include_delay, dependency_groups)
    message(stan_code)
    # Compile Stan model from generated code
    stan_model_compiled <- rstan::stan_model(model_code = stan_code)

    # parallel processing for chains
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())

    #model run
    stan_fit <- rstan::sampling(stan_model_compiled,
                         data = stan_data,
                         iter = iter,
                         seed = 1234,
                         chains = chains,
                         warmup = warmup,
                         !!!stan_arg)  # Pass additional arguments using bang bang bang to unpack list

    #model outputs

    # Model summary df
    fit_summary <- rstan::summary(stan_fit)
    stan_fit_summary <- fit_summary$summary #summary for all chains merged. #se_mean = monte carlo standard errors (from MCMC)
    stan_fit_summary_df <- as.data.frame(stan_fit_summary )

    #Prev
    if("prev" %in% rownames(stan_fit_summary_df)) { #checks exact match

      prev_rows <- grepl("^prev$", rownames(stan_fit_summary_df))
      stan_fit_summary_prev <- subset(stan_fit_summary_df, prev_rows)
      stan_fit_summary_prev_df <- stan_fit_summary_prev %>% dplyr::select(mean, '2.5%', '97.5%')

      #Function to format prevalence:
      format.prevalence <- function(summary_row) {
        # Extract relevant values and convert to percentages
        prevalence <- summary_row["mean"] * 100
        lower_ci <- summary_row["2.5%"] * 100
        upper_ci <- summary_row["97.5%"] * 100
        # Format to one decimal place
        prevalence <- format(round(prevalence, 1), nsmall = 1)
        lower_ci <- format(round(lower_ci, 1), nsmall = 1)
        upper_ci <- format(round(upper_ci, 1), nsmall = 1)

        result <- paste0("Mean prevalence: ", prevalence, "% (", lower_ci, "-", upper_ci, "%)")
        return(result)
      }

      prev_mean <- format.prevalence(stan_fit_summary_prev_df)
    }

    else {
      #Prev over time
      prev_time_rows <- grepl("^weekly_prevalence", rownames(stan_fit_summary_df))
      stan_fit_summary_prev_time <- subset(stan_fit_summary_df, prev_time_rows)
      stan_fit_summary_prev_time_df <- stan_fit_summary_prev_time %>% dplyr::select(mean, '2.5%', '97.5%') %>%
        dplyr::mutate(week = row_number()) %>%
        dplyr::rename(CI_min = '2.5%', CI_max = '97.5%', mean_prevalence = mean)

      prev_time_plot <- stan_fit_summary_prev_time_df %>%
        ggplot(aes(y=mean_prevalence, x=week)) +
        geom_point() +
        geom_ribbon(aes(ymin = CI_min, ymax = CI_max), alpha=0.3, fill="grey") +
        theme_bw()

    }


    #Sens/spec
    sens_rows <- grepl("^Se_median", rownames(stan_fit_summary_df))
    stan_fit_summary_sens_median <- stan_fit_summary_df %>%
      subset(sens_rows) %>%
      dplyr::select(mean, '2.5%', '97.5%')

    spec_rows <- grepl("^Sp_median", rownames(stan_fit_summary_df))
    stan_fit_summary_spec_median <- stan_fit_summary_df %>%
      subset(spec_rows) %>%
      dplyr::select(mean, '2.5%', '97.5%')

    if (is.null(test_names_defined)) {
      test_names <- test_names
    } else {
      test_names <- test_names_defined
    }

    if (is.null(model_name)) {
      model_name <- "unknown" }

    if (is.null(data_ID)) {
      data_ID <- "unknown" }

    #Function to format sens/spec
    sens.spec.table <- function(summary_df, data_name = NULL, test_names, model_fit, model_name=NULL) {

      # Extract rows for sensitivity and specificity
      sens_rows <- grepl("^Se_median", rownames(summary_df))
      spec_rows <- grepl("^Sp_median", rownames(summary_df))
      sens_df <- summary_df[sens_rows, ]
      spec_df <- summary_df[spec_rows, ]
      # Assign row names based on test order
      rownames(sens_df) <- test_names
      rownames(spec_df) <- test_names
      # Formatting
      sens_table <- sens_df %>%
        dplyr::mutate(Sensitivity = sprintf("%.3f (%.3f - %.3f)", mean, `2.5%`, `97.5%`)) %>%
        dplyr::select(Sensitivity)
      spec_table <- spec_df %>%
        dplyr::mutate(Specificity = sprintf("%.3f (%.3f - %.3f)", mean, `2.5%`, `97.5%`)) %>%
        dplyr::select(Specificity)
      combined_table <- cbind(sens_table, spec_table)

      # Extract model info
      iterations <- model_fit@sim$iter
      chains <- dim(stan_fit)[2]
      warmup <- model_fit@sim$warmup

      # Create gt table
      gt_table <- combined_table %>%
        gt::gt(rownames_to_stub = TRUE) %>%
        tab_header(title = "Median sens/spec (mean across chains (95%CI))",
                   subtitle = sprintf("Model: %s; Data: %s; Iterations: %d; Warmup: %d; Chains: %d",
                                      model_name, data_name,
                                      ifelse(exists("iterations"), iterations, NA),
                                      ifelse(exists("warmup"), warmup, NA),
                                      ifelse(exists("chains"), chains, NA))
        )

      return(gt_table)

    }

    sens_spec_table <- sens.spec.table(summary_df = stan_fit_summary_df,
                                       data_name = data_ID, test_names = test_names,
                                       model_fit = stan_fit,
                                       model_name = model_name)

    #Traceplots
    plot.stan.mod <- function(stan_mod, separate_chains = plot_chains) {
      plots <- list()
      available_params <- names(rstan::extract(stan_mod))
      #log posterior (`lp__`)
      if ("lp__" %in% available_params) {
        plots$traceplot_lp_warmup <- rstan::traceplot(stan_mod, pars = "lp__", inc_warmup = TRUE)
        plots$traceplot_lp <- rstan::traceplot(stan_mod, pars = "lp__")
      }
      # Check if "prev" exists before plotting
      if ("prev" %in% available_params) {
        plots$stan_dens_prev <- rstan::stan_dens(stan_mod, pars = "prev", separate_chains = separate_chains)
      }
      # other parameters to check and plot
      params_to_plot <- c("Se_mean", "Sp_mean", "Se_median", "Sp_median")

      for (param in params_to_plot) {
        if (param %in% available_params) {
          plots[[paste0("stan_dens_", param)]] <- rstan::stan_dens(stan_mod, pars = param, separate_chains = separate_chains)
        }
      }
      return(plots)
    }
    trace_plots <- plot.stan.mod(stan_fit)

    #pairs plot
    pairs_plot <- function() {
      pairs(stan_fit, pars = c("prev", "Se_median", "Sp_median"))
    }

    #mean Rhat and ESS
    mean_Rhat <- mean(stan_fit_summary[, "Rhat"], na.rm = TRUE)
    mean_ESS <- mean(stan_fit_summary[, "n_eff"], na.rm = TRUE)

    #Priors/posteriors plots
    priors.posteriors.plots <- function(stan_fit, n_samples, prior_list = list()) {
      set.seed(123)
      # Extract posteriors
      posterior_samples <- rstan::extract(stan_fit)
      priors <- prior_list
      param_names <- names(priors)
      pp_plots <- list()

      # Loop - dynamically adjusts number of posterior samples to match prior samples
      for (param in param_names) {
        if (param %in% names(posterior_samples)) {  # Check if parameter exists in Stan output
          posterior_param <- posterior_samples[[param]]
          n_post_samples <- length(posterior_param)
          # Resample priors to match posterior sample size
          priors_resampled <- sample(priors[[param]], size = n_post_samples, replace = TRUE)
          # Create data frame for plotting
          df_param <- data.frame(
            value = c(priors_resampled, posterior_param),
            type = rep(c("Prior", "Posterior"), each = n_post_samples))
          # Create density plot
          p <- ggplot(df_param, aes(x = value, fill = type)) +
            geom_density(alpha = 0.5) +
            labs(title = param, x = param, y = "Density") +
            theme_minimal()
          # Store plot in list
          pp_plots[[param]] <- p
        }
      }
      return(pp_plots)
    }

    priors.posteriors.plots_function <- function() {
      pp_plots <- priors.posteriors.plots(stan_fit, n_samples, prior_list)
      if (length(pp_plots) > 0) {
        do.call(grid.arrange, c(pp_plots, ncol = 4))
      } else {
        message("No prior/posterior plots")
      }
    }


    return(list(stan_fit = stan_fit,
                data_inputs = stan_data,
                test_data_input = data,
                stan_fit_summary_df = stan_fit_summary_df,
                prev_mean = if (exists("prev_mean")) prev_mean else NULL,
                prev_time = if (exists("stan_fit_summary_prev_time_df")) stan_fit_summary_prev_time_df else NULL,
                prev_time_plot = if (exists("prev_time_plot")) prev_time_plot else NULL,
                median_sens = stan_fit_summary_sens_median,
                median_spec = stan_fit_summary_spec_median,
                median_sens_spec_table = sens_spec_table,
                traceplots = trace_plots,
                pairs.plot_function = pairs_plot,
                mean_Rhat = mean_Rhat,
                mean_ESS = mean_ESS,
                priors.posteriors.plots_function = priors.posteriors.plots_function,
                stan_code = stan_code,
                include_time = include_time,
                include_delay = include_delay

    ) )

  }, error = function(e) {
    message("! Error in `run.LC.model()`: ", e$message)
    return(list(
      stan_fit = NULL,
      message = paste("Error in run.LC.model:", e$message)
    ))
  })
  return(result)
}

