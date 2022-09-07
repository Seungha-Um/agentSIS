#'  Create hyperparameter object for PMCMC in agent-base SIS model
#'
#'  @param mu_alpha0 mean parameter for prior on alpha intercept ; if not provided, defaults to 0
#'  @param mu_alpha1 mean parameter for prior on alpha slope ; if not provided, defaults to 0
#'  @param mu_lam0 mean parameter for prior on lambda intercept ; if not provided, defaults to 0
#'  @param mu_lam1 mean parameter for prior on lambda slope ; if not provided, defaults to 0
#'  @param mu_gam0 mean parameter for prior on gamma intercept ; if not provided, defaults to 0
#'  @param mu_gam1 mean parameter for prior on gamma slope ; if not provided, defaults to 0
#'  @param mu_rho mean parameter for prior on rho ; if not provided, defaults to 0.7
#'  @param sigma_alpha0 shape prameter for prior on alpha intercept; if not provided, defaults to 3
#'  @param sigma_alpha1 shape prameter for prior on alpha slope; if not provided, defaults to 3
#'  @param sigma_lam0 shape prameter for prior on lambda intercept; if not provided, defaults to 3
#'  @param sigma_lam1 shape prameter for prior on lambda slope; if not provided, defaults to 3
#'  @param sigma_gam0 shape prameter for prior on gamma intercept; if not provided, defaults to 3
#'  @param sigma_gam1 shape prameter for prior on gamma slope; if not provided, defaults to 3
#'  @param sigma_rho shape prameter for prior on rho; if not provided, defaults to 3
#'  @param lower_alpha0 lower bound for truncated Normal prior on alpha intercept
#'  @param lower_alpha1 lower bound for truncated Normal prior on alpha slope
#'  @param lower_lam0 lower bound for truncated Normal prior on lambda intercept
#'  @param lower_lam1 lower bound for truncated Normal prior on lambda slope
#'  @param lower_gam0 lower bound for truncated Normal prior on gamma intercept
#'  @param lower_gam1 lower bound for truncated Normal prior on gamma slope
#'  @param lower_rho lower bound for truncated Normal prior on rho
#'  @param upper_alpha0 upper bound for truncated Normal prior on alpha intercept
#'  @param upper_alpha1 upper bound for truncated Normal prior on alpha slope
#'  @param upper_lam0 upper bound for truncated Normal prior on lambda intercept
#'  @param upper_lam1 upper bound for truncated Normal prior on lambda slope
#'  @param upper_gam0 upper bound for truncated Normal prior on gamma intercept
#'  @param upper_gam1 upper bound for truncated Normal prior on gamma slope
#'  @param upper_rho upper bound for truncated Normal prior on rho
#'  @return Returns a list containing the function arguments.
Hypers <- function(mu_alpha0 = 0, mu_alpha1 = 0, mu_lam0 = 0, mu_lam1 = 0,
                   mu_gam0 = 0, mu_gam1 = 0, mu_rho = 0.7, sigma_alpha0 = 3, sigma_alpha1 = 3,
                   sigma_lam0 = 3, sigma_lam1 = 3, sigma_gam0 = 3, sigma_gam1 = 3, sigma_rho = 1,
                   lower_alpha0 = -Inf, lower_alpha1 = -Inf, lower_lam0 = -Inf, lower_lam1 = -Inf,
                   lower_gam0 = -Inf, lower_gam1 = -Inf, lower_rho = -Inf, upper_alpha0 = Inf,
                   upper_alpha1 = Inf, upper_lam0 = Inf, upper_lam1 = Inf, upper_gam0 = Inf,
                   upper_gam1 = Inf, upper_rho =Inf){
  out <- list()
  out$mu<-  c(mu_alpha0, mu_alpha1, mu_lam0, mu_lam1, mu_gam0, mu_gam1 , mu_rho)
  out$sigma = c(sigma_alpha0, sigma_alpha1, sigma_lam0, sigma_lam1, sigma_gam0, sigma_gam1, sigma_rho)
  out$lower = c(lower_alpha0, lower_alpha1, lower_lam0, lower_lam1, lower_gam0, lower_gam1, lower_rho)
  out$upper = c(upper_alpha0, upper_alpha1, upper_lam0, upper_lam1, upper_gam0, upper_gam1, upper_rho)
  return(out)
}


#'
#' create a list which provides the parameters for implementing Markov chian.
#'
#' @param num_burn The number of burn-in iterations for the chain
#' @param num_save The number of samples to collect ; num_burn + num_save iterations are run
#' @param num_print The interval at which the chain's progress is printed
#' @param step_size The step size for MHRW
#' @param update_gamma If true, gamma (recovery rate) is updated using (Truncated) Normal prior
#' @param update_lambda If true, lambda (infection rate) is updated using (Truncated) Normal prior
#' @param num_particle The number of particles in the PF algorithm
#'
#' @return Returns a list containing the function arguments.

Opt <- function(num_burn = 5000, num_save = 5000, num_print = 100, step_size = 0.2,
                update_lambda = TRUE, update_gamma = TRUE, fully_network = TRUE,
                num_particle = 100){
  out <- list()
  out$num_burn <- num_burn
  out$num_save <- num_save
  out$num_print <- num_print
  out$step_size <- step_size
  out$update_lambda <- update_lambda
  out$update_gamma <- update_gamma
  out$fully_network <- fully_network
  out$num_particle <- num_particle
  return(out)
}
#' Implement the PMCMC for the agent-based SIS model
#'
#' Runs the Markov chain for the agent-based SIS model and collects the output
#'
#' @param X N x (1+P) matrix of agents characteristic including the intercept column
#' @param y T x 1 vector of observed infected counts
#' @param hypers List of hyperparameter values obtained from Hypers function
#' @param opts List of MCMC chain settings obatined from Opts function
#' @param ini_param a set of parameters for initial values
#' @param fixed_gam If true, the recovery rate will be updated, otherwise, specify the fixed recovery rate (scalar)
#'
#' @return Returns a list with the following components:
#' \itemize{
#'   \item actual_counts : The estimated actual counts for each iteration of the chain
#'   \item actual_counts_mean : The estimated actual counts averaged over iterations
#'   \item obs_counts : The estimated reported counts for each iteration of the chain
#'   \item obs_counts_mean : The estimated reported counts averaged over iterations
#'   \item param_mat : matrix of parameter samples
#'   \item alpha0 : posterior samples of alpha0
#'   \item alpha1 : posterior samples of alpha1
#'   \item lam0 : posterior samples of lambda0
#'   \item lam1 : posterior samples of lambda1
#'   \item gam0 : posterior samples of gamma0
#'   \item gam1 : posterior samples of gamma1
#'   \item rho : posterior samples of rho
#'   \item accept_chain : acceptance indicator for MH
#'   \item acceptance_rate : acceptance rate
#' }

PMCMC_SIS <- function(X, y, opt = Opt(), hypers = NULL, ini_param = NULL, fixed_gam = NULL){

  if(is.null(hypers)){
    hypers <- Hypers()
  }

  N <- ncol(X)
  if(is.null(ini_param)) ini_param = truncnorm::rtruncnorm(1, a = hypers$lower,
                                                               b = hypers$upper, mean = hypers$mu, sd = hypers$sigma)

  if(is.null(ini_param))  ini_param[1:2] <- c(-log((N/y[1]-1)), 0)

  ini_alpha <- fun_rate(ini_param[1], ini_param[2])
  ini_lam <- fun_rate(ini_param[3], ini_param[4])

  if(is.null(fixed_gam)) {ini_gam <- fun_rate(ini_param[5], ini_param[6])
  } else { ini_gam <- rep(fixed_gam, N)}

  ini_rho <- expit(ini_param[7])
  curr_param <- ini_param
  curr_param[7] <- ini_rho

  num_particles <- opt$num_particle

  config <- list(N = N, alpha0 = ini_alpha, lambda = ini_lam,
                 gamma = ini_gam, rho = ini_rho)

  step_size <- opt$step_size

  burn_in <- opt$num_burn
  num_mcmc <- opt$num_save + burn_in;

  chain <- matrix(nrow = num_mcmc, ncol = 7)
  accept_chain <- rep(0, num_mcmc);
  curr_config <- config
  curr_lpost <- sis_bpf(y, curr_config, num_particles)$log_final_likelihood + ljacobian(curr_param);

  if(is.null(fixed_gam)) { lprior  <- lprior_fun(curr_param, hypers)
  } else { lprior <- lprior_fix_gam(curr_param, hypers) }

  curr_lpost <- curr_lpost + lprior
  sel_array <- matrix(NA,  num_mcmc, length(y))
  selected <- matrix(nrow = N, ncol = length(y))

  for(imcmc in 1 : num_mcmc){
    ## propose
    prop_param <- qkernel(curr_param, step_size);

    if(is.null(fixed_gam)) { prop_config <- update_config(curr_config, prop_param)
    } else { prop_config <- update_config_fix_gam(curr_config, prop_param) }

    pf_result <- sis_bpf(y, prop_config, num_particles)
    prop_lpost <-  pf_result$log_final_likelihood + ljacobian(prop_param);

    if(is.null(fixed_gam)) { lprior  <- lprior_fun(prop_param, hypers)
    } else { lprior <- lprior_fix_gam(prop_param, hypers) }

    prop_lpost <- prop_lpost + lprior

    ind <- sample(num_particles, 1,  prob = pf_result$weight)
    selected_new <- pf_result$particles[, , ind]

    if (log(runif(1)) < prop_lpost - curr_lpost){
      accept_chain[imcmc] <- 1;
      curr_param <- prop_param;
      curr_lpost <- prop_lpost;
      selected <- selected_new
    }

    if(is.null(fixed_gam)) { chain[imcmc, ] <- curr_param
    } else { chain[imcmc, c(1:4,7)] <- curr_param[c(1:4,7)] }

    sel_array[imcmc,] <- colSums(selected)

    if(imcmc %% opt$num_print == 0){
      cat("parameter est : ",round(chain[imcmc, ], 3), "( iteration:",imcmc,")", "\n");
    }
  }

  temp <- t(sapply(burn_in:num_mcmc, function(x) sel_array[x,] * chain[x, 7]))

  out <- list()
  out$actual_counts <- sel_array[burn_in:num_mcmc,]
  out$actual_counts_mean <- colMeans(sel_array[burn_in:num_mcmc,])
  out$obs_counts <- temp
  out$obs_counts_mean <- colMeans(temp)
  out$param_mat <- chain[burn_in:num_mcmc, ]
  out$alpha0 <- chain[burn_in:num_mcmc, 1]
  out$alpha1 <- chain[burn_in:num_mcmc, 2]
  out$lam0 <- chain[burn_in:num_mcmc, 3]
  out$lam1 <- chain[burn_in:num_mcmc, 4]
  out$gam0 <- chain[burn_in:num_mcmc, 5]
  out$gam1 <- chain[burn_in:num_mcmc, 6]
  out$rho <- chain[burn_in:num_mcmc, 7]
  out$accept_chain <- accept_chain[burn_in:num_mcmc]
  out$acceptance_rate <- mean(accept_chain[burn_in:num_mcmc])
  return(out)
}





