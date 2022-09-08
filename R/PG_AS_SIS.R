#' Implement the Particle Gibbs with ancestor sampling for the agent-based SIS model
#'
#' Runs the Markov chain for the agent-based SIS model and collects the output
#'
#' @param X N x (1+P) matrix of agents characteristic including the intercept column
#' @param y T x 1 vector of observed infected counts
#' @param hypers List of hyperparameter values obtained from Hypers function
#' @param opts List of MCMC chain settings obatined from Opts function
#' @param ini_param a set of parameters for initial values; c(alpha0, alpha1, lambda0, lambda1, gamma0, gamma1,rho)
#' @param fixed_gam If true, the recovery rate will be updated, otherwise, specify the fixed recovery rate (scalar)
#' @param rho_alpha The value of hyperparameter alpha for beta prior for reporting rate (rho) ; if not provided, defaults to 10
#' @param rho_beta The value of hyperparameter beta for beta prior for reporting rate (rho) ; if not provided, defaults to 2
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


PG_AS_SIS <- function(X, y, opt = Opt(), hypers = NULL, ini_param = NULL, fixed_gam = NULL,
                   rho_alpha = 10, rho_beta = 5){

  if(is.null(hypers)){
    hypers <- Hypers()
  }

  N <- ncol(X)
  num_obs <- length(y)

  if(is.null(ini_param)) ini_param = truncnorm::rtruncnorm(1, a = hypers$lower,
                                                           b = hypers$upper, mean = hypers$mu, sd = hypers$sigma)

  if(is.null(ini_param))  ini_param[1:2] <- c(-log((N/y[1]-1)), 0)

  ini_alpha <- fun_rate(ini_param[1], ini_param[2])
  ini_lam <- fun_rate(ini_param[3], ini_param[4])

  if(is.null(fixed_gam)) {ini_gam <- fun_rate(ini_param[5], ini_param[6])
  } else { ini_gam <- rep(fixed_gam, N)}

  ini_rho <- rbeta(1, rho_alpha, rho_beta)
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

  sel_array <- matrix(NA,  num_mcmc, length(y))
  selected <- matrix(nrow = N, ncol = length(y))

  # sample a initial reference trajectory from BPF
  BPF_result <- BPF(y, curr_config, num_particles);
  ind <- sample(num_particles, 1,  prob = BPF_result$weight)
  path_ind <- rep(NA, num_obs)
  path_ind[num_obs] <- ind

  for(j in num_obs:2){
    path_ind[j-1] <- BPF_result$ancestor_mat[j, ind];
    ind <- path_ind[j-1];
  }

  selected <- sapply(1:num_obs, function(i) BPF_result$particles[,i, path_ind[i]])

  for(imcmc in 1 : num_mcmc){
    config$fix_agent <- selected
    CBPF_AS_result <- CBPF_AS(y, config, num_particles);

    # Sample a particle proportional to weights
    ind <- sample(num_particles, 1,  prob = CBPF_AS_result$weight)
    selected <- CBPF_AS_result$particles[, , ind]

    alpha_est <- apply(selected, MARGIN = 2, FUN = function(xx) sis_get_alpha(agent_state = xx, model_config = curr_config));
    curr_lpost_beta <- sum(dbinom(selected[,1], size=1, curr_config$alpha0, log=TRUE)) +
      sum(sapply(2:length(y), function(tt) dbinom(selected[,tt], size=1, alpha_est[,tt-1],log=TRUE))) +
      lprior_PG(curr_param, hypers)

    prop_param <- qkernel(curr_param, step_size);

    if(is.null(fixed_gam)) { prop_config <- update_config(curr_config, prop_param)
    } else { prop_config <- update_config_fix_gam(curr_config, prop_param) }

    alpha_est_prop <- apply(selected, MARGIN = 2, FUN = function(xx) sis_get_alpha(agent_state = xx, model_config = prop_config));

    prop_lpost_beta <- sum(dbinom(selected[,1], size=1, prop_config$alpha0, log=TRUE)) +
      sum(sapply(2:length(y), function(tt) dbinom(selected[,tt], size=1, alpha_est_prop[,tt-1], log=TRUE))) +
      lprior_PG(prop_param, hypers)

    if (log(runif(1)) < prop_lpost_beta - curr_lpost_beta){
      accept_chain[imcmc] <- 1;
      curr_param[1:6] <- prop_param[1:6];
    }

    curr_param[7] <- rbeta(1, sum(y) + rho_alpha, sum(colSums(selected) - y) + rho_beta )

    if(is.null(fixed_gam)) { curr_config <- update_config(curr_config, curr_param)
    } else { curr_config <- update_config_fix_gam(curr_config, curr_param) }

    ## save the parameters
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


