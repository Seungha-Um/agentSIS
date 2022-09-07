#' Bootstrap Particle Filter for agent-based SIS model 
#' @param y a vector of length (T+1)
#' @param model_config a list with the following components:
#' \itemize{
#'  \item N : the number of agents
#'  \item alpha0 : N x 1 vector of initial infection rates for each agent
#'  \item lambda : N x 1 vector of infection rates for each agent
#'  \item gamma : N x 1 vector of recovery rates for each agent
#'  \item rho : reporting rate 
#'  }
#' @param num_particles The number of particles for BPF
#' @return Returns a list with the following components:
#' \itemize{
#'  \item log_incremental_likelihood : T x 1 vector of estimated loglikelihood for each time step
#'  \item log_final_likelihood : The estimated loglikelihood 
#'  \item weights : The weight for each particle
#'  \item particles : (N x T x P) list of the estimated agent states for each particles 
#'  \item ancestor_mat : (T x P) matrix of indices for sampled particles}
#' @export

sis_bpf <- function(y, model_config, num_particles){
  num_observations <- length(y)
  loglikelihood <- rep(-Inf, num_observations)
  particles <- array(NA, dim = c(model_config$N, num_observations, num_particles));
  ancestor_mat <- matrix(NA, nrow=num_observations, ncol = num_particles) 
  logweights <- logw <- rep(0, num_particles)
  logW <- rep(log(1 / num_particles), num_particles);
  ancestors <- c(1:num_particles)
  alphats <- matrix(model_config$alpha0, nrow = model_config$N, ncol = num_particles);
  xts <- (runif(num_particles * model_config$N) < alphats);
  its <- colSums(xts);
  logweights <- dbinom(x = as.numeric(y[1]), size = its, prob = model_config$rho, log = TRUE);
  loglikelihood[1] <- log(mean(exp(logweights)));
  logw <- logW + logweights;
  weights <- lw.normalize(logw);
  particles[,1,] <- xts
  ancestor_mat[1,] <- ancestors
  for (t in 2:num_observations){
    alpha_t <- apply(xts, 2, function(x) sis_get_alpha(x, model_config = model_config))
    xts_new <- matrix(NA, nrow = model_config$N, ncol = num_particles)
    xts_new <- runif(model_config$N*num_particles) < alpha_t
    xts <- xts_new; 
    its <- colSums(xts);
    logW <- rep(-log(num_particles), num_particles);
    logw <- rep(0, num_particles);
    logweights <- dbinom(x = as.numeric(y[t]), size = its, prob = model_config$rho, log = TRUE);
    loglikelihood[t] <- log(mean(exp(logweights)));
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ancestors <- sample.int(n = num_particles, prob = weights, replace = TRUE);
    xts[,] <- xts[,ancestors];
    ancestor_mat[t,] <- ancestors
    particles[,t,] <- xts
  }
  result <- list(log_incremental_likelihood = loglikelihood, 
                 log_final_likelihood = sum(loglikelihood),
                 weights=weights,
                 particles = particles, ancestor_mat=ancestor_mat)
  return(result)
}


