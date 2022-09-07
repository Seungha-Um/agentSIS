
#' Conditional Bootstrap Particle Filter with ancestor sampling for agent-based SIS model
#' @param y a vector of length (T+1)
#' @param model_config a list with the following components:
#' \itemize{
#'  \item N : the number of agents
#'  \item alpha0 : N x 1 vector of initial infection rates for each agent
#'  \item lambda : N x 1 vector of infection rates for each agent
#'  \item gamma : N x 1 vector of recovery rates for each agent
#'  \item rho : reporting rate
#'  \item fix_agent : a reference trajectory for conditional BPF
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

sis_aspf <- function(y, model_config, num_particles){
  num_observations <- length(y)
  loglikelihood <- rep(-Inf, num_observations)
  particles <- array(NA, dim = c(model_config$N, num_observations, num_particles));
  ancestor_mat <- matrix(NA, nrow=num_observations, ncol = num_particles)
  logweights <- logw <- rep(0, num_particles)
  logW <- rep(log(1 / num_particles), num_particles);
  ancestors <- c(1:num_particles)
  alphats <- matrix(model_config$alpha0, nrow = model_config$N, ncol = num_particles);
  est_agent <- matrix(NA, nrow=N, ncol = num_observations)
  xts <- (runif(num_particles * model_config$N) < alphats);
  xts[,num_particles] <- model_config$fix_agent[,1]
  alphats <- apply(xts, MARGIN = 2, FUN = function(xx) sis_get_alpha(agent_state = xx, model_config = model_config));
  its <- colSums(xts)
  logweights <- dbinom(x = as.numeric(y[1]), size = its, prob = model_config$rho, log = TRUE);
  loglikelihood[1] <- log(mean(exp(logweights)));
  logw <- logW + logweights;
  weights <- lw.normalize(logw);
  particles[,1,] <- xts
  est_agent[,1] <- colSums(apply(xts, MARGIN=1, FUN = function(x) x * weights))

  for (t in 2 : num_observations){
    sample_p <- logweights + colSums(apply(alphats, 2, function(x) dbinom(model_config$fix_agent[,t], size=1, x, log=TRUE)))
    AS_ind <- sample.int(n = num_particles, size =1 ,prob = lw.normalize(sample_p), replace = TRUE);
    ancestors <- sample.int(n = num_particles, prob = weights, replace = TRUE);
    ancestors[num_particles] <- AS_ind;
    ancestor_mat[t,] <- ancestors
    xts[,] <- xts[,ancestors];
    particles[,t-1,] <- xts
    alphats <- apply(xts, MARGIN = 2, FUN = function(xx) sis_get_alpha(agent_state = xx, model_config = model_config));
    xts <- (runif((num_particles-1) * model_config$N) < alphats[,-num_particles]);
    xts <- cbind(xts,model_config$fix_agent[,t])
    its <- colSums(xts)
    logW <- rep(-log(num_particles), num_particles);
    logw <- rep(0, num_particles);
    logweights <- dbinom(x = as.numeric(y[t]), size = its, prob = model_config$rho, log = TRUE);
    loglikelihood[t] <- log(mean(exp(logweights)));
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    #est_agent[,t] <- colSums(apply(xts, MARGIN=1, FUN = function(x) x * weights))
    if(t == num_observations) particles[,t,] <- xts
  }
  result <- list(log_incremental_likelihood = loglikelihood,
                 log_final_likelihood = sum(loglikelihood),
                 weights=weights, particles = particles, ancestor_mat=ancestor_mat)
  return(result)
}
