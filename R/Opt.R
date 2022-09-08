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
