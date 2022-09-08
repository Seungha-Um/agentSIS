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
