fun_rate <- function(b1, b2){
  expit(b1 * X[1,] + b2 * X[2,]);
}

lw.normalize <- function (lw){
  maxlw <- max(lw)
  lw <- lw - maxlw
  return(exp(lw)/sum(exp(lw)))
}

expit <- function(z){ 1/(1 + exp(-z)) }

logit <- function(z){log(z/(1-z))}

sis_get_alpha <- function (agent_state, model_config) {
  alpha <- agent_state * (1 - model_config$gamma) +
    (1 -  agent_state) * (model_config$lambda * (sum(agent_state)/model_config$N))
  return(as.vector(alpha))
}

update_config <- function(config, parameters){
  config$alpha0 <- fun_rate(parameters[1], parameters[2]);
  config$lambda <- fun_rate(parameters[3], parameters[4]);
  config$gamma <- fun_rate(parameters[5], parameters[6]);
  config$rho <- parameters[7];
  return(config);
}

update_config_fix_gam <- function(config, parameters){
  config$alpha0 <- fun_rate(parameters[1], parameters[2]);
  config$lambda <- fun_rate(parameters[3], parameters[4]);
  #config$gamma <- fun_rate(parameters[5], parameters[6]);
  config$rho <- parameters[7];
  return(config);
}

qkernel <- function(parameters, step_size){
  z <- step_size * rnorm(7);
  proposal <- rep(0, 7)
  proposal[1:6] <- parameters[1:6] + z[1:6];
  proposal[7] <- expit(logit(parameters[7]) + z[7]);
  return(proposal)
}


ljacobian <- function(parameters){
  + log(parameters[7]) + log(1 - parameters[7]);
}

lprior_fun <- function(parameters, Hypers){
  TruncatedNormal::dtmvnorm(c(parameters[1:6], logit(parameters[7])),
           mu = Hypers$mu, sigma = diag(Hypers$sigma), lb = Hypers$lower,
           ub= Hypers$upper, log = TRUE)
}

lprior_fix_gam <- function(parameters, Hypers){
  TruncatedNormal::dtmvnorm(c(parameters[1:4], logit(parameters[7])),
                     mu = Hypers$mu[c(1:4,7)], sigma = diag(Hypers$sigma[c(1:4,7)]), lb = Hypers$lower[c(1:4,7)],
                     ub = Hypers$upper[c(1:4,7)], log = TRUE)
}

lprior_PG <- function(parameters, Hypers){
  TruncatedNormal::dtmvnorm(parameters[1:6],
           mu = Hypers$mu[1:6], sigma = diag(Hypers$sigma[1:6]), lb = Hypers$lower[1:6],
           ub = Hypers$upper[1:6], log = TRUE)
}



