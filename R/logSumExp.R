#' Calculate log sum of exponentials
#'
#' @param mu is a vector of values in the log scale. They typically correspond to
#' log likelihood values for a Bayesian model.
#'
#' @return The log of the sum of exponentials
#' @details The log sum of exponentials is calculated as:
#' log-sum-exp_m=1^M mu_m = max(mu) + log sum_m^M exp(mu_m - max(mu)),
#' where mu_m are typically log likelihood values for a Bayesian model obtained
#' for iteration m of an MCMC algorithm with a total of M iterations.
#' This prevents underflow. See https://mc-stan.org/docs/stan-users-guide/
#' for more information
#'
#' @export
#'
#' @examples
#' log_sum_exp(log(5))
#'
logSumExp <- function(mu){

  maxmu <- max(mu)

  maxmu + log(sum(exp(mu - maxmu)))

}
