#' Generalised exponential-Gaussian family
#'
#' @description
#' The function \code{UMB()} defines the Generalised exponential-Gaussian distribution, a four parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "identity" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu.
#' @param tau.link defines the tau.link, with "log" link as the default for the tau.
#'
#' @details
#' The Generalised exponential-Gaussian with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}
#' has density given by
#'
#' \eqn{f(x | \mu, \sigma, \nu, \tau) = \frac{\tau}{\nu} \exp(w) \Phi \left( z - \frac{\sigma}{\nu} \right) \left[ \Phi(z) - \exp(w)  \Phi \left( z - \frac{\sigma}{\nu} \right) \right]^{\tau-1}}
#'
#' for \eqn{-\infty < x < \infty}. With \eqn{w=\frac{\mu-x}{\nu} + \frac{\sigma^2}{2\nu^2}} and \eqn{z=\frac{x-\mu}{\sigma}}
#' and \eqn{\Phi} is the cumulative function for the standard normal distribution.
#'
#' @example examples/examples_UMB.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @importFrom stats dnorm pnorm
#' @export
UMB <- function (mu.link="log") {
  
    mstats <- checklink("mu.link", "Unit Maxwell-Boltzmann",
                        substitute(mu.link), c("log", "inverse", "own"))
    
    structure(list(family=c("UMB", "Unit Maxwell-Boltzmann"),
                   parameters=list(mu=TRUE),
                   nopar=1,
                   type="Continuous",
                 
                 mu.link    = as.character(substitute(mu.link)),

                 
                 mu.linkfun    = mstats$linkfun,

                 
                 mu.linkinv    = mstats$linkinv,

                 
                 mu.dr    = mstats$mu.eta,

                 
                 # Primeras derivadas
                 
                 dldm = function(y, mu) {
                   L <- log(1 / y)
                   dldd <- (-3 / mu + (L^2) / mu^3)
                   dldd
                 },
                 
                 
                 # Segundas derivadas
                 
                 d2ldm2 = function(y, mu) {
                   L <- log(1 / y)
                   d2ldm2 <- 3 / mu^2 - (3 * L^2) / mu^4
                   return(d2ldm2)
                 },
                 

                 
                 
                 G.dev.incr = function(y, mu, ...) -2*dUMB(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pUMB", type="Continuous", y=y, mu=mu)),
                 
                 mu.initial    = expression(mu    <- rep(estim_mu_UMB(y)[1], length(y)) ),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 
                 y.valid = function(y) TRUE
  ),
  class=c("gamlss.family", "family"))
}
#'
#' estim_mu_UMB
#'
#' This function generates initial values for UMB distribution.
#'
#' @param y vector with the random sample
#' @examples
#' y <- rUMB(n=100, mu=1)
#' estim_mu_UMB(y=y)
#' @importFrom stats optim
#' @export
estim_mu_UMB <- function(y) {
  theta_hat <- sqrt(sum(log(y)^2) / (3 * length(y)))
  res <- c(mu_hat = theta_hat)
  names(res) <- c("mu_hat")
  return(res)
}
#'
#' logLik_UMB
#'
#' This is an auxiliar function to obtain the logLik for UMB.
#'
#' @param param vector with the values for mu
#' @param x vector with the data
#' @examples
#' y <- rUMB(n=100, mu=1)
#' logLik_UMB(param=c(0), x=y)
#' @importFrom stats optim
#' @export
logLik_UMB <- function(param=c(0), x){
  return(sum(dUMB(x,
                  mu    = exp(logparam[1])),
                  log=TRUE))
}

