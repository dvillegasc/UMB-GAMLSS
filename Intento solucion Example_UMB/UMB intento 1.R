#' Unit Maxwell-Boltzmann family
#'
#' @description
#' The function \code{UMB()} defines the Unit Maxwell-Boltzmann distribution, a one parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu.
#'
#' @details
#' The Unit Maxwell-Boltzmann with parameters \code{mu} has density given by
#'
#' \deqn{f(y, \mu) = \frac{\sqrt{\frac{\pi}{2}} \, \ln^{2}\left(\frac{1}{y}\right) e^{-\frac{\ln^{2}\left(\frac{1}{y}\right)}{2\theta^{2}}}}{\theta^{3} y}
#'
#' for \eqn{0 < x < 1}.
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
                 
                 mu.initial = expression(mu <- rep(estim_mu_UMB(y)[1], length(y)) ),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 
                 y.valid = function(y) all(y > 0) && all(y < 1)
        
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
  mu_hat <- sqrt(sum(log(y)^2) / (3 * length(y)))
  res <- c(mu_hat = mu_hat)
  names(res) <- "mu_hat"
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
  return(sum(dUMB(x, mu = exp(param[1]), log=TRUE)))
}

