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
#' \eqn{f(x | \mu, \mu, \nu, \tau) = \frac{\tau}{\nu} \exp(w) \Phi \left( z - \frac{\mu}{\nu} \right) \left[ \Phi(z) - \exp(w)  \Phi \left( z - \frac{\mu}{\nu} \right) \right]^{\tau-1}}
#'
#' for \eqn{-\infty < x < \infty}. With \eqn{w=\frac{\mu-x}{\nu} + \frac{\mu^2}{2\nu^2}} and \eqn{z=\frac{x-\mu}{\mu}}
#' and \eqn{\Phi} is the cumulative function for the standard normal distribution.
#'
#' @example examples/examples_UMB.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @importFrom stats dnorm pnorm
#' @export
UMB <- function (mu.link="log") {
  
  dstats <- checklink("mu.link", "Unit Maxwell-Boltzmann",
                      substitute(mu.link), c("log", "logit", "probit", "cloglog", "own"))
  
  structure(list(family=c("UMB", "Unit Maxwell-Boltzmann"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Continuous",
                 
                 mu.link = as.character(substitute(mu.link)),
                 
                 mu.linkfun = dstats$linkfun,
                 
                 mu.linkinv = dstats$linkinv,
                 
                 mu.dr = dstats$mu.eta,
                 
                 # Primeras derivadas
                 
                 dldd = function(x, mu = 1, mu = 1, nu = 1, tau = 1) {
                   L <- log(1 / x)
                   dldd <- (-3 / mu + (L^2) / mu^3)
                   dldd
                 
                 # Segundas derivadas
                 
                 
                 d2ldmdd = function(y, mu, mu, nu, tau) {
                   L <- log(1 / x)
                   dldd <- (-3 / mu + (L^2) / mu^3)
                   
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldmdv = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- 1/nu
                   p2 <- - dnorm(z - mu/nu) / (mu * pnorm(z - mu/nu))
                   k1 <-  pnorm(z - mu/nu) / nu
                   k2 <- -dnorm(z - mu/nu) / mu
                   p3 <- (tau-1) * (-dnorm(z)/mu - exp(w) * (k1 + k2))
                   p3 <- p3 / (pnorm(z) - exp(w) * pnorm(z - mu/nu))
                   dldm <- p1 + p2 + p3
                   
                   p1 <- (-1/nu)
                   p2 <- (y-mu)/nu^2 - mu^2/nu^3
                   p3 <- dnorm(q) * (mu/nu^2) / pnorm(q)
                   k1 <- pnorm(q) * p2
                   k2 <- dnorm(q) * mu/nu^2
                   p4 <- (1-tau) * exp(w) * (k1 + k2)
                   p4 <- p4 / (pnorm(z) - exp(w) * pnorm(q))
                   dldv <- p1 + p2 + p3 + p4
                   
                   d2ldmdv <- - dldm * dldv
                   d2ldmdv
                 },
                 
                 d2ldmdt = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- 1/nu
                   p2 <- - dnorm(z - mu/nu) / (mu * pnorm(z - mu/nu))
                   k1 <-  pnorm(z - mu/nu) / nu
                   k2 <- -dnorm(z - mu/nu) / mu
                   p3 <- (tau-1) * (-dnorm(z)/mu - exp(w) * (k1 + k2))
                   p3 <- p3 / (pnorm(z) - exp(w) * pnorm(z - mu/nu))
                   dldm <- p1 + p2 + p3
                   
                   dldt <- 1/tau + log(pnorm(z) - exp(w) * pnorm(q))
                   
                   d2ldmdt <- - dldm * dldt
                   d2ldmdt
                 },
                 
                 d2ldd2  = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- mu/nu^2
                   p2 <- dnorm(q) * ((mu-y)/mu^2 - 1/nu) / pnorm(q)
                   k1 <- pnorm(q) * (mu/nu^2)
                   k2 <- dnorm(q) * ((mu-y)/mu^2 - 1/nu)
                   p3 <- (tau-1) * (dnorm(z) * (mu-y)/mu^2 - exp(w) * (k1 + k2))
                   p3 <- p3 / (pnorm(z) - exp(w) * pnorm(q))
                   dldd <- p1 + p2 + p3
                   
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },
                 
                 d2ldddv = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- (-1/nu)
                   p2 <- (y-mu)/nu^2 - mu^2/nu^3
                   p3 <- dnorm(q) * (mu/nu^2) / pnorm(q)
                   k1 <- pnorm(q) * p2
                   k2 <- dnorm(q) * mu/nu^2
                   p4 <- (1-tau) * exp(w) * (k1 + k2)
                   p4 <- p4 / (pnorm(z) - exp(w) * pnorm(q))
                   dldv <- p1 + p2 + p3 + p4
                   
                   p1 <- mu/nu^2
                   p2 <- dnorm(q) * ((mu-y)/mu^2 - 1/nu) / pnorm(q)
                   k1 <- pnorm(q) * (mu/nu^2)
                   k2 <- dnorm(q) * ((mu-y)/mu^2 - 1/nu)
                   p3 <- (tau-1) * (dnorm(z) * (mu-y)/mu^2 - exp(w) * (k1 + k2))
                   p3 <- p3 / (pnorm(z) - exp(w) * pnorm(q))
                   dldd <- p1 + p2 + p3
                   
                   d2ldddv <- - dldd * dldv
                   d2ldddv
                 },
                 
                 d2ldddt = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- mu/nu^2
                   p2 <- dnorm(q) * ((mu-y)/mu^2 - 1/nu) / pnorm(q)
                   k1 <- pnorm(q) * (mu/nu^2)
                   k2 <- dnorm(q) * ((mu-y)/mu^2 - 1/nu)
                   p3 <- (tau-1) * (dnorm(z) * (mu-y)/mu^2 - exp(w) * (k1 + k2))
                   p3 <- p3 / (pnorm(z) - exp(w) * pnorm(q))
                   dldd <- p1 + p2 + p3
                   
                   dldt <- 1/tau + log(pnorm(z) - exp(w) * pnorm(q))
                   
                   d2ldddt <- - dldd * dldt
                   d2ldddt
                 },
                 
                 d2ldv2 = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- (-1/nu)
                   p2 <- (y-mu)/nu^2 - mu^2/nu^3
                   p3 <- dnorm(q) * (mu/nu^2) / pnorm(q)
                   k1 <- pnorm(q) * p2
                   k2 <- dnorm(q) * mu/nu^2
                   p4 <- (1-tau) * exp(w) * (k1 + k2)
                   p4 <- p4 / (pnorm(z) - exp(w) * pnorm(q))
                   dldv <- p1 + p2 + p3 + p4
                   
                   d2ldv2 <- - dldv * dldv
                   d2ldv2
                 },
                 
                 d2ldvdt = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   p1 <- (-1/nu)
                   p2 <- (y-mu)/nu^2 - mu^2/nu^3
                   p3 <- dnorm(q) * (mu/nu^2) / pnorm(q)
                   k1 <- pnorm(q) * p2
                   k2 <- dnorm(q) * mu/nu^2
                   p4 <- (1-tau) * exp(w) * (k1 + k2)
                   p4 <- p4 / (pnorm(z) - exp(w) * pnorm(q))
                   dldv <- p1 + p2 + p3 + p4
                   
                   dldt <- 1/tau + log(pnorm(z) - exp(w) * pnorm(q))
                   
                   d2ldvdt <- - dldv * dldt
                   d2ldvdt
                 },
                 
                 d2ldt2 = function(y, mu, mu, nu, tau) {
                   z <- (y-mu)/mu
                   w <- (mu-y)/nu + mu^2/(2*nu^2)
                   q <- z - mu/nu
                   
                   dldt <- 1/tau + log(pnorm(z) - exp(w) * pnorm(q))
                   
                   d2ldt2 <- - dldt * dldt
                   d2ldt2
                 },
                 
                 G.dev.incr = function(y, mu, mu, nu, tau, ...) -2*dUMB(y, mu, mu, nu, tau, log=TRUE),
                 rqres      = expression(rqres(pfun="pUMB", type="Continuous", y=y, mu=mu, mu=mu, nu=nu, tau=tau)),
                 
                 mu.initial    = expression(mu    <- rep(estim_mu_mu_nu_tau_UMB(y)[1], length(y)) ),
                 mu.initial = expression(mu <- rep(estim_mu_mu_nu_tau_UMB(y)[2], length(y)) ),
                 nu.initial    = expression(nu    <- rep(estim_mu_mu_nu_tau_UMB(y)[3], length(y)) ),
                 tau.initial   = expression(tau   <- rep(estim_mu_mu_nu_tau_UMB(y)[4], length(y)) ),
                 
                 mu.valid    = function(mu)    TRUE,
                 mu.valid = function(mu) all(mu > 0),
                 nu.valid    = function(nu)    all(nu > 0),
                 tau.valid   = function(tau)   all(tau > 0),
                 
                 y.valid = function(y) TRUE
  ),
  class=c("gamlss.family", "family"))
}
#'
#' estim_mu_mu_nu_tau_UMB
#'
#' This function generates initial values for UMB distribution.
#'
#' @param y vector with the random sample
#' @examples
#' y <- rUMB(n=100, mu=1, mu=1, nu=1, tau=1)
#' estim_mu_mu_nu_tau_UMB(y=y)
#' @importFrom stats optim
#' @export
estim_mu_mu_nu_tau_UMB <- function(y) {
  mod <- optim(par=c(0, 0, 0, 0),
               fn=logLik_UMB,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    =     mod$par[1],
           mu_hat = exp(mod$par[2]),
           nu_hat    = exp(mod$par[3]),
           tau_hat   = exp(mod$par[4]))
  #res <- c(0, 1, 1, 1) # esto se lo puse para que no tenga en cuenta optim
  names(res) <- c("mu_hat", "mu_hat", "nu_hat", "tau_hat")
  return(res)
}
#'
#' logLik_UMB
#'
#' This is an auxiliar function to obtain the logLik for UMB.
#'
#' @param logparam vector with the values for mu, mu, nu and tau
#' @param x vector with the data
#' @examples
#' y <- rUMB(n=100, mu=1, mu=1, nu=1, tau=1)
#' logLik_UMB(logparam=c(0, 0, 0, 0), x=y)
#' @importFrom stats optim
#' @export
logLik_UMB <- function(logparam=c(0, 0, 0, 0), x){
  return(sum(dUMB(x,
                  mu    = logparam[1],
                  mu = exp(logparam[2]),
                  nu    = exp(logparam[3]),
                  tau   = exp(logparam[4]),
                  log=TRUE)))
}