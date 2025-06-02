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


estim_mu_UMB <- function(y) {
  mu_hat <- sqrt(sum(log(y)^2) / (3 * length(y)))
  res <- c(mu_hat = mu_hat)
  names(res) <- "mu_hat"
  return(res)
}


logLik_UMB <- function(param=c(0), x){
  return(sum(dUMB(x, mu = exp(param[1]), log=TRUE)))
}

