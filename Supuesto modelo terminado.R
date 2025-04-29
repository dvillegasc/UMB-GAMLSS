#' @export
dUMB <- function(x, mu=1, sigma=1, nu=1, tau=1, log=FALSE){
  if (any(mu <= 0)) stop("parameter mu must be > 0")
  
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  
  valid <- (x > 0 & x < 1)
  
  if (any(valid)) {
    log_term <- log(1 / x[valid])
    log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(mu) -
      log(x[valid]) - (log_term^2) / (2 * mu^2)
    
    if (log) {
      out[valid] <- log_density
    } else {
      out[valid] <- exp(log_density)
    }
  }
  
  out[!valid] <- if (log) -Inf else 0
  return(out)
}

#' @export
pUMB <- function(q, mu=1, sigma=1, nu=1, tau=1, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0)) stop("parameter mu must be > 0")
  
  q <- as.numeric(q)
  out <- rep(NA_real_, length(q))
  
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  
  valid <- (q > 0 & q < 1)
  
  if (any(valid)) {
    log_term <- log(1 / q[valid])
    sqrt_log2 <- sqrt(log_term^2)
    erf_arg <- log_term / (sqrt(2) * mu)
    erf_part <- erf(erf_arg)
    
    term1 <- (sqrt_log2 * erf_part) / log_term
    term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * mu^2))) / mu
    
    cdf <- 1 - term1 + term2
    
    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)
    
    out[valid] <- cdf
  }
  
  out[!valid] <- if (log.p) -Inf else 0
  
  return(out)
}

#' @export
qUMB <- function(p, mu=1, sigma=1, nu=1, tau=1){
  if (any(mu <= 0)) stop("parameter mu must be > 0")
  
  if (any(p < 0 | p > 1)) stop("p must be in [0,1]")
  
  quantile_fn <- function(prob) {
    sapply(prob, function(pp) {
      if (pp == 0) return(0)
      if (pp == 1) return(1)
      uniroot(function(y) pUMB(y, mu = mu) - pp,
              lower = .Machine$double.eps, upper = 1 - .Machine$double.eps,
              tol = 1e-9)$root
    })
  }
  
  return(quantile_fn(p))
}

#' @export
rUMB <- function(n=1, mu=1, sigma=1, nu=1, tau=1){
  if (any(mu <= 0)) stop("parameter mu must be > 0")
  
  u <- runif(n)
  return(qUMB(u, mu = mu))
}


curve(function(x) dUMB(x, mu=3), from=0, to=1, col="blue", lwd=2, ylab="Density")
curve(dbeta(x, 2, 5), add=TRUE, col="red", lty=2)
legend("topright", legend=c("UMB (mu=3)", "Beta(2,5)"), col=c("blue", "red"), lty=1:2, bty="n")










UMB <- function(mu.link="log") {
  
  mstats <- checklink("mu.link", "UMB", substitute(mu.link), c("log", "identity"))
  
  structure(
    list(
      family = c("UMB"),
      parameters = list(mu=TRUE),
      nopar = 1,
      type = "Continuous",
      
      mu.link = as.character(substitute(mu.link)),
      mu.linkfun = mstats$linkfun,
      mu.linkinv = mstats$linkinv,
      mu.dr = mstats$mu.eta,
      
      dldm = function(y, mu) {
        log_term <- log(1 / y)
        (log_term^2 - mu^2) / (mu^3)
      },
      
      d2ldm2 = function(y, mu) {
        log_term <- log(1 / y)
        (3 * mu^2 - 2 * log_term^2) / (mu^4)
      },
      
      G.dev.incr = function(y, mu, ...) -2 * log(dUMB(y, mu=mu)),
      
      rqres = expression(
        rqres(pfun="pUMB", type="Continuous", y=y, mu=mu)
      ),
      
      mu.initial = expression(mu <- rep(1, length(y))),
      
      mu.valid = function(mu) all(mu > 0),
      
      y.valid = function(y) all(y > 0 & y < 1),
      
      mean = function(mu, ...) NA, # No disponible
      variance = function(mu, ...) NA # No disponible
    ),
    class = c("gamlss.family", "family")
  )
}


library(gamlss)

set.seed(123)
y <- rUMB(200, mu=3)

fit <- gamlss(y ~ 1, family=UMB)
summary(fit)

