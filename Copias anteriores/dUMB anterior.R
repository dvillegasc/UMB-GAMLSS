dUMB <- function(y, theta, log = FALSE) {
  if (any(theta <= 0)) stop("theta must be > 0")
  if (any(y <= 0 | y >= 1)) return(if (log) -Inf else 0)
  
  log_term <- log(1 / y)
  log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(theta) - log(y) - (log_term^2) / (2 * theta^2)
  
  
  #pdensity <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * theta^2))) / ((theta^3)* y)
  #log_density <- log(pdensity)
  
  
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}


pUMB <- function(y, theta, lower.tail = TRUE, log.p = FALSE) {
  if (any(theta <= 0)) stop("theta must be > 0")
  if (any(y <= 0 | y >= 1)) return(if (log.p) -Inf else 0)
  
  log_term <- log(1 / y)
  sqrt_log2 <- sqrt(log_term^2)
  erf_arg <- log_term / (sqrt(2) * theta)
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  erf_part <- erf(erf_arg)
  
  term1 <- (sqrt_log2 * erf_part) / log_term
  term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * theta^2))) / theta
  
  cdf <- 1 - term1 + term2
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}


qUMB <- function(p, theta, lower.tail = TRUE, log.p = FALSE, tol = 1e-9) {
  if (any(theta <= 0)) stop("theta must be > 0")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  
  quantile_fn <- function(prob) {
    sapply(prob, function(pp) {
      if (pp == 0) return(0)
      if (pp == 1) return(1)
      uniroot(function(y) pUMB(y, theta) - pp,
              lower = .Machine$double.eps, upper = 1 - .Machine$double.eps,
              tol = tol)$root
    })
  }
  
  return(quantile_fn(p))
}


rUMB <- function(n, theta) {
  if (theta <= 0) stop("theta must be > 0")
  u <- runif(n)
  qUMB(u, theta)
}


integrate(function(y) dUMB(y, theta=1), lower=0, upper=1)


integrate(function(y) dUMB(y, theta = 3), lower = 0, upper = 1)




sapply(c(1, 2, 3, 5, 10), function(th) {
  integrate(function(y) dUMB(y, theta = th), lower = 0, upper = 1)$value
})





curve(function(x) dUMB(x, theta = 3), from = 0, to = 1, col = "blue", lwd = 2, ylab = "Density")
curve(dbeta(x, 2, 5), add = TRUE, col = "red", lty = 2)
legend("topright", legend = c("UMB (theta=3)", "Beta(2,5)"), col = c("blue", "red"), lty = 1:2, bty = "n")

dUMB(y, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE)

pUMB(y, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)

qUMB(p, mu = 0, sigma = 1, nu = 1, tau = 1)

rUMB(n = 1, mu = 0, sigma = 1, nu = 1, tau = 1) 












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
y <- rUMB(200, theta=1)

fit <- gamlss(y ~ 1, family=UMB)
summary(fit)


get.dist.family("")
