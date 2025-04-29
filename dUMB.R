dUMB <- function(y, theta, log = FALSE) {
  if (any(theta <= 0)) stop("theta must be > 0")
  if (any(y <= 0 | y >= 1)) return(if (log) -Inf else 0)
  
  log_term <- log(1 / y)
  log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(theta) - log(y) - (log_term^2) / (2 * theta^2)
  
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
  sqrt_ln2 <- sqrt(log_term^2)
  erf_arg <- log_term / (sqrt(2) * theta)
  erf_part <- erf(erf_arg)
  
  term1 <- sqrt_ln2 * erf_part / log_term
  term2 <- sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * theta^2)) / theta
  
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


integrate(dUMB, lower=0, upper=1, theta=1) 
