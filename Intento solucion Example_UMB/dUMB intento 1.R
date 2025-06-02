


dUMB <- function(x, mu = 1, log = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  if (any(x <= 0 | x >= 1)) stop("x must be in (0, 1)")
  
  log_term <- log(1 / x)
  res <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(mu) - log(x) - (log_term^2) / (2 * mu^2)
  
  if(log)
    return(res)
  else
    return(exp(res))
}


pUMB <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) {
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  
  # Begin aux fun
  auxfun <- function(q) {
    log_term <- log(1 / q)
    erf_arg <- log_term / (sqrt(2) * mu)
    erf_part <- erf(erf_arg)
    term1 <- abs(log_term) * erf_part / log_term
    term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * mu^2))) / mu
    res <- 1 - term1 + term2
    res
  }
  # End aux fun
  
  cdf <- ifelse(q<=0, 0,
                ifelse(q >=1, 1, auxfun(q)))
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  return(cdf)
}



qUMB <- function(p, mu, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0 | p > 1)) stop("p must be in the interval (0, 1)")
  if (any(mu < 0))        stop("mu must be positive")
  
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  
  aux_fun <- function(x, p, mu){
    res <- p - pUMB(x, mu)
    return(res)
  }
  res <- uniroot(aux_fun, interval=c(0, 1),
                 p=p, mu=mu)$root
  
  return(res)
}
qUMB <- Vectorize(qUMB)



rUMB <- function(n = 1, mu = 1) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  u <- runif(n)
  # Si mu es un solo valor, replicamos. Si es vector, usamos tal cual.
  if (length(mu) == 1) mu <- rep(mu, n)
  sapply(1:n, function(i) qUMB(u[i], mu[i]))
}


