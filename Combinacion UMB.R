dUMB <- function(y, theta=1, log=FALSE){
  if (any(theta <= 0))
    stop("parameter sigma has to be positive!")
  
  
  
  log_term <- log(1 / y)
  pdf <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * theta^2))) / ((theta^3)* y)
  #log_density <- log(pdensity)
  
  if(log)
    result <- pdf
  else
    result <- exp(pdf)
  # Too small to compute, zero it
  result[is.nan(result)] <- 0
  return(result)
}
#' @importFrom stats pnorm
#' @export
#' @rdname dUMB
 

pUMB <- function(x, theta=1, sigma=1, nu=1, tau=1, lower.tail=TRUE, log.p=FALSE){
  if (any(theta <= 0))
    stop("parameter sigma has to be positive!")
  
  log_term <- log(1 / y)
  sqrt_log2 <- sqrt(log_term^2)
  erf_arg <- log_term / (sqrt(2) * theta)
  erf_part <- erf(erf_arg)
  
  term1 <- (sqrt_log2 * erf_part) / log_term
  term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * theta^2))) / theta
  
  cdf <- 1 - term1 + term2
  
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- cdf
  else cdf <- exp(cdf)
  cdf
}
#' @importFrom gamlss.dist qexGAUS
#' @export
#' @rdname dUMB
 

qUMB <- function(p, theta, lower.tail = TRUE, log.p = FALSE, tol = 1e-9) {
  if (any(theta <= 0)) stop("theta thetast be > 0")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0 | p > 1)) stop("p thetast be in [0, 1]")
  
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


#' @importFrom stats runif
#' @importFrom gamlss.dist qexGAUS
#' @export
#' @rdname dUMB
rUMB <- function(n, theta) {
  if (theta <= 0) stop("theta thetast be > 0")
  u <- runif(n)
  qUMB(u, theta)
}




integrate(dUMB, lower=0, upper=1, theta=1) 



curve(dUMB(y, mu=7),
      from=0, to=10, col="red", las=1, ylab="f(x)")

curve(dUMB(y, mu=7),
      add=TRUE, col="blue3")

curve(dUMB(y, mu=7),
      add=TRUE, col="green4")

legend("topleft", col=c("red", "blue3", "green4"), lty=1, bty="n",
       legend=c("mu=7",
                "mu=7",
                "mu=7"))
