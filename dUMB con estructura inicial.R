#-----------------------Construccion de las funciones---------------------------------------

#Esta distribución solo tiene un parametro de escala (σ).
#La distribución es continua.
#La funcion de enlace para el unico parametro que tenemos es "Log".



## The probability density function

dUMB <- function(x, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE) {
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  
  sapply(x, function(xi) {
    if (xi <= 0 || xi >= 1) return(if (log) -Inf else 0)
    
    log_term <- log(1 / xi)
    log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(sigma) - log(xi) - (log_term^2) / (2 * sigma^2)
    
    if (log) return(log_density)
    else return(exp(log_density))
  })
}


## The cumulative distribution function

pUMB <- function(q, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  
  sapply(q, function(qi) {
    if (qi <= 0 || qi >= 1) return(if (log.p) -Inf else 0)
    
    log_term <- log(1 / qi)
    erf_arg <- log_term / (sqrt(2) * sigma)
    erf_part <- erf(erf_arg)
    
    term1 <- (sqrt(log_term^2) * erf_part) / log_term
    term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * sigma^2))) / sigma
    
    cdf <- 1 - term1 + term2
    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)
    
    return(cdf)
  })
}


## The quantile function

qUMB <- function(p, mu = 0, sigma = 1, nu = 1, tau = 1) {
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  
  quantile_fn <- function(prob) {
    sapply(prob, function(pp) {
      if (pp == 0) return(0)
      if (pp == 1) return(1)
      uniroot(function(y) pUMB(y, sigma = sigma) - pp,
              lower = .Machine$double.eps, upper = 1 - .Machine$double.eps,
              tol = 1e-9)$root
    })
  }
  
  return(quantile_fn(p))
}


## The random function

rUMB <- function(n = 1, mu = 0, sigma = 1, nu = 1, tau = 1) {
  if (sigma <= 0) stop("parameter sigma must be positive!")
  u <- runif(n)
  qUMB(u, sigma = sigma)
}



#-------------------------Verificacion----------------------------------------

integrate(function(x) dUMB(x, sigma = 1), lower = 0, upper = 1)


integrate(dUMB, lower=0, upper=1, mu=-6, sigma=3, nu=4) #Esta es la forma de verificar que habia en la guia


sapply(c(1, 2, 3, 5, 10), function(th) {
  integrate(function(x) dUMB(x, sigma = 1 ), lower = 0, upper = 1)$value
})


#-----------------------Ejemplos de dGEG---------------------------------------

dUMB(x, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE)

pUMB(q, mu = 0, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)

qUMB(p, mu = 0, sigma = 1, nu = 1, tau = 1)

rUMB(n = 1, mu = 0, sigma = 1, nu = 1, tau = 1)


#-----------------------Graficas---------------------------------------


## The probability density function
curve(dUMB(x, mu=7, sigma= 0.10, nu=0.05, tau=1),
      from=0, to=1, col="blue", las=1, ylab="f(x)")

curve(dUMB(x, mu=7, sigma= 0.25, nu=0.05, tau=1),
      from=0, to=1, col="red", las=1, ylab="f(x)")

curve(dUMB(x, mu=7, sigma= 0.50, nu=0.05, tau=0.5),
      add=TRUE, col="green")

curve(dUMB(x, mu=7, sigma= 1, nu=0.05, tau=1),
      from=0, to=1, col="pink", las=1, ylab="f(x)")

curve(dUMB(x, mu=7, sigma= 2, nu=0.05, tau=0.1),
      add=TRUE, col="black")

legend("topright", col=c("red", "blue3", "green4"), lty=1, bty="n",
       legend=c("mu=7, sigma=0.25, nu=0.05, tau=1",
                "mu=7, sigma=0.25, nu=0.05, tau=0.5",
                "mu=7, sigma=0.25, nu=0.05, tau=0.1"))



## The cumulative distribution function
curve(pUMB(x, mu=7, sigma=0.10, nu=0.05, tau=1), from=0, to=1,
      col="black", las=1, ylab="F(x)")

curve(pUMB(x, mu=7, sigma=0.25, nu=0.05, tau=0.5), from=0, to=1,
      add=TRUE, col="green", las=1)

curve(pUMB(x, mu=7, sigma=0.50, nu=0.05, tau=0.5), from=0, to=1,
      add=TRUE, col="pink3", las=1)

curve(pUMB(x, mu=7, sigma=0.70, nu=0.05, tau=0.5), from=0, to=1,
      add=TRUE, col="red", las=1) 

curve(pUMB(x, mu=7, sigma=0.90, nu=0.05, tau=0.5), from=0, to=1,
      add=TRUE, col="blue3", las=1)




## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qUMB(p, mu=7, sigma=0.25, nu=0.05, tau=1), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pUMB(x, mu=7, sigma=0.25, nu=0.05, tau=1), from=0, add=TRUE, col="red")



## The random function
hist(rUMB(n=1000, mu=7, sigma=0.25, nu=0.05, tau=1), freq=FALSE, xlab="x",
     ylim=c(0, 6), las=1, main="")
curve(dUMB(x, mu=7, sigma=0.25, nu=0.05, tau=1), add=TRUE, col="red")


#------------------------------Derivadas--------------------------------

# Libreria base
library(gamlss)

# Mi funcion de densidad dUMB
dUMB <- function(x, mu = 0, sigma = 1, nu = 1, tau = 1, log = FALSE) {
  if (any(sigma <= 0)) stop("parameter sigma must be positive!")
  
  sapply(x, function(xi) {
    if (xi <= 0 || xi >= 1) return(if (log) -Inf else 0)
    
    log_term <- log(1 / xi)
    log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(sigma) - log(xi) - (log_term^2) / (2 * sigma^2)
    
    if (log) return(log_density)
    else return(exp(log_density))
  })
}

# Derivada  respecto a sigma (dldd)
dldd_manual <- function(x, mu = 1, sigma = 1, nu = 1, tau = 1) {
  L <- log(1 / x)
  return(-3 / sigma + (L^2) / sigma^3)
}



# Derivada computacional respecto a sigma (dldd)
dldd_compu <- function(x, mu, sigma, nu, tau) {
  dm <- gamlss::numeric.deriv(expr = dUMB(x, mu, sigma, nu, tau, log = TRUE),
                              theta = "sigma",
                              delta = 1e-04)
  dldd <- as.vector(attr(dm, "gradient"))
  return(dldd)
}

# Prueba con valores válidos de la distribución UMB (x > 0)
x     <- c(0.1, 0.2, 0.5, 0.9)  # Valores > 0
mu    <- 1
sigma <- 1.5
nu    <- 2
tau   <- 1

manual <- dldd_manual(x=x, mu=mu, sigma=sigma, nu=nu, tau=tau)
compu  <- dldd_compu(x=x, mu=mu, sigma=sigma, nu=nu, tau=tau)

# Comparación
comparison <- data.frame(x, manual, compu, difference = manual - compu)
print(comparison)



#-----------------------------------------------------------



