#-----------------------Construccion de las funciones---------------------------------------

#Esta distribución solo tiene un parametro de escala (σ).
#La distribución es continua.
#La funcion de enlace para el unico parametro que tenemos es "Log".



## The probability density function

dUMB <- function(x, mu = 1, log = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  
  sapply(x, function(xi) {
    if (xi <= 0 || xi >= 1) return(if (log) -Inf else 0)
    
    log_term <- log(1 / xi)
    log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(mu) - log(xi) - (log_term^2) / (2 * mu^2)
    
    if (log) return(log_density)
    else return(exp(log_density))
  })
}


## The cumulative distribution function

pUMB <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  
  sapply(q, function(qi) {
    if (qi <= 0 || qi >= 1) return(if (log.p) -Inf else 0)
    
    log_term <- log(1 / qi)
    erf_arg <- log_term / (sqrt(2) * mu)
    erf_part <- erf(erf_arg)
    
    term1 <- (sqrt(log_term^2) * erf_part) / log_term
    term2 <- (sqrt(2 / pi) * log_term * exp(-log_term^2 / (2 * mu^2))) / mu
    
    cdf <- 1 - term1 + term2
    if (!lower.tail) cdf <- 1 - cdf
    if (log.p) cdf <- log(cdf)
    
    return(cdf)
  })
}


## The quantile function

qUMB <- function(p, mu = 1) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  
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


## The random function

rUMB <- function(n = 1, mu = 1) {
  if (mu <= 0) stop("parameter mu must be positive!")
  u <- runif(n)
  qUMB(u, mu = mu)
}



#-------------------------Verificacion----------------------------------------

integrate(function(x) dUMB(x, mu = 1), lower = 0, upper = 1)


integrate(dUMB, lower=0, upper=1, mu= 3) #Esta es la forma de verificar que habia en la guia


sapply(c(1, 2, 3, 5, 10), function(th) {
  integrate(function(x) dUMB(x, mu = 1 ), lower = 0, upper = 1)$value
})


#-----------------------Ejemplos de dGEG---------------------------------------

dUMB(x, mu = 1, log = FALSE)

pUMB(q = 0.3 , mu = 1, lower.tail = TRUE, log.p = FALSE)

qUMB(p= 0.5, mu = 1)

rUMB(n = 3, mu = 0.3)


#-----------------------Graficas---------------------------------------


## The probability density function
curve(dUMB(x, mu= 0.10),
      from=0, to=1, col="blue", las=1, ylab="Pdf of the UMB",
      ylim = c(0,15), add=TRUE, )

curve(dUMB(x, mu= 0.25),
      from=0, to=1, col="red", las=1, 
      ylim = c(0,15), add=TRUE)

curve(dUMB(x, mu= 0.50), las=1, 
      ylim = c(0,15),
      add=TRUE, col="green")

curve(dUMB(x, mu= 1), las=1, 
      from=0, to=1, col="violetred1", las=1, 
      ylim = c(0,15), add=TRUE)

curve(dUMB(x, mu= 2), ylim = c(0,15),
      add=TRUE, col="black")

legend("topright", col=c("blue3", "red","green", "violetred1", "black"), lty=1, bty="n",
       legend=c("mu=0.10",
                "mu=0.25",
                "mu=0.50",
                "mu=1.0",
                "mu=2.0"))



## The cumulative distribution function
curve(pUMB(x, mu=0.10), from=0, to=1,
      col="black", las=1, ylab="Cdf of the UMB")

curve(pUMB(x, mu=0.25), from=0, to=1,
      add=TRUE, col="green", las=1)

curve(pUMB(x, mu=0.50), from=0, to=1,
      add=TRUE, col="violetred1", las=1)

curve(pUMB(x, mu=0.70), from=0, to=1,
      add=TRUE, col="red", las=1) 

curve(pUMB(x, mu=0.90), from=0, to=1,
      add=TRUE, col="blue3", las=1)


legend("topleft", col=c("black","green","violetred1", "red", "blue3"), lty=1, bty="n",
       legend=c("mu=0.10",
                "mu=0.25",
                "mu=0.50",
                "mu=0.70",
                "mu=0.9"))


## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qUMB(p, mu=0.25), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pUMB(x, mu=0.25), from=0, add=TRUE, col="red")



## The random function
hist(rUMB(n=1000,mu=0.25), freq=FALSE, xlab="x",
     ylim=c(0, 4.5), las=1, main="")
curve(dUMB(x, mu=0.25), add=TRUE, col="red")


#------------------------------Derivadas--------------------------------

# Libreria base
library(gamlss)

# Mi funcion de densidad dUMB
dUMB <- function(x, mu = 1, log = FALSE) {
  if (any(mu <= 0)) stop("parameter mu must be positive!")
  
  sapply(x, function(xi) {
    if (xi <= 0 || xi >= 1) return(if (log) -Inf else 0)
    
    log_term <- log(1 / xi)
    log_density <- 0.5 * log(2 / pi) + 2 * log(log_term) - 3 * log(mu) - log(xi) - (log_term^2) / (2 * mu^2)
    
    if (log) return(log_density)
    else return(exp(log_density))
  })
}

# Derivada  respecto a mu (dldd)
dldd_manual <- function(x, mu = 1) {
  L <- log(1 / x)
  return(-3 / mu + (L^2) / mu^3)
}



# Derivada computacional respecto a mu (dldd)
dldd_compu <- function(x, mu) {
  dm <- gamlss::numeric.deriv(expr = dUMB(x, mu, log = TRUE),
                              theta = "mu",
                              delta = 1e-04)
  dldd <- as.vector(attr(dm, "gradient"))
  return(dldd)
}

# Prueba con valores válidos de la distribución UMB (x > 0)
x     <- c(0.1, 0.2, 0.5, 0.9)  # Valores > 0
mu    <- 1.5


manual <- dldd_manual(x=x, mu=mu)
compu  <- dldd_compu(x=x, mu=mu)

# Comparación
comparison <- data.frame(x, manual, compu, difference = manual - compu)
print(comparison)



#-----------------------------------------------------------



