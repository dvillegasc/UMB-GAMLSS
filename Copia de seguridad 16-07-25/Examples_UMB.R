UMB(
  mu.link = "log"
)


# Example 1 - without covariates ------------------------------------------

n <- 100

# The true parameters are:
true_mu <- 0.30

true_theta <- c(true_mu)

# Graphing the pdf
curve(dUMB(x, mu=true_mu),
      ylab="Density", xlab="X", las=1,
      from=0, to=1, lwd=3, col="tomato")


# Simulating a random sample
y <- rUMB(n=n, mu=true_mu)

# Estimating paramaters
library(gamlss)
mod <- gamlss(y ~ 1, family=UMB,
              control=gamlss.control(n.cyc=100, trace=TRUE))

# Vector with the estimated results
res <- c(mu_hat=coef(mod, what="mu"))

# Comparing true vector and estimated vector
round(cbind(true_theta, with_UMB=res), digits=2)

# Histogram, estimated density and true density
MASS::truehist(y, ylab="Density", col="gray", las=1)

curve(dUMB(x, mu=res[1]),
      add=TRUE, col="blue2", lwd=2)

curve(dUMB(x, mu=true_theta[1]),
      add=TRUE, col="green4", lwd=2)

legend("topright", lwd=2, bty="n",
       legend=c("with UMB", "true density"),
       col=c("blue2", "green4"))



# Example 2 - with covariates ---------------------------------------------
n <- 500

# The true parameters are:
b0_mu <-  10
b1_mu <-  20

# The true theta vector
true_theta <- c(b0_mu, b1_mu)

# Simulating covariates
x1 <- runif(n, min=0, max=1)

# Simulating a random sample
y <- rUMB(n=n,
          mu    =    exp( b0_mu    + b1_mu    * x1))

# The dataframe
datos <- data.frame(y=y, x1=x1)

# Estimating paramaters
# Using gamlss with our proposal
mod <- gamlss(y ~ x1,
              family=UMB,
              data= datos,
              control=gamlss.control(n.cyc=10000, trace=TRUE))

# To obtain the estimated parameters
param <- unlist(coefAll(mod))

# Comparing true vector and estimated vector
res <- cbind(true_theta, with_gamlss= param)
round(res, digits=2)

