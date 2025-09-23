library(gamlss)

# Vamos a definir el directorio de trabajo donde
# se almacenaran el script, simulaciones y graficos

#setwd("")

#-----------------------------------------------------------------
#------------------ EJEMPLO CON LA DISTRIBUCION ------------------
#----------------------- Unit Maxwell-Boltzmann ------------------
#---------------------------- UMB --------------------------------
#-----------------------------------------------------------------

## En la distribucion UMB
# mu > 0,         por tanto usaremos funcion de enlace log

# Creando las funciones de enlace inversas
# No vamos a crear log_inv porque esa funcion ya existe y 
# llama exp( )

# The parameters ----------------------------------------------------------
true_mu    <- 1.74 # Cumple que mu > 0

# Useful functions to the simulation study --------------------------------

# Funcion para obtener mu_hat para un valor fijo de n
simul_one <- function(size) {
  y <- rUMB(n=size, mu=true_mu)
  y[y <= 0] <- 1e-8
  y[y >= 1] <- 1 - 1e-8
  mod <- NULL
  mod <- gamlss(y~1, sigma.fo=~1, family='UMB',
                control=gamlss.control(n.cyc=2500, trace=FALSE))
  res <- c(      exp(coef(mod, what='mu')))
  res
}

# Super function to simulate and write the estimated parameters
simul <- function(n) {
  result <- t(replicate(n=nrep, expr=simul_one(size=n)))
  result <- cbind(result, n)
  write(x=t(result), file='simul_without_cov.txt', 
        ncol=2, append=TRUE)
}

# Code to generate the simulations given n --------------------------------

# Aqui se definen los valores de tamano muestral n
# Luego se define el numero de repeticiones
n <- seq(from=20, to=300, by=20)
nrep <- 1000

values <- expand.grid(n=n)
values
apply(values, 1, simul)


# Plots -------------------------------------------------------------------
dt <- read.table('simul_without_cov.txt', 
                 col.names=c('mu_hat', 'n'))

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Numero de observaciones por cada n
num <- dt %>% group_by(n) %>% count()
mean(num$nn)
min(num$nn)

# Para obtener la metricas
res <- dt %>% 
  drop_na() %>% 
  group_by(n) %>% 
  summarise(mean_mu=mean(mu_hat),
            bias_mu=true_mu - mean(mu_hat),
            mse_mu=mean((true_mu - mu_hat)^2))
# Mean -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mean_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(mu))) +
  geom_line(y=true_mu, col='red', lty='dashed')


mean1 <- grid.arrange(p1, nrow = 1)
mean1
ggsave(filename="mean1.pdf", 
       plot=mean1, 
       width=10, height=4)

# Bias -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=bias_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(Bias~hat(mu)))

bias1 <- grid.arrange(p1, nrow = 1)
bias1
ggsave(filename="bias1.pdf", 
       plot=bias1, 
       width=10, height=4)

# MSE -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mse_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE~hat(mu)))

mse1 <- grid.arrange(p1, nrow = 1)
mse1
ggsave(filename="mse1.pdf", 
       plot=mse1, 
       width=10, height=4)
