# Cargar librerías necesarias
library(gamlss)
library(MASS) # para truehist

# Definir la familia (si no está cargada en otro archivo)
UMB(mu.link = "log")

# ---------------------
# Ejemplo 1: sin covariables
# ---------------------

n <- 100
true_mu <- 0.3
true_theta <- c(true_mu)

# Graficar la densidad verdadera
curve(dUMB(x, mu = true_mu),
      ylab = "Density", xlab = "X", las = 1,
      from = 0.001, to = 0.999, lwd = 3, col = "tomato")

# Simular muestra aleatoria
y <- rUMB(n = n, mu = true_mu)

# Ajustar el modelo sin covariables
mod <- gamlss(y ~ 1, family = UMB,
              control = gamlss.control(n.cyc = 1000, trace = TRUE))

# Estimar parámetros
res <- c(mu_hat = coef(mod, what = "mu"))

# Comparación
round(cbind(true_theta, with_UMB = res), digits = 2)

# Histograma y densidades
truehist(y, ylab = "Density", col = "gray", las = 1)
curve(dUMB(x, mu = res[1]), add = TRUE, col = "blue2", lwd = 2)
curve(dUMB(x, mu = true_theta[1]), add = TRUE, col = "green4", lwd = 2)
legend("topright", lwd = 2, bty = "n",
       legend = c("with UMB", "true density"),
       col = c("blue2", "green4"))

# ---------------------
# Ejemplo 2: con covariables
# ---------------------

n <- 500
b0_mu <- 0.1
b1_mu <- 0.2
true_theta <- c(b0_mu, b1_mu)

# Covariable
x1 <- runif(n, min = 0, max = 1)

# Simular muestra
mu_vals <- b0_mu + b1_mu * x1
y <- rUMB(n = n, mu = mu_vals)

# Data frame
datos <- data.frame(y = y, x1 = x1)

# Ajuste del modelo con covariable
mod <- gamlss(y ~ x1,
              family = UMB,
              data = datos,
              control = gamlss.control(n.cyc = 10000, trace = TRUE))

# Obtener parámetros estimados
param <- coef(mod)

# Comparar
res <- cbind(true_theta, with_gamlss = param)
round(res, digits = 2)

