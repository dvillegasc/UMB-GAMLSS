library(testthat)
library(gamlss)  # si estás probando en un entorno GAMLSS
source("C:/Users/davil/Desktop/A.GAMLSS/UMB-GAMLSS/Distribucion UMB/UMB correcto.R")  # o el archivo donde tengas tu código UMB

# 1. Verifica que la densidad integra a 1
test_that("dUMB integrates to 1", {
  res <- integrate(function(x) dUMB(x, mu = 1), lower = 1e-6, upper = 1)
  expect_true(abs(res$value - 1) < 1e-4)
})

# 2. Verifica que la densidad es siempre positiva
test_that("dUMB is always positive", {
  x <- seq(0.01, 0.99, length.out = 100)
  expect_true(all(dUMB(x, mu = 1) > 0))
})

# 3. Comparar derivada manual vs numérica para sigma (mu en tu caso)
test_that("Manual vs numeric derivative of log-density matches", {
  x <- c(0.2, 0.4, 0.7)
  mu <- 1
  manual <- dldm(x, mu)
  
  # Derivada computacional
  dldm_compu <- function(x, mu) {
    dm <- numericDeriv(expr = quote(dUMB(x, mu, log=TRUE)), theta = "mu", rho = list(x = x, mu = mu))
    as.vector(attr(dm, "gradient"))
  }
  computa <- dldm_compu(x, mu)
  
  expect_equal(manual, computa, tolerance = 1e-4)
})

# 4. Verifica que la función pUMB está entre 0 y 1
test_that("pUMB is a valid CDF", {
  x <- seq(0.01, 0.99, length.out = 100)
  p <- pUMB(x, mu = 1)
  expect_true(all(p >= 0 & p <= 1))
  expect_true(is.unsorted(p) == FALSE)
})

# 5. Verifica que pUMB(qUMB) ≈ identidad
test_that("qUMB is inverse of pUMB", {
  p <- seq(0.01, 0.99, length.out = 100)
  expect_equal(pUMB(qUMB(p, mu = 1), mu = 1), p, tolerance = 1e-4)
})

# 6. Verifica que rUMB genera datos dentro de (0, 1)
test_that("rUMB generates values in (0, 1)", {
  x <- rUMB(1000, mu = 1)
  expect_true(all(x > 0 & x < 1))
})

# 7. Verifica estimación inicial de mu
test_that("estim_mu_UMB returns a positive estimate", {
  y <- rUMB(1000, mu = 1.5)
  mu_init <- estim_mu_UMB(y)[1]
  expect_true(mu_init > 0)
  expect_equal(mu_init, 1.5, tolerance = 0.2)
})

# 8. Verifica que gamlss() puede usar la familia UMB
test_that("GAMLSS model with UMB runs", {
  y <- rUMB(500, mu = 1)
  fit <- gamlss(y ~ 1, family = UMB())
  expect_s3_class(fit, "gamlss")
})


devtools::test()

