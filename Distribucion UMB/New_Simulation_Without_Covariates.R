library(gamlss)

# To perform the simulation -----------------------------------------------
if (!dir.exists("Simuls")) {
  dir.create("Simuls")
}


library("parSim")

parSim(
  ### SIMULATION CONDITIONS
  n = c(50, 100, 500, 1000),
  mu = c(0.25, 0.5, 1, 2),
  
  reps = 1000,                         # repetitions
  write = TRUE,                       # Writing to a file
  name = "Simuls/sim_without_covariates_07",  # Name of the file
  nCores = 1,                         # Number of cores to use
  
  expression = {
    # True parameter values
    y <- rUMB(n=n, mu)
    
    mod <- try(gamlss(y~1, family=UMB,
                      control=gamlss.control(n.cyc=1000, trace=FALSE)),
               silent=TRUE)
    
    if (class(mod)[1] == "try-error") {
      mu_hat    <- NA
    }
    else {
      mu_hat    <- exp(coef(mod, what="mu"))
    }
    
    # Results list:
    Results <- list(
      mu_hat = mu_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

archivos <- list.files(pattern = "^sim_without_cov.*\\.txt$", 
                       path="Simuls",
                       full.names = TRUE)

archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)

datos$case <- with(datos, ifelse(mu==0.25, 1, 
                                 ifelse(mu==0.5, 2,
                                        ifelse(mu==1, 3,
                                               ifelse(mu==2, 4, 5)))))

datos$case <- as.factor(datos$case)

# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.10 # percentage of values to be trimmed

dat <- datos %>% group_by(n, mu, case) %>% 
  summarise(nobs = n(),
            mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
            mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE),
            bias_mu = mean(mu_hat-mu, trim=trim, na.rm=TRUE),
  )

dat

# Plots
p1 <- ggplot(dat, aes(x=n, y=bias_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", mu))) +
  ylim(min(dat$bias_mu), 0.0019)

p1

p2 <- ggplot(dat, aes(x=n, y=mse_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", mu)))

p2

ggsave(filename="Figs/bias_mse_simul1.pdf", width=12, height=6,
       plot=p1+p2)


# Tables

trim <- 0.10 # percentage of values to be trimmed

dat <- datos %>% group_by(n, mu) %>% 
  summarise(nobs = n(),
            mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
            ab_mu = mean(abs(mu_hat-mu), trim=trim, na.rm=TRUE),
            mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE)
  )

dat


dat |> filter(mu == 1) |> 
  select(mean_mu, ab_mu, mse_mu) -> a
a[, -1]

library(xtable)
xtable(a[, -1])
