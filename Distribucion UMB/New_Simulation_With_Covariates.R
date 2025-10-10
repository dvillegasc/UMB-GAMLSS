require(gamlss)


# To perform the simulation -----------------------------------------------
if (!dir.exists("Simuls")) {
  dir.create("Simuls")
}


library("parSim")

gendat <- function(n) {
  x1 <- runif(n)
  mu    <- exp(0.33 + 0.11 * x1) # 1.5 approximately
  y <- rUMB(n=n, mu=mu)
  data.frame(y=y, x1=x1)
}

parSim(
  ### SIMULATION CONDITIONS
  n = c(30, 50, 100, 200),
  
  reps = 1000,                      # repetitions
  write = TRUE,                    # Writing to a file
  name = "Simuls/sim_with_covariates_nA",  # Name of the file
  nCores = 1,                      # Number of cores to use
  
  expression = {
    # True parameter values
    dat <- gendat(n=n)
    
    mod <- gamlss(y~x1, family=UMB, data=dat,
                  control=gamlss.control(n.cyc=1000, trace=FALSE))
    
    beta_0_hat  <- coef(mod, what="mu")[1]
    beta_1_hat  <- coef(mod, what="mu")[2]
    
    # Results list:
    Results <- list(
      beta_0_hat = beta_0_hat,
      beta_1_hat = beta_1_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

archivos <- list.files(pattern = "^sim_with_cov.*\\.txt$", 
                       path="Simuls",
                       full.names = TRUE)
archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)


prop.table(table(datos$error == TRUE))


# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.10

dat <- datos %>% group_by(n) %>% 
  summarise(nobs = n(),
            
            bias_b0 = mean(beta_0_hat - (0.33), trim=trim, na.rm=TRUE),
            bias_b1 = mean(beta_1_hat - (0.11), trim=trim, na.rm=TRUE),
            
            mse_b0 = mean((beta_0_hat - (0.33))^2, trim=trim, na.rm=TRUE),
            mse_b1 = mean((beta_1_hat - (0.11))^2, trim=trim, na.rm=TRUE)
            
  )

dat

# Legend and colores
leyenda <- c(expression(hat(beta)[0]), 
             expression(hat(beta)[1]))

colores <- c("#F8766D", "#00BA38")


d <- pivot_longer(data=dat, 
                  cols=c("bias_b0", "bias_b1"),
                  names_to="Estimator",
                  values_to="value")

# Plots
p1 <- ggplot(d, aes(x=n, y=value, colour=Estimator)) +
  geom_line() + 
  ylab("Bias") + 
  scale_color_manual(labels=leyenda,
                     values=colores)

p1

d <- pivot_longer(data=dat, 
                  cols=c("mse_b0", "mse_b1"),
                  names_to="Estimator",
                  values_to="value")

p2 <- ggplot(d, aes(x=n, y=value, colour=Estimator)) +
  geom_line() + 
  ylab("MSE") + 
  scale_color_manual(labels=leyenda,
                     values=colores)

p2

ggsave(filename="Figs/bias_mse_simul2nA.pdf", width=12, height=6,
       plot=p1+p2)


