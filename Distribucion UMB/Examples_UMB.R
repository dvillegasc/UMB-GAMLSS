# Example 1
# Generating some random values with a known mu and fitting the model.

# Seed for reproducibility
set.seed(1234)

# Generate 300 random values from the UMB distribution.
y <- rUMB(n=300, mu=0.5)

#Fit the model using gamlss.
require(gamlss)
mod1 <- gamlss(y~1, family=UMB,
               control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod1, what="mu"))


# Example 2
# Generating random values for a regression model

# A function to simulate a data set with Y ~ UMB, where mu depends on x1

gendat <- function(n) {
  x1 <- runif(n)
  mu <- exp(0.5 - 2 * x1)
  y <- rUMB(n= n, mu= mu)
  data.frame(y= y, x1= x1)
}

# Seed for reproducibility
set.seed(1234)

# Generate the data
dat <- gendat(n= 200)

mod2 <- gamlss(y~ x1, family= UMB, data= dat,
               control= gamlss.control(n.cyc=5000, trace=FALSE))

summary(mod2)


# Example 3
# The first dataset measured the concentration of air pollutant CO in Alberta, 
# Canada from the Edmonton Central (downtown) Monitoring Unit (EDMU) station during 1995. 
# Measurements are listed for the period 1976–1995.
# Taken from Bicer et al. (2024) page 12.


data1 <- c(0.19,0.20,0.20,0.27,0.30,
           0.37,0.30,0.25,0.23,0.23,
           0.26,0.23,0.19,0.21,0.20,
           0.22,0.21,0.25,0.25,0.19) 


mod3 <- gamlss(data1 ~ 1, family= UMB, mu.link= "log",
               control= gamlss.control(n.cyc= 5000, trace= FALSE))


# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod3, what="mu"))

# Extraction of the log likelihood
logLik(mod3)

# Example 4
# The second data set measured air quality monitoring of the 
# annual average concentration of the pollutant benzo(a)pyrene (BaP). 
# The data were obtained from the Edmonton Central (downtown) Monitoring Unit (EDMU) location
# in Alberta, Canada, in 1995.
# Taken from Bicer et al. (2024) page 12.


data2 <- c(0.22,0.20,0.25,0.15,0.38,
           0.18,0.52,0.27,0.27,0.27,
           0.13,0.15,0.24,0.37,0.20) 


mod4 <- gamlss(data2 ~ 1, family= UMB, mu.link= "log",
               control= gamlss.control(n.cyc= 5000, trace= FALSE))



# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod4, what="mu"))

# Extraction of the log likelihood
logLik(mod4)


# Replicating figure 5 from Bicer et al. (2024)
# Hist and estimated pdf of Data-I and Data-II
mu1 <- 0.8452875 
mu2 <- 0.8593051 

par(mfrow = c(1, 2))

# Data-I
hist(data1, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 10),
     main = "Histogram of Data-I",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dUMB(x, mu = mu1), add = TRUE, 
      col = "blue", lwd = 2)

legend("topright", legend = c("UMB"), 
       col = c("blue"), lwd = 2, bty = "n")

# Data-II
hist(data2, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 6),
     main = "Histogram of Data-II",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dUMB(x, mu = mu2), add = TRUE, 
      col = "blue", lwd = 2)

legend("topright", 
       legend = c("UMB"), 
       col = c("blue"), 
       lwd = 2, 
       bty = "n")

par(mfrow = c(1, 1))

# Example 5
# The third dataset measured the concentration of sulphate 
# in Calgary from 31 different periods during 1995.
# Taken from Bicer et al. (2024) page 13.


data3 <- c(0.048,0.013,0.040,0.082,0.073,0.732,0.302,0.728,0.305,0.322,
           0.045,0.261,0.192,0.357,0.022,0.143,0.208,0.104,0.330,0.453,
           0.135,0.114,0.049,0.011,0.008,0.037,0.034,0.015,0.028,0.069,0.029) 


mod5 <- gamlss(data3 ~ 1, family= UMB, mu.link= "log",
               control= gamlss.control(n.cyc= 5000, trace= FALSE))



# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod5, what="mu"))

# Extraction of the log likelihood
logLik(mod5)

# Example 6
# The fourth dataset measured the concentration of pollutant CO in Alberta, Canada
# from the Calgary northwest (residential) monitoring unit (CRMU) station during 1995.
# Measurements are listed  for the period 1976-95.
# Taken from Bicer et al. (2024) page 13.


data4 <- c(0.16,0.19,0.24,0.25,0.30,0.41,0.40,
           0.33,0.23,0.27,0.30,0.32,0.26,0.25,
           0.22,0.22,0.18,0.18,0.20,0.23) 


mod6 <- gamlss(data4 ~ 1, family= UMB, mu.link= "log",
               control= gamlss.control(n.cyc= 5000, trace= FALSE))



# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod6, what="mu"))

# Extraction of the log likelihood
logLik(mod6)


# Replicating figure 6 from Bicer et al. (2024)
# Hist and estimated pdf of Data-III and Data-IV
mu3 <- 1.582003
mu4 <- 0.8161202 

par(mfrow = c(1, 2))

# Data-III
hist(data3, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 10),
     main = "Histogram of Data-III",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dUMB(x, mu = mu3), add = TRUE, 
      col = "blue", lwd = 2)

legend("topright", legend = c("UMB"), 
       col = c("blue"), lwd = 2, bty = "n")

# Data-IV
hist(data4, freq = FALSE,
     xlim = c(0, 1.0), ylim = c(0, 6),
     main = "Histogram of Data-IV",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")

curve(dUMB(x, mu = mu4), add = TRUE, 
      col = "blue", lwd = 2)

legend("topright", 
       legend = c("UMB"), 
       col = c("blue"), 
       lwd = 2, 
       bty = "n")

par(mfrow = c(1, 1))


# David Villegas Ceballos
