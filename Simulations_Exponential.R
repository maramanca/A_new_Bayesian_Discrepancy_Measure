library(invgamma)
library(EnvStats)
library(fbst)

# Increasing the memory available to R processes
memory.limit()
memory.limit(size=10^13)

set.seed(1)

# BDM
delta_H = function(m1, theta_H, theta, alpha, delta){
  
  if(m1 < theta_H){
    out = 1-2*(1-pinvgamma(theta_H, shape = alpha, rate = delta))
  }
  else {
    out = 1-2*(pinvgamma(theta_H, shape = alpha, rate = delta))
  }
}

# JEFFREYS reference function
djeffr <- function(x){
  out = 1/x
  return(out)
}

theta_values <- seq(.01, 5, .001)
theta_H <- 1.2 # hypothesis

S <- 50000     # number of simulations
n <- 10        # sample size
D <- 50000     # number of posterior draws

mea <- rep(NA, S)     # mean of the sample
alpha <- 0            # first parameter posterior
delta <- rep(NA, S)   # second parameter posterior
m1 <- rep(NA, S)      # median

x <- matrix(NA, n, S)      # sample
theta <- matrix(NA, D, S)  # posterior draws

e <- rep(NA, S)            # e-value with reference prop to 1
e_g <- rep(NA, S)          # e-value with Jeffreys' prior as reference
I <- rep(NA, S)            # Integral for the BDM
d_H <- rep(NA, S)          # BDM

for (i in 1:S) {
  x[,i] <- rexp(n = n, rate = 1/theta_H)     
  mea[i] <- mean(x[,i])
  
  alpha <- n
  delta[i] <- mea[i]*n
  
  theta[,i] <- rinvgamma(n = D, shape = alpha, rate = delta[i])
  
  # e-values
  e[i] = fbst(theta[,i], 
              nullHypothesisValue=theta_H, 
              dimensionTheta=1, 
              dimensionNullset=1, 
              dim = 1, gridSize = 1000)$eValue
  
  e_g[i] = fbst(theta[,i], 
                nullHypothesisValue = theta_H,
                dimensionTheta = 1, 
                dimensionNullset = 1,
                dim = 1,
                FUN = djeffr)$eValue
  
  # BDM
  m1[i] <- qinvgamma(0.5, shape = alpha, rate = delta[i])
  
  d_H[i] = delta_H(m1[i], 
                   theta_H = theta_H,
                   theta = theta_values, 
                   alpha = alpha, 
                   delta = delta[i])
  print(i)
}

# Summary for each evidence measure
summary(e)
summary(e_g)
summary(d_H)

# False positive rate for threshold = 0.90, 0.95, 0.99
c(round(sum(e>=0.9)/S,3),
  round(sum(e_g>=0.9)/S,3),
  round(sum(d_H>=0.9)/S,3))

c(round(sum(e>=0.95)/S,3),
  round(sum(e_g>=0.95)/S,3),
  round(sum(d_H>=0.95)/S,3))

c(round(sum(e>=0.99)/S,3),
  round(sum(e_g>=0.99)/S,3),
  round(sum(d_H>=0.99)/S,3))

# Plot of a given chain
x11()
plot(density(theta[,5000]))

# Tables in order to appreciate the differences 
# between the e-value and the BDM
# r prop 1                          # (e,delta_H)
n_11 = sum(e<0.95 & d_H<0.95)/S     # AA
n_22 = sum(e>=0.95 & d_H>=0.95)/S   # RR
n_12 = sum(e<0.95 & d_H>=0.95)/S    # AR
n_21 = sum(e>=0.95 & d_H<0.95)/S    # RA

# AA AR
# RA RR
matrix(data = round(c(n_11,n_12,n_21,n_22),3),2,2, 
       byrow = TRUE)

# r prop Jeffreys
n_11 = sum(e_g<0.95 & d_H<0.95)/S
n_22 = sum(e_g>=0.95 & d_H>=0.95)/S
n_12 = sum(e_g<0.95 & d_H>=0.95)/S
n_21 = sum(e_g>=0.95 & d_H<0.95)/S

# AA AR
# RA RR
matrix(data = round(c(n_11,n_12,n_21,n_22),3),2,2, 
       byrow = TRUE)

