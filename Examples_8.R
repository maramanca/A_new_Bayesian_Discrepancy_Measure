library(ggplot2)
library(invgamma)
library(NormalGamma)
library(ContourFunctions)
library(rgl)
library(latex2exp)
library(pracma)
library(EnvStats)
library(fbst)

theta <- seq(.01, 2.3, .001)

alpha_A <- 6      # n
delta_A <- 7.2    # t_n
alpha_B <- 12      # n
delta_B <- 14.4    # t_n
alpha_C <- 24     # n
delta_C <- 28.8   # t_n

#### fbst by hand for theta_H = 2.4 (flat reference) ####
set.seed(1)

theta_H <- 2.4

# [A]
theta_prime_A <- theta[which(round(dinvgamma(theta, shape = alpha_A, rate = delta_A),3)==
                                 round(dinvgamma(theta_H, shape = alpha_A, rate = delta_A),3))][1]

Ev_A <- pinvgamma(theta_H, shape = alpha_A, rate = delta_A) -
  pinvgamma(theta_prime_A, shape = alpha_A, rate = delta_A)
Ev_A <- round(Ev_A,3)

# [B]
theta_prime_B <- theta[which(round(dinvgamma(theta, shape = alpha_B, rate = delta_B),2)==
                                 round(dinvgamma(theta_H, shape = alpha_B, rate = delta_B),2))][1]

Ev_B <- pinvgamma(theta_H, shape = alpha_B, rate = delta_B) -
  pinvgamma(theta_prime_B, shape = alpha_B, rate = delta_B) 
Ev_B <- round(Ev_B,3)
  
# [C]
theta_prime_C <- theta[which(round(dinvgamma(theta, shape = alpha_C, rate = delta_C),3)==
                                   round(dinvgamma(theta_H, shape = alpha_C, rate = delta_C),3))][1]

theta_prime_C

Ev_C <- pinvgamma(theta_H, shape = alpha_C, rate = delta_C) -
  pinvgamma(theta_prime_C, shape = alpha_C, rate = delta_C)
Ev_C <- round(Ev_C,3)

c(Ev_A, Ev_B, Ev_C)

#### fbst by hand for theta_H = 0.7 (flat reference) ####
set.seed(1)

theta_H <- 0.7

# [A]
theta_prime_A <- theta[which(round(dinvgamma(theta, shape = alpha_A, rate = delta_A),3)==
                               round(dinvgamma(theta_H, shape = alpha_A, rate = delta_A),3))][2]

Ev_A <- - pinvgamma(theta_H, shape = alpha_A, rate = delta_A) +
  pinvgamma(theta_prime_A, shape = alpha_A, rate = delta_A)
Ev_A <- round(Ev_A,3)

# [B]
theta_prime_B <- theta[which(round(dinvgamma(theta, shape = alpha_B, rate = delta_B),3)==
                               round(dinvgamma(theta_H, shape = alpha_B, rate = delta_B),3))][2]

Ev_B <- - pinvgamma(theta_H, shape = alpha_B, rate = delta_B) +
  pinvgamma(theta_prime_B, shape = alpha_B, rate = delta_B) 
Ev_B <- round(Ev_B,3)

# [C]
theta_prime_C <- theta[which(round(dinvgamma(theta, shape = alpha_C, rate = delta_C),3)==
                               round(dinvgamma(theta_H, shape = alpha_C, rate = delta_C),3))][2]

theta_prime_C

Ev_C <- - pinvgamma(theta_H, shape = alpha_C, rate = delta_C) +
  pinvgamma(theta_prime_C, shape = alpha_C, rate = delta_C)
Ev_C <- round(Ev_C,3)

c(Ev_A, Ev_B, Ev_C)




#### fbst package (flat reference) ####
set.seed(1)
sample_A <- rinvgamma(n = 1000000, shape = alpha_A, rate = delta_A)
sample_B <- rinvgamma(n = 1000000, shape = alpha_B, rate = delta_B)
sample_C <- rinvgamma(n = 1000000, shape = alpha_C, rate = delta_C)

# FLAT reference function
resFlatSim_A = fbst(
  posteriorDensityDraws = sample_A, 
  nullHypothesisValue = theta_H,
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1)
round(resFlatSim_A$eValue,3) 

resFlatSim_B = fbst(
  posteriorDensityDraws = sample_B, 
  nullHypothesisValue = theta_H, 
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1)
round(resFlatSim_B$eValue,3) 

resFlatSim_C = fbst(
  posteriorDensityDraws = sample_C, 
  nullHypothesisValue = theta_H, 
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1)
round(resFlatSim_C$eValue,3)

round(c(resFlatSim_A$eValue, 
        resFlatSim_B$eValue, 
        resFlatSim_C$eValue), 3)

# Plots
x11()
par(mfrow = c(1,3))
plot(resFlatSim_A)
plot(resFlatSim_B)
plot(resFlatSim_C)

#### fbst package (JEFFREYS reference) ####
# JEFFREYS reference function
djeffr <- function(x){
  out = 1/x
  return(out)
}

resFlatSim_A = fbst(
  posteriorDensityDraws = sample_A, 
  nullHypothesisValue = theta_H,
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1,
  FUN = djeffr)
resFlatSim_A$eValue

resFlatSim_B = fbst(
  posteriorDensityDraws = sample_B, 
  nullHypothesisValue = theta_H, 
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1,
  FUN = djeffr)
resFlatSim_B$eValue

resFlatSim_C = fbst(
  posteriorDensityDraws = sample_C, 
  nullHypothesisValue = theta_H, 
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1,
  FUN = djeffr)
resFlatSim_C$eValue

# Results
round(c(resFlatSim_A$eValue, 
        resFlatSim_B$eValue, 
        resFlatSim_C$eValue), 3)

# Plots
x11()
par(mfrow = c(1,3))
plot(resFlatSim_A)
plot(resFlatSim_B)
plot(resFlatSim_C)



#### BDM ####
set.seed(1)

# Posterior distributions
g_post_A <- dinvgamma(x = theta, shape = alpha_A, rate = delta_A)
g_post_B <- dinvgamma(x = theta, shape = alpha_B, rate = delta_B)
g_post_C <- dinvgamma(x = theta, shape = alpha_C, rate = delta_C)

# Medians
m1_A <- qinvgamma(0.5, shape = alpha_A, rate = delta_A)
m1_B <- qinvgamma(0.5, shape = alpha_B, rate = delta_B)
m1_C <- qinvgamma(0.5, shape = alpha_C, rate = delta_C)

# Discrepancy measure
delta_H = function(m1, theta_H, theta, alpha, delta){
  
  if(m1 < theta_H){
    out = 1-2*(1-pinvgamma(theta_H, shape = alpha, rate = delta))
  }
  else {
    out = 1-2*(pinvgamma(theta_H, shape = alpha, rate = delta))
  }
}

delta_H_A <- delta_H(m1_A, theta_H = theta_H,
                     theta = theta, alpha = alpha_A, delta = delta_A)

delta_H_B <- delta_H(m1_A,theta_H = theta_H,
                     theta = theta, alpha = alpha_B, delta = delta_B)

delta_H_C <- delta_H(m1_A,theta_H = theta_H,
                     theta = theta, alpha = alpha_C, delta = delta_C)

round(c(delta_H_A, delta_H_B, delta_H_C), 3)

