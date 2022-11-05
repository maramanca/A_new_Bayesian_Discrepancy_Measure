#### Example 9 - Wald ####

# Function that computes the BDM 
d_H = function(Int){
  
  if(Int < 0.5){
    out = 1-2*Int
  }
  else {
    out = 1-2*(1-Int)
  }
}

n <- 8
mea <- 4.2
nu_0 <- 5

# Hypotheses
mu_H1 <- 2.5
mu_H2 <- 12

# Log-posterior distribution
log_posterior = function(theta,n,mea,nu_0){
  
  d.1 = (-3/2)*log(theta)
  d.2 = (-(n*mea)/2)*(nu_0/(theta^2))
  d.3 = n*(nu_0/theta)

  out = d.1+d.2+d.3
  
  return(out)
}

### Metropolis-Hastings algorithm ###
set.seed(5)
MH = function(S, B, theta0, n, mea, nu_0, sd, thin){
  
  R = B+S
  mat = rep(NA, R)
  acc = rep(0,R)
  
  theta_chain = theta0
  
  for (i in 1:R) {
    # Sampling from the proposal
    l_theta_star = rnorm(1, mean = log(theta_chain), sd = sd)
    theta_star = exp(l_theta_star)
    
    # Compute r
    l1 = log_posterior(theta = theta_star, n = n, mea = mea, nu_0 = nu_0)
    l2 = log_posterior(theta = theta_chain, n = n, mea = mea, nu_0 = nu_0)
    r = min(1, exp(l1 - l2))
    
    # Draw U from the uniform distribution in (0,1)
    u = runif(1)
    # Accept or reject theta.s
    if(r >= u) {
      theta_chain = theta_star
      acc[i] = 1
    }
    mat[i] = theta_chain
  }
  m1 <- mat[(B+1):R]
  m2 <- m1[seq(1, length(m1), by = thin)]
  
  acc1 <- acc[(B+1):R]
  acc2 <- acc1[seq(1, length(acc1), by = thin)]
  
  # Output
  list(values = m2, acc_rate=sum(acc2)/length(acc2))
}

# Starting values
theta0 = 4.80668
sd = 0.6

# Simulation
S = 320000
B = 100000

posterior_sample = MH(S = S, B = B, theta0 = theta0, n = n, mea = mea,
                      nu_0 = nu_0, sd = sd , thin = 3)

# Acceptance rate
posterior_sample$acc_rate

# Chain
mu <- posterior_sample$values

# Traceplot
x11()
plot(mu, type="l")

# Histogram
mean=mean(mu)
median=median(mu)
x11()
{hist(mu,3000, freq=F,xlab="mu",main="a) Histogram of mu")
  abline(v=mean, col="red",lwd=2)
  abline(v=median, col="green",lwd=2,lty=3)}

x11()
# ACF
acf(mu,main="ACF for mu",lag.max = 10000)

# Cumulative sums
x11()
plot(cumsum(mu) / seq_along(mu), type = 'l', xlab = 't',
     ylab = 'X_t', main = 'Ergodic Mean plot mu', col="red")



#### delta_H with simulation ####
N = length(mu)

I1 = sum(mu < mu_H1)/(N)
I1
delta_H_1 = round(d_H(Int = I1),3)
delta_H_1

I2 = sum(mu < mu_H2)/(N)
I2
delta_H_2 = round(d_H(Int = I2),3)
delta_H_2

# Results
c(delta_H_1, delta_H_2)

#### delta_H with integral ####
g_posterior <- function(param1){
  out = sqrt(1/param1^3)*
    exp(((-nu_0*n*mea)/(2*param1^2))+
          ((n*nu_0)/param1))
}

# Normalizing constant
c <- integrate(f=g_posterior, lower=0.0000001, upper=50000)$value
c

int <- (1/c) * integrate(f=g_posterior, lower=0.0001, upper=mu_H1)$value
int
delta_H_1 <- round(d_H(int),3)
delta_H_1

int2 <- (1/c) * integrate(f=g_posterior, lower=0.0001, upper=mu_H2)$value
int2
delta_H_2 <- round(d_H(int2),3)
delta_H_2



#### fbst package (flat reference) ####
library(fbst)
resFlat1 = fbst(mu, 
               nullHypothesisValue=mu_H1, 
               dimensionTheta=1, 
               dimensionNullset=1, 
               dim = 1, gridSize = 1000)

resFlat2 = fbst(mu, 
                nullHypothesisValue=mu_H2, 
                dimensionTheta=1, 
                dimensionNullset=1, 
                dim = 1, gridSize = 1000)

# Results
c(round(resFlat1$eValue,3),
  round(resFlat2$eValue,3))

#### fbst package (JEFFREYS reference) ####
# JEFFREYS reference function
djeffr <- function(x){
  out = (1/x)^(1/3)
  return(out)
}

resFlatSim_A = fbst(
  posteriorDensityDraws = mu, 
  nullHypothesisValue = mu_H1,
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1,
  FUN = djeffr)
resFlatSim_A$eValue

resFlatSim_B = fbst(
  posteriorDensityDraws = mu, 
  nullHypothesisValue = mu_H2, 
  dimensionTheta = 1, 
  dimensionNullset = 1,
  dim = 1,
  FUN = djeffr)
resFlatSim_B$eValue

# Results
round(c(resFlatSim_A$eValue, 
        resFlatSim_B$eValue), 3)