library(ggplot2)
library(invgamma)
library(NormalGamma)
library(ContourFunctions)
library(rgl)
library(latex2exp)
library(pracma)
library(EnvStats)
library(grid)
library(gridExtra)
library(e1071)

#### Example 3 ####
rm(list = ls())

x <- c(0.8, 1.1, 1.2, 1.4, 1.8, 2, 4, 5, 8)
n <- length(x)

# Statistics
x_mean <- mean(x)
g <- geoMean(x, na.rm = FALSE)

# Hypothesis
alpha_H <- 2.5
mu_H <- 6
sigma2_H <- 2

# Function necessary to evaluate the integral
g_posterior = function(alpha, beta){
  
  out = (((beta^alpha/gamma(alpha))*(g^alpha)*exp(-x_mean*beta))^n)*
    (sqrt(alpha*psigamma(alpha, 1) - 1)/beta)
  
  return(out)
}

# Function necessary to make the plot
r = function(x){
  out = (((x[2]^x[1]/gamma(x[1]))*(g^x[1])*exp(-x_mean*x[2]))^n)*
    (sqrt((x[1]*psigamma(x[1],1))-1)/x[2])
  
  return(out)
}

# Normalizing constant
c = 1/integral2(fun=g_posterior, xmin=0.01, xmax=100, 
                ymin=0.01, ymax=100)$Q

#### ALPHA
int_ALPHA <- integral2(fun = g_posterior, xmin = alpha_H, 
                       xmax = 100, ymin = 0.01,
                       ymax = 100)$Q*c
delta_H_ALPHA <- round(1 - 2*int_ALPHA, 3)
delta_H_ALPHA

# First Plot
p1 <- gcf(r, lines_only=TRUE, bar=T, xlim = c(0.01, 4), ylim = c(0.01, 2))+
  geom_vline(xintercept = alpha_H, col="red3")+
  labs(title=TeX("Hypothesis $\\alpha = 2.5"), x=TeX("$\\alpha"),
       y=TeX("$\\beta"),legend="Density")+
  theme_classic()+
  annotate(geom='text', x=0.6, y=1.5, label=TeX("$\\Theta_{a}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.5, y=0.2, label=TeX("$\\Theta_{b}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")


#### MU
beta_MU <- function(x1) x1/mu_H
int_MU <- integral2(fun=g_posterior, xmin=0.01, xmax=4, 
                    ymin=0.01, ymax=beta_MU)$Q*c
delta_H_MU <- round(1-2*int_MU, 3)
delta_H_MU

# Second Plot
alpha_es <- seq(0.01, 4, 0.01)
beta_es <- seq(0.01, 2, 0.01)

G <- expand.grid(alpha_es,beta_es)

p2 <- gcf(r, lines_only=TRUE, bar=T, xlim = c(0.01, 4), ylim = c(0.01, 2))+
  stat_function(data=G, fun = function(Var1) Var1/mu_H, col="red3")+
  labs(title=TeX("Hypothesis $\\mu = 6"), x=TeX("$\\alpha"),
       y=TeX("$\\beta"),legend="Density")+
  annotate(geom='text', x=0.6, y=1.5, label=TeX("$\\Theta_{c}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.5, y=0.2, label=TeX("$\\Theta_{d}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")


#### SIGMA
beta_SIGMA2 <- function(alpha) sqrt(alpha/sigma2_H)
int_SIGMA2 <- 1 - integral2(fun=g_posterior, xmin=0.01, xmax=10, 
                            ymin=0.01, ymax=beta_SIGMA2)$Q*c
delta_H_SIGMA2 <- round(1-2*int_SIGMA2, 3)
delta_H_SIGMA2

# Third Plot
p3 <- gcf(r, lines_only=TRUE, bar=T, xlim = c(0.01, 4), ylim = c(0.01, 2))+
  stat_function(data=G, fun = function(Var1) (sqrt(Var1/sigma2_H)), col="red3")+
  labs(title=TeX("Hypothesis $\\sigma^2 = 2"), x=TeX("$\\alpha"),
       y=TeX("$\\beta"),legend="Density")+
  theme_classic()+
  annotate(geom='text', x=0.6, y=1.5, label=TeX("$\\Theta_{e}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.5, y=0.2, label=TeX("$\\Theta_{f}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")

# Plot
x11()
grid.arrange(p1, p2, p3, nrow = 1)

c(delta_H_ALPHA, delta_H_MU, delta_H_SIGMA2)



#### Example 4 ####
rm(list = ls())

x_mean <- 17
s_square <- 1.6

# [A]
n <- 10
# [B]
m <- 40

# [A]  
eta1 <- x_mean          # = eta2
nu1 <- n
alpha1 <- (n-1)/2
beta1 <- (n*s_square)/2
# [B]
nu2 <- m
alpha2 <- (m-1)/2
beta2 <- (m*s_square)/2

mu <- seq(15, 19, 0.01)
phi <- seq(0.01, 1.5, 0.01)
psi_H <- 0.1

# [A]
G <- expand.grid(mu,phi)
G1 <- G[which(G$Var2 < 1/((psi_H^2)*(G$Var1)^2)),]

# Function necessary to make the plot
r <- function(x){
  out = (sqrt((nu1*x[2])/(2*pi)))*exp((-0.5*nu1*x[2])*(x[1]-eta1)^2)*
    (beta1^alpha1/gamma(alpha1))*(x[2])^(alpha1-1)*exp(-beta1*x[2])
  return(out)
}

p1 <- gcf(r, lines_only=TRUE, bar=T, xlim = c(15, 19), ylim = c(0, 1.5))+
  stat_function(data=G, fun = function(Var1) (1/((0.1^2)*(Var1^2))), col="red3")+
  labs(title="[A] n = 10", x=TeX("$\\mu"),
       y=TeX("$\\phi"),legend="Density")+
  annotate(geom='text', x=18.5, y=0.1, label=TeX("$\\Theta_b", output='character'), parse=TRUE)+
  annotate(geom='text', x=15.2, y=1.3, label=TeX("$\\Theta_a", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(15,19,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1.5,0.2), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")

# Function necessary to evaluate the integral
g_post <- function(x1,x2){
  out = (sqrt((nu1*x2)/(2*pi)))*exp((-0.5*nu1*x2)*((x1-eta1)^2))*
    ((beta1^alpha1)/gamma(alpha1))*((x2)^(alpha1-1))*exp(-beta1*x2)
  return(out)
}

# 3D plot
data <- data.frame(x=G1$Var1, y=G1$Var2, z=g_post(G1$Var1,G1$Var2))
plot3d(subset(data,T,c(x,y,z)),col="red3")

# BDM
ymax <- function(x1) 1/((psi_H^2)*(x1^2))
delta_H_A <- 1-2*(integral2(fun=g_post, xmin=15, xmax=19, 
                            ymin=0.01, ymax=ymax)$Q)
delta_H_A <- round(delta_H_A, 3)
delta_H_A

# [B]
# Function necessary to make the plot
r_B <- function(x){
  out = (sqrt((nu2*x[2])/(2*pi)))*exp((-0.5*nu2*x[2])*(x[1]-eta1)^2)*
    ((beta2^alpha2)/gamma(alpha2))*(x[2])^(alpha2-1)*exp(-beta2*x[2])
  return(out)
}

p2 <- gcf(r_B, lines_only=TRUE, bar=T, xlim = c(15, 19), ylim = c(0, 1.5))+
  stat_function(data=G, fun = function(Var1) (1/((0.1^2)*(Var1^2))), col="red3")+
  labs(title="[B] n = 40", x=TeX("$\\mu"),
       y=TeX("$\\phi"),legend="Density")+
  annotate(geom='text', x=18.5, y=0.1, label=TeX("$\\Theta_b", output='character'), parse=TRUE)+
  annotate(geom='text', x=15.5, y=1.3, label=TeX("$\\Theta_a", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(15,19,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1.5,0.2), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")

# Function necessary to evaluate the integral
g_post_B <- function(x1,x2){
  out = (sqrt((nu2*x2)/(2*pi)))*exp((-0.5*nu2*x2)*((x1-eta1)^2))*
    ((beta2^alpha2)/gamma(alpha2))*((x2)^(alpha2-1))*exp(-beta2*x2)
  return(out)
}

# 3D plot
data_B <- data.frame(x=G1$Var1, y=G1$Var2, z=g_post_B(G1$Var1,G1$Var2))
plot3d(subset(data_B,T,c(x,y,z)),col="red3")

# BDM
delta_H_B <- 1-2*(integral2(fun=g_post_B, xmin=15, xmax=19, 
                            ymin=0.01, ymax=ymax)$Q)
delta_H_B <- round(delta_H_B, 3)
delta_H_B

# Results
c(delta_H_A, delta_H_B)

# Plot
x11()
grid.arrange(p1, p2, nrow = 1)



#### Example 5 ####
rm(list = ls())

# Function that computes the BDM 
d_H = function(Int){
  
  if(Int < 0.5){
    out = 1-2*Int
  }
  else {
    out = 1-2*(1-Int)
  }
}

Data = c(1.01, 1.11, 1.13, 1.15, 1.16, 1.17, 
         1.17, 1.20, 1.52, 1.54, 1.54, 1.57, 
         1.64, 1.73, 1.79, 2.09, 2.09, 2.57, 
         2.75, 2.93, 3.19, 3.54, 3.57, 5.11, 5.62)

# Sample skewness
skewness(Data)

n <- length(Data)
t <- sum(Data)
m <- sum(1/Data)

# Hypothesis
gamma_H = 2

#### Computation delta_H with integral
g_posterior <- function(param1,param2){
  out = sqrt(((param2)^(n-1))/param1^3)*
    exp(-param2*((t/(2*param1^2))-(n/param1)+0.5*m))
}

# Normalizing constant
c <- integral2(fun=g_posterior, xmin=1, xmax=4,
               ymin=0.01, ymax=16)$Q
c

# Integral
ymax <- function(x1) (x1*9)/((gamma_H)^2)
int <- (1/c) * integral2(fun=g_posterior, xmin=1, xmax=4, 
                         ymin=0.01, ymax=ymax)$Q
int

# BDM
delta_H <- round(d_H(int),3)
delta_H

# Function necessary to make the plot
r <- function(x){
  out = sqrt(((x[2])^(n-1))/x[1]^3)*
    exp(-x[2]*((t/(2*x[1]^2))-(n/x[1])+0.5*m))
}

mu_es <- seq(1, 4, 0.01)
nu_es <- seq(0.01, 16, 0.01)
G <- expand.grid(mu_es,nu_es)
G1 <- G[which(G$Var2 < (3/(gamma_H))^2*G$Var1),]

# Plot
x11()
gcf(r, lines_only=TRUE, bar=T, xlim = c(1, 4), ylim = c(0.01, 16))+
  stat_function(data=G, fun = function(Var1) (3/(gamma_H))^2*Var1, col = "red3")+
  labs(x=TeX("$\\mu"),
       y=TeX("$\\nu"),legend="Density")+
  annotate(geom='text', x=1.5, y=12, label=TeX("$\\Theta_{a}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.2, y=4, label=TeX("$\\Theta_{b}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.95, y=8, label=TeX("$\\Theta_{H}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(1,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,16,2), limits = c(0,16)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")

