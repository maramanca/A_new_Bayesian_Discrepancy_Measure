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

#### Example 6 ####
rm(list = ls())

x <- c(1294, 1279, 1274, 1264, 1263, 
       1254, 1251, 1251, 1248, 1240, 
       1232, 1220, 1218, 1210)

y <- c(1284, 1272, 1256,
       1254, 1242, 1274,
       1264, 1256, 1250)

m <- length(x)
n <- length(y)

s1 <- sqrt(sum(x^2)/m-(mean(x))^2)
s2 <- sqrt(sum(y^2)/n-(mean(y))^2)

# Hyperparameters
eta1 <- mean(x)
nu1 <- m
alpha1 <- (m-1)/2
beta1 <- (m*s1^2)/2
eta2 <- mean(y)
nu2 <- n
alpha2 <- (n-1)/2
beta2 <- (n*s2^2)/2

# [A] Test on the means
mu1 <- seq(1230, 1275, 0.01)
mu2 <- seq(1230, 1275, 0.01)
G <- expand.grid(mu1,mu2)
G1 <- G[which(G$Var1>G$Var2),]

# Function necessary to make the plot
rA <- function(x){
  out=((1/beta((2*alpha1)/2,1/2))*sqrt(((nu1*alpha1)/beta1)/(2*alpha1))*(1+(((nu1*alpha1)/beta1)/(2*alpha1))*(x[1]-eta1)^2)^(-(2*alpha1)/2-1/2))*
    ((1/beta((2*alpha2)/2,1/2))*sqrt(((nu2*alpha2)/beta2)/(2*alpha2))*(1+(((nu2*alpha2)/beta2)/(2*alpha2))*(x[2]-eta2)^2)^(-(2*alpha2)/2-1/2))
  return(out)
}
 
# First Plot
p1 <- gcf(rA, lines_only=TRUE, bar=T, xlim = c(1230, 1275), ylim = c(1245, 1275))+
  stat_function(data=G, fun = function(Var1) (Var1), col = "red3")+
  labs(title="[A] Comparison between means", x=TeX("$\\mu_{1}"),
       y=TeX("$\\mu_{2}"),legend="Density")+
  annotate(geom='text', x=1232, y=1273, label=TeX("$\\Theta_a", output='character'), parse=TRUE)+
  annotate(geom='text', x=1273, y=1232, label=TeX("$\\Theta_b", output='character'), parse=TRUE)+
  annotate(geom='text', x=1269, y=1271, label=TeX("$\\Theta_{H_A}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(1230,1275,5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(1245,1275,5)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")

# Function necessary to evaluate the integral
g_post_A <- function(x1,x2){
  out=((1/beta((2*alpha1)/2,1/2))*sqrt(((nu1*alpha1)/beta1)/(2*alpha1))*(1+(((nu1*alpha1)/beta1)/(2*alpha1))*(x1-eta1)^2)^(-(2*alpha1)/2-1/2))*
    ((1/beta((2*alpha2)/2,1/2))*sqrt(((nu2*alpha2)/beta2)/(2*alpha2))*(1+(((nu2*alpha2)/beta2)/(2*alpha2))*(x2-eta2)^2)^(-(2*alpha2)/2-1/2))
  return(out)
}

# BDM
ymax <- function(x1) x1
delta_H_A <- round(1-2*(integral2(fun=g_post_A, xmin=1230, 
                                  xmax=1275, ymin=1230, 
                                  ymax=ymax)$Q),3)
delta_H_A


# [B] Test on the precisions
phi1 <- seq(0, 0.015, 0.0001)
phi2 <- seq(0, 0.015, 0.0001)
Grid <- expand.grid(phi1,phi2)
Grid1 <- Grid[which(Grid$Var1>Grid$Var2),]

r_B <- function(phi){
  out = dgamma(x = phi[1], shape = alpha1, scale = 1/beta1)*
    dgamma(x = phi[2], shape = alpha2, scale = 1/beta2)
  return(out)
}

# Second Plot
p2 <- gcf(r_B, lines_only=TRUE, bar=T, xlim = c(0, 0.008), ylim = c(0, 0.015))+
  stat_function(data=Grid, fun = function(Var1) (Var1), col = "red3")+
  labs(title="[B] Comparison between precisions", x=TeX("$\\phi_{1}"),
       y=TeX("$\\phi_{2}"),legend="Density")+
  annotate(geom='text', x=0.001, y=0.014, label=TeX("$\\Theta_c", output='character'), parse=TRUE)+
  annotate(geom='text', x=0.007, y=0.001, label=TeX("$\\Theta_d", output='character'), parse=TRUE)+
  annotate(geom='text', x=0.0075, y=0.008, label=TeX("$\\Theta_{H_B}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,0.008,0.001)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,0.015,0.001)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")

# Function necessary to evaluate the integral
g_post_B <- function(phi1,phi2){
  out = dgamma(x = phi1, shape = alpha1, scale = 1/beta1)*
    dgamma(x = phi2, shape = alpha2, scale = 1/beta2)
  return(out)
}

# BDM
ymax <- function(x1) x1
delta_H_A <- round(1-2*(integral2(fun=g_post_B, xmin=0, 
                                  xmax=0.015, ymin=0, 
                                  ymax=ymax)$Q),3)
delta_H_A

# Plot
x11()
grid.arrange(p1, p2, nrow = 1)


#### Example 7 ####
rm(list = ls())

# Sample 1
x <- c(0.8, 1.1, 1.2, 1.4, 1.8, 2, 4, 5, 8)
n1 <- length(x)
x_mean <- mean(x)
g1 <- geoMean(x, na.rm = FALSE)

# Sample 2
y <- c(0.25, 0.65, 0.81, 0.86, 0.87, 0.95, 1.06, 1.32, 2.12, 3.90, 4.03, 6.86)
n2 <- length(y)
y_mean <- mean(y)
g2 <- geoMean(y, na.rm = FALSE)

# Function necessary to evaluate the integral
g_posterior <- function(alpha1, alpha2){
  out = ((gamma(n1*alpha1)*gamma(n2*alpha2))/((gamma(alpha1)^n1)*gamma(alpha2)^n2))*
    (((g1/(n1*x_mean))^(n1*alpha1))*((g2/(n2*y_mean))^(n2*alpha2)))*
    sqrt(alpha1*psigamma(alpha1, 1) - 1)*
    sqrt(alpha2*psigamma(alpha2, 1) - 1)
  return(out)
}

# Normalizing constant
c <- 1/integral2(fun = g_posterior, xmin = 0.001, xmax = 9, 
                 ymin = 0.001, ymax = 9)$Q

# Integral
y_max <- function(x1) x1
int <- integral2(fun = g_posterior, xmin = 0.001, xmax = 9, 
                 ymin = 0.001, ymax = y_max)$Q*c

# BDM
delta_H <- round(1-2*(1-int), 3)
delta_H

# Plot
alpha1_es <- seq(0.01, 10, 0.01)
alpha2_es <- seq(0.01, 10, 0.01)

G_alpha <- expand.grid(alpha1_es,alpha2_es)

# Function necessary to make the plot
r <- function(x){
  prop = ((gamma(n1*x[1])*gamma(n2*x[2]))/((gamma(x[1])^n1)*gamma(x[2])^n2))*
    (((g1/(n1*x_mean))^(n1*x[1]))*((g2/(n2*y_mean))^(n2*x[2])))*
    sqrt(x[1]*psigamma(x[1], 1) - 1)*
    sqrt(x[2]*psigamma(x[2], 1) - 1)
  
  return(prop)
}

# Plot
x11()
gcf(r, lines_only=TRUE, bar=T, xlim = c(0.01, 4), ylim = c(0.01, 4))+
  stat_function(data=G_alpha, fun = function(Var1) Var1, col = "red3")+
  labs(title=TeX("Hypothesis $\\alpha_1 = \\alpha_2"), x=TeX("$\\alpha_1"),
       y=TeX("$\\alpha_2"),legend="Density")+
  annotate(geom='text', x=0.5, y=3.5, label=TeX("$\\Theta_{b}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.5, y=0.5, label=TeX("$\\Theta_{a}", output='character'), parse=TRUE)+
  annotate(geom='text', x=3.9, y=3.7, label=TeX("$\\Theta_{H}", output='character'), parse=TRUE)+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,4,0.5), limits = c(0,4)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5), legend.position = "none")


