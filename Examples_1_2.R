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


#### Example 1 ####
rm(list = ls())

theta <- seq(.01, 4, .01)
theta_H <- 2.4

alpha_A <- 6      # n
delta_A <- 7.2    # t_n
alpha_B <- 12     # n
delta_B <- 14.4   # t_n
alpha_C <- 24     # n
delta_C <- 28.8   # t_n

# Posterior distributions
g_post_A <- dinvgamma(x = theta, shape = alpha_A, rate = delta_A)
g_post_B <- dinvgamma(x = theta, shape = alpha_B, rate = delta_B)
g_post_C <- dinvgamma(x = theta, shape = alpha_C, rate = delta_C)

# Medians
m1_A <- qinvgamma(0.5, shape = alpha_A, rate = delta_A)
m1_B <- qinvgamma(0.5, shape = alpha_B, rate = delta_B)
m1_C <- qinvgamma(0.5, shape = alpha_C, rate = delta_C)

# Single plots
data_A <- data.frame(theta=theta, g_posterior=g_post_A)
p1 <- ggplot(data_A, aes(x=theta,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
  geom_vline(xintercept=theta_H,size = 0.4,col="red3")+
  labs(title="[A] n = 6", x=TeX("$\\theta"),
       y=TeX("$g_1(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

data_B <- data.frame(theta=theta, g_posterior=g_post_B)
p2 <- ggplot(data_B, aes(x=theta,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  geom_vline(xintercept=theta_H,size = 0.4,col="red3")+
  labs(title="[B] n = 12", x=TeX("$\\theta"),
       y=TeX("$g_1(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

data_C <- data.frame(theta=theta, g_posterior=g_post_C)
p3 <- ggplot(data_C, aes(x=theta,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_C,size = 0.4,col="blue3", linetype="twodash")+
  geom_vline(xintercept=theta_H,size = 0.4,col="red3")+
  labs(title="[C] n = 24", x=TeX("$\\theta"),
       y=TeX("$g_1(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

# Labels: median and H
m1 <- TeX("$m_1$")
H <- TeX("$\\theta_H$")

# Plot
x11()
grid.arrange(p1, p2, p3, nrow = 1)
grid.text(m1, gp = gpar(cex=1), x=0.15, y = 0.82)
grid.text(H, gp = gpar(cex=1), x=0.23, y = 0.82)
grid.text(m1, gp = gpar(cex=1), x=0.48, y = 0.82)
grid.text(H, gp = gpar(cex=1), x=0.56, y = 0.82)
grid.text(m1, gp = gpar(cex=1), x=0.81, y = 0.82)
grid.text(H, gp = gpar(cex=1), x=0.895, y = 0.82)

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

delta_H_B <- delta_H(m1_B,theta_H = theta_H,
                     theta = theta, alpha = alpha_B, delta = delta_B)

delta_H_C <- delta_H(m1_C,theta_H = theta_H,
                     theta = theta, alpha = alpha_C, delta = delta_C)

round(c(delta_H_A, delta_H_B, delta_H_C),3)


#### The BDM trend when varying the MLE and the sample size n 
alpha <- seq(1, 50, 1)      # n

#### Plot 1 - Left
MLE1 <- c(0.8,1.2,1.6)  
delta1 <- matrix(NA, length(alpha), length(MLE1)) #t_n 
for (i in 1:length(MLE1)) {
  delta1[,i] <- MLE1[i]*alpha  # t_n
} 

m1_1 <- matrix(NA, length(alpha), length(MLE1)) 
BDM_1 <- matrix(NA, length(alpha), length(MLE1)) 
for (i in 1:length(alpha)) {
  for (j in 1:length(MLE1)) {
    # Median
    m1_1[i,j] <- qinvgamma(0.5, 
                           shape = alpha[i], 
                           rate = delta1[i,j])
    #BDM
    BDM_1[i,j] <- delta_H(m1_1[i,j], 
                          theta_H = theta_H,
                          alpha = alpha[i], 
                          delta = delta1[i,j])    
  }
}

data <- data.frame(n = rep(alpha, length(MLE1)), 
                   delta_H = c(BDM_1[,1], 
                               BDM_1[,2], 
                               BDM_1[,3]),
                   MLE1 = c(rep("0.8", length(alpha)),
                            rep("1.2", length(alpha)),
                            rep("1.6", length(alpha))))

p1 <- ggplot(data, aes(x = n,y = delta_H, color = MLE1))+
  geom_line(aes(color = MLE1))+
  geom_hline(yintercept = 0.95, size = 0.8, col = "red3")+
  scale_x_continuous(breaks = seq(0,50,2))+
  scale_y_continuous(breaks = seq(0,1,0.05), 
                     limits = c(0,1))+
  theme_classic()+
  scale_colour_hue(h = c(180, 270))+
  labs(title = "", x = TeX("n"),
       colour = "MLE",
       y = TeX("$\\delta_H"))


#### Plot 2 - Right 
MLE2 <- c(3.2,3.6,4)
delta2 <- matrix(NA, length(alpha), length(MLE1)) #t_n 
for (i in 1:length(MLE2)) {
  delta2[,i] <- MLE2[i]*alpha  # t_n
} 

m1_2 <- matrix(NA, length(alpha),length(MLE2)) 
BDM_2 <- matrix(NA, length(alpha), length(MLE2)) 
for (i in 1:length(alpha)) {
  for (j in 1:length(MLE2)) {
    # Median
    m1_2[i,j] <- qinvgamma(0.5, shape = alpha[i], rate = delta2[i,j])
    #BDM
    BDM_2[i,j] <- delta_H(m1_2[i,j], theta_H = theta_H,
                          alpha = alpha[i], delta = delta2[i,j])    
  }
}

data2 <- data.frame(n = rep(alpha, length(MLE2)), 
                    delta_H=c(BDM_2[,1], BDM_2[,2], BDM_2[,3]),
                    MLE2=c(rep("3.2", length(alpha)),
                           rep("3.6", length(alpha)),
                           rep("4", length(alpha))))
p2 <- ggplot(data2, aes(x = n, y = delta_H, color = MLE2))+
  geom_line(aes(color = MLE2))+
  geom_hline(yintercept = 0.95,size = 0.8, col = "red3")+
  scale_x_continuous(breaks = seq(0,50,2))+
  scale_y_continuous(breaks = seq(0,1,0.05), limits = c(0,1))+
  scale_color_hue(h = c(0, 90))+
  theme_classic()+
  labs(title = "", 
       colour = "MLE",
       x = TeX("n"),
       y = TeX("$\\delta_H"))

# Plot
x11()
grid.arrange(p1, p2, nrow = 1)



#### Example 2 ####
rm(list = ls())

tau <- seq(.01, 1, .00001)

x <- 7640       # successes number (M)
n <- 14928      # births total number 

tau1 <- 1/2
tau2 <- 13/25
tau3 <- 1050/2050
tau4 <- 23/45

alpha_star <- 1+x      # alpha_star = alpha+x
beta_star <- 1+n-x     # beta_star = beta+n-x

# Posterior distribution
g_post <- dbeta(x = tau, shape1 = alpha_star, 
                shape2 = beta_star)

# Plot
data <- data.frame(tau=tau, g_posterior=g_post)
ggplot(data, aes(x=tau,y=g_posterior))+
  geom_line(size = 0.4, col="blue3")+
  labs(title="", x=TeX("$\\tau"),
       y=TeX("$g_1(\\tau | \\alpha^* ,\\beta^*)"))+
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

# Median
m <- qbeta(0.5, shape1 = alpha_star, shape2 = beta_star)

# Discrepancy measure
delta_H = function(m, tau_H, tau, alpha, beta){
  
  if(m < tau_H){
    out = 1-2*(1-pbeta(tau_H, shape1 = alpha, shape2 = beta))
  }
  else {
    out = 1-2*(pbeta(tau_H, shape1 = alpha, shape2 = beta))
  }
}

delta_H_1 <- delta_H(m, tau_H = tau1, tau = tau, alpha = alpha_star, 
                     beta = beta_star)
delta_H_2 <- delta_H(m, tau_H = tau2, tau = tau, alpha = alpha_star, 
                     beta = beta_star)
delta_H_3 <- delta_H(m, tau_H = tau3, tau = tau, alpha = alpha_star, 
                     beta = beta_star)
delta_H_4 <- delta_H(m, tau_H = tau4, tau = tau, alpha = alpha_star, 
                     beta = beta_star)

round(c(delta_H_1, delta_H_2, delta_H_3, delta_H_4),3)