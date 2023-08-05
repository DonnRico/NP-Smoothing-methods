
####################################################################
#                         SCRIPT LSTAT2150                         #
#   Warnauts Aymeric 87031800             UCLouvain-LIDAM          #
#                         Von Sachs Rainer                         #
####################################################################


#####################################
# Packages installation and loading #
#####################################

library(EnvStats)
library(ggplot2)
theme_set(theme_bw())
library(latex2exp)
library(UsingR)
library(mixtools)
library(kableExtra)

#############################
# Generation of the mixture #
#############################

mixture <- function(alpha) {
  set.seed(87031800)
  w1 <- alpha
  w2 <- 1 - alpha
  w <- sample(c(1,2), size = 100, replace = TRUE, prob = c(w1,w2))
  X <- (w==1)*rnorm(50, 2, 1) + (w==2)*rnorm(50,6,1)
  return(X)
}


#########
# Task1 #
#########

#####################################
#    Histogram density estimator    #
#####################################

for (i in c(0,0.05,0.25,0.5,0.75,0.95,1)){
  X <- mixture(i)
  mix <- as.data.frame(X)
  p <- ggplot(mix, aes(x=X)) + geom_histogram(aes(y = ..density..), colour = 1,
                                              fill = "white", alpha = 0.4) + 
    xlim(0, 11) + 
    ylim(0,1) +
    geom_density(color = "red", fill = "red", alpha = 0.15) + 
    annotate("text", x = 1, y = 1, label = paste('alpha =',i)) + 
    labs(title = 
           TeX("Density of the Normal mixture ($\\mu_1 = 2, \\sigma_1 =1$) and 
               ($\\mu_2 = 6, \\sigma_2 =1$)"), 
         x = "sample")
  
  print(p)
}

#####################################
#      Kernel density estimator     #
#####################################


h <- 0.3   ## change bandwidth if you want till optimal is not determined
for (i in c(0,0.05,0.25,0.5,0.75,0.95,1)){
  png(filename = paste(i, "kernel.png"), width = 7, height = 6, units = 'in',
      res=600)
  X <- mixture(i)
  x <- seq(from = min(X), to = max(X), length = 100)  ## Points on the domain
  prDensEst <- function(x, X, h, K) mean(K((x - X)/h))/h 
  ## standardization of the kernel
  
  u2Knorm <- function(u) u^2 * dnorm(u) ## Gaussian kernel
  u2Kunif <- function(u) u^2 * (abs(u) <= 1) * 0.5 ## Uniform kernel
  u2Kepan <- function(u) u^2 * 0.75 * (1 - u^2) * (abs(u) <= 1) ## Epan kernel
  
  IntKnorm <- integrate(u2Knorm, -Inf, Inf)$value
  IntKunif <- integrate(u2Kunif, -Inf, Inf)$value
  IntKepan <- integrate(u2Kepan, -Inf, Inf)$value
  
  h_norm <- h/IntKnorm^0.5
  h_unif <- h/IntKunif^0.5
  h_epan <- h/IntKepan^0.5
  
  Knorm <- function(u) dnorm(u)
  Kunif <- function(u) (abs(u) <= 1) * 0.5
  Kepan <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)
  
  fnorm <- sapply(x, function(x) prDensEst(x, X, h_norm, Knorm)) 
  funif <- sapply(x, function(x) prDensEst(x, X, h_unif, Kunif))
  fepan <- sapply(x, function(x) prDensEst(x, X, h_epan, Kepan))
  
  plot(density(X, n = 100, from = min(X), to = max(X), 
               bw = h, kernel = "gaussian"),
       main = "Kernel Density, Implemented Kernel",
       col = 1, lwd = 2, xlim = c(-2,10), ylim = c(0,0.6))
  lines(x, fnorm, col = 2)
  lines(x, funif, col = 3)
  lines(x, fepan, col = 4)
  legend(7, 0.6, legend=c("Gaussian kernel", "Epan kernel", 
                          "Uniform kernel", "R density"),
         col=c("red", "blue", "green", "black"),lty=1, cex=0.8)
  legend(-2.2,0.6, legend = paste("alpha =",i))
  rug(X)
  dev.off()
}

par(mfrow = c(1, 2))
for (i in c(0.05,0.25,0.5,0.75,0.95,1)){
  X <- mixture(i)
  EM <- normalmixEM(X, lambda = i, mu = c(2, 6), sigma = c(1,1))
  plot(EM, density=TRUE, 
       main2= TeX("EM Density of the Normal mixture ($\\mu_1 = 2,
                  \\sigma_1 =1$) and ($\\mu_2 = 6, \\sigma_2 =1$)"), 
       xlab2=paste("alpha = ", i))
}

#####################################
#      Parametric EM generation     #
#####################################

X <- mixture(0.05) ## Change alpha to visualize EM estimations
EM <- normalmixEM(X, lambda = 0.25, mu = c(2, 6), sigma = c(1,1))


#####################
# Bias (parametric) #
#####################

EM$mu  ## perform this at each EM estimation for different alpha values
EM$mu - c(2,6)

EM$sigma
EM$sigma - c(1, 1)


#####################################
#  Plot comparison for each alpha   #
#####################################


X <- mixture(0.05) ## Change alpha to perfom this for different mixtures
True <- density(X, n = 100)
plot(True)
data <- cbind(True$x, True$y)
colnames(data) <- c("TrueVal", "TrueDens")
data


for (i in c(0.05,0.25,0.5,0.75,0.95,1)){
  png(filename = paste(i, "comparison.png"), width = 7, height = 6, units = 'in',
      res=600)
  X <- mixture(i)
  True <- density(X, n = 100)
  x <- seq(from = min(X), to = max(X), length = 100) ## Points on the domain
  prDensEst <- function(x, X, h, K) mean(K((x - X)/h))/h 
  ## standardization of the kernel
  
  u2Knorm <- function(u) u^2 * dnorm(u) ## Gaussian kernel
  u2Kunif <- function(u) u^2 * (abs(u) <= 1) * 0.5 ## Uniform kernel
  u2Kepan <- function(u) u^2 * 0.75 * (1 - u^2) * (abs(u) <= 1)
  
  IntKnorm <- integrate(u2Knorm, -Inf, Inf)$value
  IntKunif <- integrate(u2Kunif, -Inf, Inf)$value
  IntKepan <- integrate(u2Kepan, -Inf, Inf)$value
  
  h_norm <- h/IntKnorm^0.5
  h_unif <- h/IntKunif^0.5
  h_epan <- h/IntKepan^0.5
  
  Knorm <- function(u) dnorm(u)
  Kunif <- function(u) (abs(u) <= 1) * 0.5
  Kepan <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)
  
  fnorm <- sapply(x, function(x) prDensEst(x, X, h_norm, Knorm)) 
  funif <- sapply(x, function(x) prDensEst(x, X, h_unif, Kunif))
  fepan <- sapply(x, function(x) prDensEst(x, X, h_epan, Kepan))
  
  EM <- normalmixEM(X, lambda = i, mu = c(2, 6), sigma = c(1,1))
  w1 <- EM$lambda[1]
  w2 <- EM$lambda[2]
  w <- sample(c(1,2), size = 100, replace = TRUE, prob = c(w1,w2))
  X <- (w==1)*rnorm(50, EM$mu[1], EM$sigma[1]) + (w==2)*rnorm(50,EM$mu[2],
                                                              EM$sigma[2])
  EM <- density(X, n = 100)
  plot(EM, col = 4, lwd = 2, type = "l",main = " ", xlab=paste("alpha = ", i),
       xlim = c(-1,9), ylim = c(0,0.8))
  lines(True, col = 1, lwd = 2, type = "h")
  lines(x, fnorm, col=2, lwd = 2)
  legend(-0.4, 0.8, legend=c("Gaussian kernel", "EM estimation", "True density"),
         col=c(2, 4, "black"),lty=c(1,1,5), lwd = 2, cex = 1)
  dev.off()
}


#########################################
# Bias for parametric and nonparametric #
#########################################

X <- mixture(0.05) ## Chaege alpah to perform for different mixtures
x <- seq(from = min(X), to = max(X), length = 100)
True <- density(X, n = 100)
True
prDensEst <- function(x, X, h, K) mean(K((x - X)/h))/h
u2Knorm <- function(u) u^2 * dnorm(u)
h_norm <- h/IntKnorm^0.5
Knorm <- function(u) dnorm(u)
fnorm <- sapply(x, function(x) prDensEst(x, X, h_norm, Knorm)) 

bias<-fnorm- True$y
mean(bias)
mse <- mean((fnorm- True$y)^2)
mse
## Task 2

## local bw

for (i in c(0.25, 0.75)) {
  X <- mixture(i)
  h<- bw.nrd(X)
  print(h)
}

################################################
# Monte Carlo Simulation implemetation for MSE #
################################################

rm(list = ls())
K <- 1E4

mixture <- function(alpha, n) {
  w1 <- alpha
  w2 <- 1 - alpha
  w <- sample(c(1,2), size = n, replace = TRUE, prob = c(w1,w2))
  X <- (w==1)*rnorm(50, 2, 1) + (w==2)*rnorm(50,6,1)
  return(X)
}

n.seq <- c(25, 50, 75, 100, 200, 500, 1000)
prDensEst <- function(x, X, h, K) mean(K((x - X)/h))/h ## Gaussian
Knorm <- function(u) dnorm(u)

Result <- function(n) {
  
  H <- matrix(nrow = K, ncol = 1) ## to stack the bw
  EM <- matrix(nrow = K, ncol = 4)
  SE <- matrix(nrow=K, ncol=8)
  Kern <- matrix(nrow = K, ncol = 1)
  
  for (k in 1:K) {
    set.seed(k)
    X <- mixture(0.25, n)
    
    h <- bw.nrd(X) ## bw normal rule of thumb
    
    par <- normalmixEM(X, lambda = 0.25, mu = c(2, 6), sigma = c(1,1))
    ## EM parametric
    w1 <- par$lambda[1]
    w2 <- par$lambda[2]
    w <- sample(c(1,2), size = n, replace = TRUE, prob = c(w1,w2))
    Xpar <- (w==1)*rnorm(n/2, par$mu[1], par$sigma[1]) + 
      (w==2)*rnorm(n/2,par$mu[2],par$sigma[2])
    denspar <- density(Xpar, n = n)
    
    x <- seq(from = min(X), to = max(X), length = n) ## nonparam gaussian kernel
    fhat<- sapply(x, function(x) prDensEst(x, X, h, Knorm))
    
    SE[k, ] <- c((fhat[1] - density(X, n = n)$y[1])^2, 
                 (fhat[round(n/4)] - density(X, n = n)$y[round(n/4)])^2, 
                 (fhat[round((3*n)/4)] - density(X, n = n)$y[round((3*n)/4)])^2, 
                 (fhat[n] - density(X, n = n)$y[n])^2,
                 (denspar$y[1] - density(X, n = n)$y[1])^2, 
                 (denspar$y[round(n/4)] - density(X, n = n)$y[round(n/4)])^2, 
                 (denspar$y[round((3*n)/4)] - 
                    density(X, n = n)$y[round((3*n)/4)])^2, 
                 (denspar$y[n] - density(X, n = n)$y[n])^2)
    H[k, ] <- h
    EM[k, ] <- c(par$mu[1], par$mu[2], par$sigma[1], par$sigma[2])
    Kern[k, ] <- fhat[5]
  }
  
  return(c(colMeans(SE), colMeans(H), colMeans(EM), colMeans(Kern)))
  
}

MSE <- sapply(n.seq, function(n.seq) Result(n.seq))
table <- MSE
rownames(table) <- c("K MSE at x = 1", "K MSE at x = n/4", 
                     "K MSE at x = 3n/4", "K MSE at x = n",
                     "EM MSE at x = 1", "EM MSE at x = n/4", 
                     "EM MSE at x = 3n/4", "EM MSE at x = n", 
                     "Estimated SR bw", "mu1", "mu2", "sigma1", 
                     "sigma2", "kern")
colnames(table) <- c("n = 25", "n = 50", "n = 75", "n = 100", "n = 200", 
                     "n = 500", "n = 1000")

kbl(table, caption = "Monte Carlo simulation for 0.75") %>% kable_material()
MSE

par(mfrow = c(1, 1))

png(filename = paste("MSEK0.75.png"), width = 7, height = 6,
    units = 'in', res=600)
plot(n.seq, MSE[1,], type="l", main="MSE of kernel density estimator", 
     xlab = "sample size",
     ylab = "MSE", ylim = c(0, 9.350343e-03))
lines(n.seq, MSE[2,], col = 2)
lines(n.seq, MSE[3,], col = 3)
lines(n.seq, MSE[4,], col = 4)
abline(h = 0, lty = 2)
legend("topright", c("MSE at x = 1", "MSE at x = n/4", "MSE at x = 3n/4", 
                     "MSE at x = n"),
       col = c(1,2,3,4), lty = c(1,1,1,1), bty = "n")
dev.off()

png(filename = paste("MSEEM0.75.png"), width = 7, height = 6, units = 'in',
    res=600)
plot(n.seq, MSE[5,], type="l", main="MSE of EM density estimator", 
     xlab = "sample size", 
     ylab = "MSE", ylim = c(0, 1.436243e-02))
lines(n.seq, MSE[6,], col = 2)
lines(n.seq, MSE[7,], col = 3)
lines(n.seq, MSE[8,], col = 4)
abline(h = 0, lty = 2)
legend("topright", c("MSE at x = 1", "MSE at x = n/4",
                     "MSE at x = 3n/4", "MSE at x = n"),
       col = c(1,2,3,4), lty = c(1,1,1,1), bty = "n")
dev.off()




#################################################
# Monte Carlo Simulation implemetation for Bias #
#################################################


rm(list = ls())
K <- 1E4

mixture <- function(alpha, n) {
  w1 <- alpha
  w2 <- 1 - alpha
  w <- sample(c(1,2), size = n, replace = TRUE, prob = c(w1,w2))
  X <- (w==1)*rnorm(50, 2, 1) + (w==2)*rnorm(50,6,1)
  return(X)
}

n.seq <- c(25, 50, 75, 100, 200, 500, 1000)
prDensEst <- function(x, X, h, K) mean(K((x - X)/h))/h ## Gaussian
Knorm <- function(u) dnorm(u)

Result <- function(n) {
  
  B <- matrix(nrow=K, ncol=8)
  
  for (k in 1:K) {
    set.seed(k)
    X <- mixture(0.25, n)
    
    h <- bw.nrd(X) ## bw normal rule of thumb
    
    par <- normalmixEM(X, lambda = 0.25, mu = c(2, 6), sigma = c(1,1))
    ## EM parametric
    w1 <- par$lambda[1]
    w2 <- par$lambda[2]
    w <- sample(c(1,2), size = n, replace = TRUE, prob = c(w1,w2))
    Xpar <- (w==1)*rnorm(n/2, par$mu[1], par$sigma[1]) +
      (w==2)*rnorm(n/2,par$mu[2],par$sigma[2])
    denspar <- density(Xpar, n = n)
    
    x <- seq(from = min(X), to = max(X), length = n) ## nonparam gaussian kernel
    fhat<- sapply(x, function(x) prDensEst(x, X, h, Knorm))
    
    B[k, ] <- c((fhat[1] - density(X, n = n)$y[1]), 
                (fhat[round(n/4)] - density(X, n = n)$y[round(n/4)]), 
                 (fhat[round((3*n)/4)] - density(X, n = n)$y[round((3*n)/4)]),
                (fhat[n] - density(X, n = n)$y[n]),
                 (denspar$y[1] - density(X, n = n)$y[1]), 
                (denspar$y[round(n/4)] - density(X, n = n)$y[round(n/4)]), 
                 (denspar$y[round((3*n)/4)] - 
                    density(X, n = n)$y[round((3*n)/4)]), 
                (denspar$y[n] - density(X, n = n)$y[n]))
  }
  
  return(colMeans(B))
  
}

Bias <- sapply(n.seq, function(n.seq) Result(n.seq))
Bias
table <- Bias
rownames(table) <- c("K Bias at x = 1", "K Bias at x = n/4", 
                     "K Bias at x = 3n/4", "K Bias at x = n", 
                     "EM Bias at x = 1", "EM Bias at x = n/4",
                     "EM Bias at x = 3n/4", "EM Bias at x = n")
colnames(table) <- c("n = 25", "n = 50", "n = 75", "n = 100", 
                     "n = 200", "n = 500", "n = 1000")

kbl(table, caption = "Monte Carlo simulation for 0.75") %>% kable_material()

par(mfrow = c(1, 1))

png(filename = paste("BiasK0.75.png"), width = 7, height = 6, units = 'in', res=600)
plot(n.seq, Bias[1,], type="l", main="Bias of kernel density estimator", 
     xlab = "sample size",
     ylab = "Bias", ylim = c(0, 9.350343e-02))
lines(n.seq, Bias[2,], col = 2)
lines(n.seq, Bias[3,], col = 3)
lines(n.seq, Bias[4,], col = 4)
abline(h = 0, lty = 2)
legend("topright", c("Bias at x = 1", "Bias at x = n/4", "Bias at x = 3n/4",
                     "Bias at x = n"),
       col = c(1,2,3,4), lty = c(1,1,1,1), bty = "n")
dev.off()

png(filename = paste("BiasEM0.75.png"), width = 7, height = 6, units = 'in',
    res=600)
plot(n.seq, Bias[5,], type="l", main="Bias of EM density estimator", 
     xlab = "sample size",
     ylab = "Bias", ylim = c(-4.886804e-02, 0.0388959))
lines(n.seq, Bias[6,], col = 2)
lines(n.seq, Bias[7,], col = 3)
lines(n.seq, Bias[8,], col = 4)
abline(h = 0, lty = 2)
legend("topright", c("Bias at x = 1", "Bias at x = n/4", "Bias at x = 3n/4",
                     "Bias at x = n"),
       col = c(1,2,3,4), lty = c(1,1,1,1), bty = "n")
dev.off()

##########################################
# Analytical computation of the variance #
##########################################

VAR <- MSE[c(1,2,3,4,5,6,7,8),] - (Bias[c(1,2,3,4,5,6,7,8),])^2
VAR
table <- VAR
rownames(table) <- c("K Var at x = 1", "K Var at x = n/4",
                     "K Var at x = 3n/4", "K Var at x = n", 
                     "EM Var at x = 1", "EM Var at x = n/4", 
                     "EM Var at x = 3n/4", "EM Var at x = n")
colnames(table) <- c("n = 25", "n = 50", "n = 75", "n = 100", "n = 200", 
                     "n = 500", "n = 1000")
kbl(table, caption = "Monte Carlo simulation for 0.25") %>% kable_material()


##################################
# Display Bias variance trade-off #
##################################

png(filename = paste("VarK0.25.png"), width = 7, height = 6, units = 'in', 
    res=600)
plot(n.seq, VAR[1,], type="l", main="Variance of kernel density estimator",
     xlab = "sample size", 
     ylab = "Var", ylim = c(0, 0.0027385))
lines(n.seq, VAR[2,], col = 2)
lines(n.seq, VAR[3,], col = 3)
lines(n.seq, VAR[4,], col = 4)
abline(h = 0, lty = 2)
legend("topright", c("Var at x = 1", "Var at x = n/4", "Var at x = 3n/4",
                     "Var at x = n"),
       col = c(1,2,3,4), lty = c(1,1,1,1), bty = "n")
dev.off()

png(filename = paste("VarEM0.25.png"), width = 7, height = 6, units = 'in',
    res=600)
plot(n.seq, VAR[5,], type="l", main="Variance of EM density estimator", 
     xlab = "sample size", 
     ylab = "Var", ylim = c(0, 0.0128495))
lines(n.seq, VAR[6,], col = 2)
lines(n.seq, VAR[7,], col = 3)
lines(n.seq, VAR[8,], col = 4)
abline(h = 0, lty = 2)
legend("topright", c("Var at x = 1", "Var at x = n/4", "Var at x = 3n/4",
                     "Var at x = n"),
       col = c(1,2,3,4), lty = c(1,1,1,1), bty = "n")
dev.off()

png(filename = paste("traoff0.251.png"), width = 7, height = 6, units = 'in', 
    res=600)
hesti <- c(0.8434, 0.8611, 0.8030, 0.7643, 0.6732, 0.5673, 0.4972) * 10
x1 <- as.data.frame(cbind(Bias[3,]^2 * 10^3, MSE[3,] * 10^3, VAR[3,] * 10^3,
                          hesti))
colnames(x1) <- c("Bias", "MSE", "VAR", "h")
x1
ggplot(x1, aes(x = h, y = Bias)) + 
  geom_line(color = "red") + geom_line(aes(y=MSE), color = "blue") +
  geom_line(aes(y=VAR), color = "green") + geom_vline(xintercept = 5.84) + 
  labs(title = 
         TeX("Pointwise Bias-Variance trade-off $\\alpha = 0.25$ at 
             $\\x = \\frac{3n}{4}$"), 
       x = "Bandwidth h*E-1", y = "Bias^2, MSE and Variance * E-3") + 
  annotate(geom="text", x=8, y=2.1, label="VAR",color="green", size = 7) +
  annotate(geom="text", x=7.6, y=4.4, label="BIAS",color="red", size = 7) +
  annotate(geom="text", x=7.8, y=6.4, label="MSE",color="blue", size = 7)
dev.off()

png(filename = paste("traoff0.252.png"), width = 7, height = 6, units = 'in', 
    res=600)
hesti <- c(0.8434, 0.8611, 0.8030, 0.7643, 0.6732, 0.5673, 0.4972) * 10
x1 <- as.data.frame(cbind(Bias[2,]^2 * 10^5, MSE[2,] * 10^5, VAR[2,] * 10^5,
                          hesti))
colnames(x1) <- c("Bias", "MSE", "VAR", "h")
x1
ggplot(x1, aes(x = h, y = Bias)) + geom_line(color = "red") +
  geom_line(aes(y=MSE), color = "blue") +
  geom_line(aes(y=VAR), color = "green") + 
  labs(title = TeX("Pointwise Bias-Variance trade-off $\\alpha = 0.25$ at 
                   $\\x = \\frac{n}{4}$"), 
       x = "Bandwidth h*E-1", y = "Bias^2, MSE and Variance * E-5") + 
  annotate(geom="text", x=7.6, y=50, label="VAR",color="green", size = 7) +
  annotate(geom="text", x=8, y=16, label="BIAS",color="red", size = 7) +
  annotate(geom="text", x=7.6, y=77, label="MSE",color="blue", size = 7)
dev.off()







