#question 2
install.packages("maxLik")
library(maxLik)
# Load the data
load("dataex2.Rdata")
data <- dataex2
# Define the negative log-likelihood function
neg_log_likelihood <- function(mu, x, r, sigma2) {
  n <- length(x)
  phi <- dnorm(x, mean = mu, sd = sqrt(sigma2))
  Phi <- pnorm(x, mean = mu, sd = sqrt(sigma2))
  
  logL <- sum(r * log(phi) + (1 - r) * log(Phi))
  return(-logL)
}

# Known sigma^2
sigma2_known <- 1.5^2

# Initial guess for mu
mu_init <- mean(data$X, na.rm = TRUE)

# Find the maximum likelihood estimate of mu
mle_result <- optim(
  mu_init, 
  neg_log_likelihood, 
  x = data$X, 
  r = data$R, 
  sigma2 = sigma2_known
)

# Display the MLE of mu
mle_result$par




#question 4
# Load data
load("dataex4.Rdata")

# Define log_fun function
log_fun <- function(beta, data = dataex4) {
  beta_0 <- beta[1]
  beta_1 <- beta[2]
  
  y_obs <- which(!is.na(data$Y))
  y_mis <- which(is.na(data$Y))
  n <- rep(0, length(data$Y)) 
  
  # Compute observed data
  for (i in y_obs) {
    n[i] <- data$Y[i] * (beta_0 + beta_1 * data$X[i]) - 
      log(1 + exp(beta_0 + beta_1 * data$X[i]))
  }
  
  return(sum(n)) 
}

# Set initial values
initial <- c(0.5, 0.7) 

# Estimate parameters using maxLik
estimates <- maxLik(log_fun, 
                    start = initial,
                    method = 'BFGS', 
                    data = dataex4) 
estimates


#question 5
load("dataex5.Rdata")
em.mixture.two <- function(y, theta0, eps){ 
  n <- length(y)
  theta <- theta0
  p <- theta[1]
  lambda <- theta[2]
  mu <- theta[3] 
  
  diff <- 1
  while(diff > eps){
    
    theta.old <- theta
    #E-step 
    ptilde1 <- p*lambda*y^(-lambda-1)
    ptilde2 <- (1-p)*mu*y^(-mu-1)
    ptilde <- ptilde1/(ptilde2+ptilde1)
    
    #M-step
    p <- mean(ptilde)
    lambda <- sum(ptilde)/sum(ptilde*log(y))
    mu <- sum(1 - ptilde)/sum((1 - ptilde)*log(y))
    
    theta <- c(p, lambda, mu)
    diff <- sum(abs(theta - theta.old))
  }
  
  return(theta)
}
theta0 <- c(0.3, 0.3, 0.4)
eps <- 0.0001
res <- em.mixture.two(y = dataex5, theta0, eps)

p <- res[1]
lambda <- res[2]
mu <- res[3]
p; lambda; mu

#calculate the interquartile range (IQR)
q1 <- quantile(dataex5, 0.25)
q3 <- quantile(dataex5, 0.75)
IQR <- q3 - q1
#calculate the range for x-axis and y-axis
den <- density(dataex5)
y_range <- c(0, max(den$y))
x_range <- c(min(dataex5), q3 + 1.5*IQR)
# Create the histogram with the specified x-axis range
hist(dataex5,
     freq=FALSE, 
     breaks = "FD", 
     xlim = x_range,ylim=y_range, 
     xlab = "Vaule", 
     ylab = "Density", 
     main = "Histogram of data with estimate density")
#add curve of the density function
curve(p*(lambda*x^(-lambda-1)) 
      + (1 - p)*(mu*x^(-mu-1)),
      add = TRUE,
      lwd = 2, 
      col = "red")