#################################################################################
# R code to estimate Markov chain dynamics with RCS data
# Input -> See Table 1 on Okabe and Nogiwa (20XX) for input data structure.
# y: N-vector of PID observations 
# x: (N x Nx) matrix of individual characteristic observations
# Ps_init: 3-vector of initical PID distribution
#
# Author: Tomohito Okabe, Hitotsubashi University, email: t-okabe@ier.hit-u.ac.jp
# For details, see Okabe, T. and Nogiwa D.,
# ``Estimation of Unobserved Dynamics of Individual Partisanship: 
# A Bayesian Approach", XXXXX, XX (20XX): - .
#################################################################################

library(ggplot2)
library(rstan)
library(reshape2)
library(plyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # Parallel computation 
set.seed(123)

# Input data
y <- scan("PID_sample.txt") # MUST BE CHANGED. PID responses, N-vector  
N <- length(y)

period <- scan("PERIOD_sample.txt") # Time period index, N-vector

N_x <- 5 # MUST BE CHANGED. Place the number of explanatory variables (individual characteristics)
x <- matrix(scan("CHARACTERISTICS_sample.txt"), N, Nx) # MUST BE CHANGED. Individual characteristics, 
													   # N x Nx matrix
x <- cbind(rep(1,N), x) # Add ones for constant terms
D <- ncol(x)

Ps_init <- c(0.4,0.2,0.4) # MUST BE CHANGED. Place the initial distribution.

K <- 3 # Total number of PID alternatives

# Set all input as a set of data
data <- list(N=N, K=K, D=D, y=y, period=period, x=x, Ps_init=Ps_init )

# Stan Computation
# MUST BE CHANGED. Set computation parameters
fit <- stan(file='MCMC_logit.stan', data=data, iter=10000, warmup=5000, thin=2, chains=4) 

# Save result
save.image("result.rdata")

# Show result 
print(fit)
# summary(fit)$summary


# Plot for convergence
# traceplot(object=fit, pars="phi") 


