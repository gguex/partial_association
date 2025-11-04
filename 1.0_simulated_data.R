# ----------------------------------------------------
# ----------------------------------------------------
# Test the simulated data generation process
# ----------------------------------------------------
# ----------------------------------------------------

# -------------------------------------------------
# Source and libraries
# -------------------------------------------------

source("local_functions.R")
library(MASS)

# -------------------------------------------------
# Parameters for the experiment
# -------------------------------------------------

set.seed(1234)
n = 500
p1 = 5
p2 = 5
r1 = 0.6
r2 = 0.6
r12 = 0.1

# -------------------------------------------------
# Code
# -------------------------------------------------

### Create the fixed quantities

# Create the weights 
f = runif(n)
f = f/sum(f)

# The diagonal matrix of weights
Pi = diag(f)

# The centering matrix
H = diag(n) - outer(rep(1, n), f)

### Computations

# Create the covariance matrix
Cov_mat = generate_cov2dataset(p1, p2, r1, r2, r12)

# Generate the dataset
data_sim = mvrnorm(n = n, mu = rep(0, 10), Sigma = Cov_mat)

# Split the dataset
X = data_sim[, 1:p1]
Y = data_sim[, (p1+1):(p1+p2)]


# Compute the kernels 
X_s = sqrt(f)*t(t(X) - colSums(f*X))
K_X = 0.5 * X_s %*% t(X_s)
Y_s = sqrt(f)*t(t(Y) - colSums(f*Y))
K_Y = 0.5 * Y_s %*% t(Y_s)







