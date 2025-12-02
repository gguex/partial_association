# ----------------------------------------------------
# ----------------------------------------------------
# Compute Cov and moments for 2 sets of variables
# ----------------------------------------------------
# ----------------------------------------------------

# -------------------------------------------------
# Source and libraries
# -------------------------------------------------

source("local_functions.R")
library(MASS)
library(parallel)

n_cores = detectCores() - 2

# -------------------------------------------------
# Parameters for the experiment
# -------------------------------------------------

set.seed(1234)
n = 100
p1 = 5
p2 = 5
n_test = 1000

sd_E = 1
m_A = 0

sd_A_vec = seq(0, 1, length.out=11)

# -------------------------------------------------
# Code
# -------------------------------------------------

# Function to compute 
c_results = function(n, p1, p2, m_A, sd_A, sd_E){
  
  # Create the weights 
  f = runif(n)
  f = f/sum(f)
  
  # The diagonal matrix of weights and the sqrt
  Pi_sqrt = diag(sqrt(f))
  
  # The centering matrix
  H = diag(n) - outer(rep(1, n), f)
  
  # Generate dataset
  data_mat = gen_2datasets(n, f, p1, p2, m_A, sd_A, sd_E)
  
  # Split the dataset
  X = data_mat[, 1:p1]
  Y = data_mat[, (p1+1):(p1+p2)]
  
  # Compute the euclidean squared dissimilarity matrices 
  D_X = as.matrix(dist(X))^2
  D_Y = as.matrix(dist(Y))^2
  
  # Compute the kernels
  K_X = -0.5 * Pi_sqrt %*% H %*% D_X %*% t(H) %*% Pi_sqrt
  K_Y = -0.5 * Pi_sqrt %*% H %*% D_Y %*% t(H) %*% Pi_sqrt
  
  # Compute the cross-covariance and theoretical moments
  C_res = compute_C(K_X, K_Y)
  
  return(C_res)
}


# Create empty dataset
df_all_res = data.frame()

for(sd_A in sd_A_vec){
  
  cat("Running for m_A =", m_A, "; sd_A =", sd_A, "sd_E =", sd_E, "\n")

  # Compute the results in parallel
  res = mclapply(1:n_test, 
                 function(x) c_results(n, p1, p2, m_A, sd_A, sd_E), 
                 mc.cores=n_cores)
  
  df_res = as.data.frame(apply(t(simplify2array(res)), 2, unlist))
  df_res$n = n
  df_res$p1 = p1
  df_res$p2 = p2
  df_res$m_A = m_A
  df_res$sd_A = sd_A
  df_res$sd_E = sd_E
  
  df_all_res = rbind(df_all_res, df_res)
}

# Save the results
write.csv(df_all_res, paste0("results_csv/nexp_2k.csv"), 
          row.names=F)

