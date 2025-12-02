# ----------------------------------------------------
# ----------------------------------------------------
# Compute partial Cov and moments for 3 sets of variables
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
p3 = 5
n_test = 1000

sd_E = 1
m_A = 0
m_B = 0

sd_B_vec = c(0, 0.25, 0.5)
sd_A_vec = seq(0, 1, length.out=11)

# -------------------------------------------------
# Code
# -------------------------------------------------

# Function to compute 
c_results = function(n, p1, p2, p3, m_A, sd_A, m_B, sd_B, sd_E){
  
  # Create the weights 
  f = runif(n)
  f = f/sum(f)
  
  # The diagonal matrix of weights and the sqrt
  Pi_sqrt = diag(sqrt(f))
  
  # The centering matrix
  H = diag(n) - outer(rep(1, n), f)
  
  # Generate the dataset
  data_mat = gen_3datasets(n, f, p1, p2, p3, m_A, sd_A, m_B, sd_B, sd_E)
  
  # Split the dataset
  X = data_mat[, 1:p1]
  Y = data_mat[, (p1+1):(p1+p2)]
  Z = data_mat[, (p1+p2+1):(p1+p2+p3)]
  
  # Compute the euclidean squared dissimilarity matrices 
  D_X = as.matrix(dist(X))^2
  D_Y = as.matrix(dist(Y))^2
  D_Z = as.matrix(dist(Z))^2
  
  # Compute the kernels
  K_X = -0.5 * Pi_sqrt %*% H %*% D_X %*% t(H) %*% Pi_sqrt
  K_Y = -0.5 * Pi_sqrt %*% H %*% D_Y %*% t(H) %*% Pi_sqrt
  K_Z = -0.5 * Pi_sqrt %*% H %*% D_Z %*% t(H) %*% Pi_sqrt
  
  # Compute the cross-covariance and theoretical moments
  C_res2 = compute_C(K_X, K_Y)
  C_res3 = compute_partial_C(K_X, K_Y, K_Z)
  
  return(c(C_res2, C_res3))
}

# Create empty dataset
df_all_res = data.frame()

for(sd_B in sd_B_vec){
  for(sd_A in sd_A_vec){
    
    cat("Running for m_A =", m_A, "; sd_A =", sd_A, "m_B =", m_B, "; sd_B =", 
        sd_B, "sd_E =", sd_E, "\n")
      
    res = mclapply(1:n_test, 
                   function(x) c_results(n, p1, p2, p3, 
                                         m_A, sd_A, m_B, sd_B, sd_E), 
                   mc.cores=n_cores)
      
    df_res = as.data.frame(apply(t(simplify2array(res)), 2, unlist))
    df_res$n = n
    df_res$p1 = p1
    df_res$p2 = p2
    df_res$p3 = p3
    df_res$m_A = m_A
    df_res$sd_A = sd_A
    df_res$m_B = m_B
    df_res$sd_B = sd_B
    df_res$sd_E = sd_E
      
    df_all_res = rbind(df_all_res, df_res)
      
  }
}
  
# Save the results
write.csv(df_all_res, paste0("results_csv/nexp_3k.csv"), 
          row.names=F)