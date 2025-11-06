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

# -------------------------------------------------
# Parameters for the experiment
# -------------------------------------------------

set.seed(1234)

n = 500
p1 = 5
p2 = 5
p3 = 5
n_test_out = 10
n_test_in = 10
r_w_vec = seq(0, 1, 0.1)
r_b1_vec = seq(0, 1, 0.1)
r_b2_vec = seq(0, 1, 0.1)
n_cores = detectCores() - 2

# -------------------------------------------------
# Code
# -------------------------------------------------

# Function to compute 
c_results = function(n, p1, p2, p3, Sigma, Pi_sqrt, H){
  # Generate the dataset
  data_sim = mvrnorm(n = n, mu = rep(0, p1+p2+p3), Sigma = Sigma)
  
  # Split the dataset
  X = data_sim[, 1:p1]
  Y = data_sim[, (p1+1):(p1+p2)]
  Z = data_sim[, (p1+p2+1):(p1+p2+p3)]
  
  # Compute the euclidean squared dissimilarity matrices 
  D_X = as.matrix(dist(X))^2
  D_Y = as.matrix(dist(Y))^2
  D_Z = as.matrix(dist(Z))^2
  
  # Compute the kernels
  K_X = -0.5 * Pi_sqrt %*% H %*% D_X %*% t(H) %*% Pi_sqrt
  K_Y = -0.5 * Pi_sqrt %*% H %*% D_Y %*% t(H) %*% Pi_sqrt
  K_Z = -0.5 * Pi_sqrt %*% H %*% D_Z %*% t(H) %*% Pi_sqrt
  
  # Compute the cross-covariance and theoretical moments
  C_res = compute_partial_C(K_X, K_Y, K_Z)
  
  return(C_res)
}

for(r_w in r_w_vec){
  cat("Running for rw =", r_w, "\n")
  r1 = r_w
  r2 = r_w
  r3 = r_w
  
  df_all_res = data.frame()
  
  ### Computations
  
  for(r_b1 in r_b1_vec){
    for(r_b2 in r_b2_vec){
      cat("Running for r_b1 =", r_b1,", r_b2 =", r_b2, "\n")
      for(i in 1:n_test_out){
        # Create the weights 
        f = runif(n)
        f = f/sum(f)
        
        # The diagonal matrix of weights and the sqrt
        Pi = diag(f)
        Pi_sqrt = diag(sqrt(f))
        
        # The centering matrix
        H = diag(n) - outer(rep(1, n), f)
        
        # Create the covariance matrix
        R_mat = gen_cor_3datasets(p1, p2, p3, r1, r2, r3, r_b1, r_b2, r_b2)
        
        res = mclapply(1:n_test_in, 
                       function(x) c_results(n, p1, p2, p3, R_mat, Pi_sqrt, H), 
                       mc.cores=n_cores)
        
        df_res = as.data.frame(apply(t(simplify2array(res)), 2, unlist))
        df_res$n = n
        df_res$p1 = p1
        df_res$p2 = p2
        df_res$p3 = p3
        df_res$r1 = r1
        df_res$r2 = r2
        df_res$r3 = r3
        df_res$r12 = r_b1
        df_res$r13 = r_b2
        df_res$r23 = r_b2
        
        df_all_res = rbind(df_all_res, df_res)
      }
    }
  }
  
  # Save the results
  write.csv(df_all_res, paste0("results_csv/res_3kernels_rw_", r_w, ".csv"), 
            row.names=F)
}