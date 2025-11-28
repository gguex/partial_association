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
r_w_vec = c(0, 0.5, 0.9)
r_c_mat = matrix(c(0, 0.1, 
                   0, 0.35,
                   0, 0.6), nrow=3, ncol=2, byrow=TRUE)

n_r_b = 10

# -------------------------------------------------
# Code
# -------------------------------------------------

# Function to compute 
c_results = function(n, p1, p2, p3, R_mat){
  
  # Create the weights 
  f = runif(n)
  f = f/sum(f)
  
  # The diagonal matrix of weights and the sqrt
  Pi_sqrt = diag(sqrt(f))
  
  # The centering matrix
  H = diag(n) - outer(rep(1, n), f)
  
  # Generate the dataset
  data_sim = diag(1/sqrt(f))%*%mvrnorm(n, mu=rep(0, p1+p2+p3), Sigma=R_mat)
  
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

for(i in 1:length(r_w_vec)){
  
  r_w = r_w_vec[i]
  df_all_res = data.frame()
  
  ### Computations
  r_c_vec = r_c_mat[i, ]
  
  for(r_c in r_c_vec){
    
    r_b_max = (1 + (p1 - 1)*r_w)/p1
    r_b_vec = seq(0, r_b_max, length.out=n_r_b)
    
    for(r_b in r_b_vec){
      
      cat("Running for r_w =", r_w, "; r_c =", r_c,"; r_b =", r_b, "\n")
      
      # Create the covariance matrix
      R_res = gen_cor_3datasets(p1, p2, p3, r_w, r_w, r_w, r_b, r_c, r_c)
      R_mat = R_res$R_mat
      semidef_pos = R_res$semidef_pos
      
      res = mclapply(1:n_test, 
                     function(x) c_results(n, p1, p2, p3, R_mat), 
                     mc.cores=n_cores)
      
      df_res = as.data.frame(apply(t(simplify2array(res)), 2, unlist))
      df_res$n = n
      df_res$p1 = p1
      df_res$p2 = p2
      df_res$p3 = p3
      df_res$r1 = r_w
      df_res$r2 = r_w
      df_res$r3 = r_w
      df_res$r12 = r_b
      df_res$r13 = r_c
      df_res$r23 = r_c
      df_res$semidef_pos = 1*semidef_pos
      
      df_all_res = rbind(df_all_res, df_res)
      
    }
  }
  
  # Save the results
  write.csv(df_all_res, paste0("results_csv/new_3k_rw_", r_w, ".csv"), 
            row.names=F)
}