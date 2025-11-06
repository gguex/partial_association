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

# -------------------------------------------------
# Parameters for the experiment
# -------------------------------------------------

set.seed(1234)
n = 100
p1 = 5
p2 = 5
n_test = 1000
r_w_vec = c(0, 0.5, 0.9)
r_b_vec = seq(0, 1, 0.1)
n_cores = detectCores() - 2

# -------------------------------------------------
# Code
# -------------------------------------------------

# Function to compute 
c_results = function(n, R_mat){
  
  # Create the weights 
  f = runif(n)
  f = f/sum(f)
  
  # The diagonal matrix of weights and the sqrt
  Pi_sqrt = diag(sqrt(f))
  
  # The centering matrix
  H = diag(n) - outer(rep(1, n), f)
  
  # Generate the dataset
  data_sim = mvrnorm(n, mu=rep(0, p1+p2), Sigma=R_mat)
  
  # Split the dataset
  X = data_sim[, 1:p1]
  Y = data_sim[, (p1+1):(p1+p2)]
  
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

for(r_b in r_b_vec){
  cat("Running for r_b =", r_b, "\n")
  
  df_all_res = data.frame()
  
  ### Computations
  
  for(r_w in r_w_vec){
    
    cat("Running for r_w =", r_w, "\n")
    
    # Create the covariance matrix
    R_res = gen_cor_2datasets(p1, p2, r_w, r_w, r_b)
    R_mat = R_res$R_mat
    semipos_def = R_res$semidef
    
    res = mclapply(1:n_test, 
                   function(x) c_results(n, R_mat), 
                   mc.cores=n_cores)
    
    df_res = as.data.frame(apply(t(simplify2array(res)), 2, unlist))
    df_res$n = n
    df_res$p1 = p1
    df_res$p2 = p2
    df_res$r1 = r_w
    df_res$r2 = r_w
    df_res$r12 = r_b
    df_res$semipos_def = 1*semipos_def
    
    df_all_res = rbind(df_all_res, df_res)
    
  }
  
  # Save the results
  write.csv(df_all_res, paste0("results_csv/res_2kernels_rb_", r_b, ".csv"), 
            row.names=F)
}

