# ----------------------------------------------------
# ----------------------------------------------------
# Functions for the projects
# ----------------------------------------------------
# ----------------------------------------------------

# --------------------------------------------
# Generate the covariance matrix 
# for 2 sets of variables 
# --------------------------------------------

generate_cov2dataset = function(p1, p2, r1, r2, r12, 
                                s1_min=0, s1_max=1, s2_min=0, s2_max=1){
  
  # Create correlation matrix
  R1 = matrix(r1, nrow = p1, ncol = p1)
  diag(R1) = 1
  R2 = matrix(r2, nrow = p2, ncol = p2)
  diag(R2) = 1
  R12 = matrix(r12, nrow = p1, ncol = p2)
  
  R = rbind( cbind(R1, R12), cbind(t(R12), R2) )
  
  # Generate variances
  s1 = runif(p1, s1_min, s1_max)
  s2 = runif(p2, s2_min, s2_max)
  
  # Create the covariance matrix and return it
  Cov_m = R * outer(c(s1, s2), c(s1, s2))
  
  # Make sure its is positive definite
  eigs = eigen(Cov_m)
  if(any(eigs$values < 0)){
    cat("Warning: covariance matrix is not positive definite. \n")
  }
  eigs$values[eigs$values < 0] = 0
  Cov_m = eigs$vectors %*% diag(eigs$values) %*% t(eigs$vectors)
  
  return(Cov_m)
}

# -------------------------------------------------
# Compute the cross-covariance matrix C_XY and
# theoretical moments for two kernels
# -------------------------------------------------

compute_CV = function(K_X, K_Y) {
  
  # Get n
  n = nrow(K_X)
  
  # Compute C_XY
  C_XY = sum(diag(K_X %*% K_Y))
  
  # Compute the z score
  tr_K_X = sum(diag(K_X))
  tr_K_Y = sum(diag(K_Y))
  E_C_XY = tr_K_X * tr_K_Y / (n-1)
  delta2_X = sum(diag(K_X%*%K_X)) 
  delta2_Y = sum(diag(K_Y%*%K_Y))
  v_X = tr_K_X^2 / delta2_X 
  v_Y = tr_K_Y^2 / delta2_Y
  Var_C_XY = 2 * delta2_X * delta2_Y * (n - 1 - v_X) * (n - 1 - v_Y) / 
    ( (n-2)*(n-1)^2*(n+1) )
  
  return(list(C_XY = C_XY, 
              E_C_XY = E_C_XY, 
              Var_C_XY = Var_C_XY))
}

