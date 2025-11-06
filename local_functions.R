# ----------------------------------------------------
# ----------------------------------------------------
# Functions for the projects
# ----------------------------------------------------
# ----------------------------------------------------

# --------------------------------------------
# Generate the correlation matrix 
# for 2 sets of variables 
# --------------------------------------------

gen_cor_2datasets = function(p1, p2, r1, r2, r12){
  
  # Create correlation matrix
  R1 = matrix(r1, nrow = p1, ncol = p1)
  diag(R1) = 1
  R2 = matrix(r2, nrow = p2, ncol = p2)
  diag(R2) = 1
  R12 = matrix(r12, nrow = p1, ncol = p2)
  
  # The correlation matrix
  R = rbind( cbind(R1, R12), 
             cbind(t(R12), R2) )
  
  # Make sure its is positive definite
  eigs = eigen(R)
  if(any(eigs$values < 0)){
    cat("Warning: the correlation matrix is not positive definite. \n")
  }
  eigs$values[eigs$values < 0] = 0
  R = eigs$vectors %*% diag(eigs$values) %*% t(eigs$vectors)
  
  return(R)
}

# --------------------------------------------
# Generate the correlation matrix 
# for 3 sets of variables 
# --------------------------------------------

gen_cor_3datasets = function(p1, p2, p3, r1, r2, r3, r12, r13, r23){
  
  # Create correlation matrix
  R1 = matrix(r1, nrow = p1, ncol = p1)
  diag(R1) = 1
  R2 = matrix(r2, nrow = p2, ncol = p2)
  diag(R2) = 1
  R3 = matrix(r3, nrow = p3, ncol = p3)
  diag(R3) = 1
  R12 = matrix(r12, nrow = p1, ncol = p2)
  R13 = matrix(r13, nrow = p1, ncol = p3)
  R23 = matrix(r23, nrow = p2, ncol = p3)
  
  R = rbind( cbind(R1, R12, R13), 
             cbind(t(R12), R2, R23), 
             cbind(t(R13), t(R23), R3) )
  
  # Make sure its is positive definite
  eigs = eigen(R)
  if(any(eigs$values < 0)){
    cat("Warning: correlation matrix is not positive definite. \n")
  }
  eigs$values[eigs$values < 0] = 0
  R = eigs$vectors %*% diag(eigs$values) %*% t(eigs$vectors)
  
  return(R)
}

# -------------------------------------------------
# Compute the dissimalarity matrix C_XY and
# theoretical moments for two kernels
# -------------------------------------------------

compute_C = function(K_X, K_Y) {
  
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


# -------------------------------------------------
# Compute the 2 partial dissimalarity covariances 
# C_XY_Z and nC_XY_Z and residuals as well as 
# theoretical moments for two kernels vs a third
# -------------------------------------------------

compute_partial_C = function(K_X, K_Y, K_Z) {
  
  # Get n
  n = nrow(K_X)
  
  # Compute traces
  tr_K_X = sum(diag(K_X))
  tr_K_Y = sum(diag(K_Y))
  tr_K_Z = sum(diag(K_Z))
  
  # Compute delta2
  delta2_X = sum(diag(K_X%*%K_X)) 
  delta2_Y = sum(diag(K_Y%*%K_Y))
  delta2_Z = sum(diag(K_Z%*%K_Z))
  
  # Compute C_XY, C_XZ, C_YZ
  C_XY = sum(diag(K_X %*% K_Y))
  C_XZ = sum(diag(K_X %*% K_Z))
  C_YZ = sum(diag(K_Y %*% K_Z))
  
  # Compute the centered covariances 
  Cc_XY = C_XY - tr_K_X * tr_K_Y / (n-1)
  Cc_XZ = C_XZ - tr_K_X * tr_K_Z / (n-1)
  Cc_YZ = C_YZ - tr_K_Y * tr_K_Z / (n-1)
  
  # Cc_ZZ
  Cc_ZZ = delta2_Z - tr_K_Z^2 / (n-1)
  
  # effective dimensionality
  v_X = tr_K_X^2 / delta2_X 
  v_Y = tr_K_Y^2 / delta2_Y
  v_Z = tr_K_Z^2 / delta2_Z
  
  # mean effective dimensionality
  m_v_X = n - 1 - v_X
  m_v_Y = n - 1 - v_Y
  m_v_Z = n - 1 - v_Z
  
  # Kappa
  k = 1 / ((n-2)*(n-1)*(n+1))
  
  # Compute the partial covariance and theoretical moments
  init_c_xy_z = (Cc_XZ * Cc_YZ) / Cc_ZZ
  C_XY_Z = init_c_xy_z + tr_K_X * tr_K_Y / (n-1)
  E_C_XY = tr_K_X * tr_K_Y / (n-1)
  Var_C_XY = 4 * k^2 * delta2_X * delta2_Y * m_v_X * m_v_Y 
  
  # Compute the residual covariance and theoretical moments
  C_XYrZ = Cc_XY - init_c_xy_z
  E_C_XYrZ = 0
  Var_C_XYrZ = (n^2 - n - 4) * Var_C_XY / 2
  
  # Compute the naive partial covariance and theoretical moments
  nC_XY_Z = (C_XZ * C_YZ) / C_ZZ
  E_nC_XY_Z = tr_K_X * tr_K_Y * v_Z / (n-1)^2
  Var_nC_XY_Z = 2 * k * delta2_X * delta2_Y * m_v_Z / (n-1)^3 *
    (2*(n-1)*k*m_v_X*m_v_Y*m_v_Z + m_v_X*v_Y*v_Z + v_X*m_v_Y*v_Z)
  
  # Compute the naive residual covariance and theoretical moments
  nC_XYrZ = C_XY - nC_XY_Z
  E_nC_XYrZ = tr_K_X * tr_K_Y * m_v_Z / (n-1)^2
  Var_nC_XYrZ = 2 * k * delta2_X * delta2_Y / (n-1)^3 *
    (2*(n-1)*k*m_v_X*m_v_Y*m_v_Z*(m_v_Z - 2*(n-1)) + (n-1)^2*m_v_X*m_v_Y + 
       (v_X*m_v_Y + m_v_X*v_Y)*v_Z*m_v_Z)
  
  return(list(C_XY_Z = C_XY_Z,
              E_C_XY_Z = E_C_XY,
              Var_C_XY_Z = Var_C_XY,
              C_XYrZ = C_XYrZ,
              E_C_XYrZ = E_C_XYrZ,
              Var_C_XYrZ = Var_C_XYrZ,
              nC_XY_Z = nC_XY_Z,
              E_nC_XY_Z = E_nC_XY_Z,
              Var_nC_XY_Z = Var_nC_XY_Z,
              nC_XYrZ = nC_XYrZ,
              E_nC_XYrZ = E_nC_XYrZ,
              Var_nC_XYrZ = Var_nC_XYrZ))
}



