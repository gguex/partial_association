# ----------------------------------------------------
#
# Functions for the projects
#
# ----------------------------------------------------

library(MASS)

# Function to generate 2 datasets 

generate_cov2dataset = function(n1, n2, r1, r2, r12, 
                                s1_min=0, s1_max=1, s2_min=0, s2_max=1) {
  
  # Create correlation matrix
  R1 = matrix(r1, nrow = n1, ncol = n1)
  diag(R1) = 1
  R2 = matrix(r2, nrow = n2, ncol = n2)
  diag(R2) = 1
  R12 = matrix(r12, nrow = n1, ncol = n2)
  
  R = rbind( cbind(R1, R12), cbind(t(R12), R2) )
  
  # Generate variances
  s1 = runif(n1, s1_min, s1_max)
  s2 = runif(n2, s2_min, s2_max)
  
  Cov = R * outer(c(s1, s2), c(s1, s2))
  
  return(Cov)
  
}