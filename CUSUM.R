##########################
# Change Point Detection #
##########################
source("utility.R")

# Contains functions required to compute CUSUM Statistics
# And current styles of CUSUM Statistics considered

diff_frobenius <- function(A, B) sum((A - B)^2)^0.5

Hetero_PCA_test <- function(Y, r, tmax = 20, vartol = 1e-6){
  
  N_t = Y
  r = min(c(r, dim(N_t)))
  diag(N_t) = 0
  U_t = matrix(NA, nrow = dim(Y)[1], r)
  t = 1
  approx = -1
  
  while(t <= tmax){ # Stop criterion: convergence or maximum number of iteration reached
    temp = svd(N_t)
    U_t = temp$u[,1:r]
    V_t = temp$v[,1:r]
    if (r > 1){
      tilde_N_test_t = U_t %*% diag(temp$d[1:r]) %*% t(V_t)
    }
    else{
      tilde_N_test_t = temp$d[1] * U_t %*% t(V_t)
    }
    
    N_test_new = N_t
    diag(N_test_new) = diag(tilde_N_test_t)
    N_t = N_test_new
    svector = diag(tilde_N_test_t)
    if (abs(sum(svector^2) - approx) > vartol){
      t = t+1
      approx = sum(svector^2)
    }
    else {
      break
    }
  }
  return(U_t)
}

Tensor_Hetero_PCA_test <- function(Y, r, tmax = 20){
  
  p = dim(Y)
  d = length(p)
  U_0 = list()
  for (i in 1:d){
    MY = k_unfold(Y, i)@data
    MY_Hetero_PCA = Hetero_PCA_test(MY %*% t(MY), r[i], tmax)
    U_0 = c(U_0, list(MY_Hetero_PCA))
  }
  return(U_0)
}

estimate_thpca <- function(Y.tensor, hat.rank, ub = 1, lb = 0, tmax = 20){
  U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank, tmax)
  P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
  P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
  P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
  Y.hat = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  P_hat  = Y.hat@data
  
  P_hat[P_hat > ub]  = ub - 1e-10
  P_hat[P_hat < lb]  = lb + 1e-10
  
  return(P_hat)
}

####################
# CUSUM Statistics #
####################

CUSUM_frobenius <- function(obj, s, e, t, rank, verbose = TRUE) {
  # Frobenius norm hat{P}^{s,t}, hat{P}^{s+1,t}
  # hat{P}^{a,b} = TH-PCA((b-a)^{-1}  \sum_{u =a}^b A(u) ,(d,d,m))
  if (verbose) {print(paste0("s = ", s, ", e = ", e, ", t = ", t, "."))}
  
  sum_s_t  <- (1/(t-s)) * as.tensor( apply(obj[(s+1):t, , , , drop = FALSE], c(2, 3, 4), sum) )
  sum_t_e <- (1/(e-t)) * as.tensor( apply(obj[(t+1):e, , , , drop = FALSE], c(2, 3, 4), sum) )
  
  P_s_t  <- estimate_thpca(sum_s_t, rank, tmax = 20)
  P_t_e <- estimate_thpca(sum_t_e, rank, tmax = 20)
  
  return(diff_frobenius(P_s_t, P_t_e))
}

CUSUM_layer <- function(obj, s, e, t, rank, verbose = TRUE) {
  # Layer-wise Frobenius norm 
  # max_l hat{P}^{s,t}_(l), hat{P}^{s+1,t}_(l)
  # hat{P}^{a,b} = TH-PCA((b-a)^{-1}  \sum_{u =a}^b A(u) ,(d,d,m))
  if (verbose) {print(paste0("s = ", s, ", e = ", e, ", t = ", t, "."))}
  
  sum_s_t  <- (1/(t-s)) * as.tensor( apply(obj[(s+1):t, , , , drop = FALSE], c(2, 3, 4), sum) )
  sum_t_e <- (1/(e-t)) * as.tensor( apply(obj[(t+1):e, , , , drop = FALSE], c(2, 3, 4), sum) )
  
  P_s_t  <- estimate_thpca(sum_s_t, rank, tmax = 20)
  P_t_e  <- estimate_thpca(sum_t_e, rank, tmax = 20)
  
  frobenius_diffs <- sapply(1:dim(obj)[4] , function(l) {
    diff_frobenius(P_s_t[, , l], P_t_e[, , l])
  })
  
  return(max(frobenius_diffs))
}

CUSUM_frob_SBS <- function(obj, s, e, t, rank, verbose = TRUE) {
  # Frobenius norm Using Weighting suggested in SBS
  # || sqrt((e-t)/((e-s)(t-s))) hat{P}^{s,t}, sqrt((t-s)/((e-s)(e-t))) hat{P}^{s+1,t}||
  # Using Weighting suggested in SBS
  # hat{P}^{a,b} = TH-PCA(\sum_{u =a}^b A(u) ,(d,d,m))
  if (verbose) {print(paste0("s = ", s, ", e = ", e, ", t = ", t, "."))}
  
  sum_s_t  <- as.tensor( apply(obj[(s+1):t, , , , drop = FALSE], c(2, 3, 4), sum) )
  sum_t_e <- as.tensor( apply(obj[(t+1):e, , , , drop = FALSE], c(2, 3, 4), sum) )
  
  P_s_t  <- sqrt((e-t)/(e-s)/(t-s)) * estimate_thpca(sum_s_t, rank, tmax = 20)
  P_t_e <- sqrt((t-s)/(e-s)/(e-t)) * estimate_thpca(sum_t_e, rank, tmax = 20)
  
  return(diff_frobenius(P_s_t, P_t_e))
}

CUSUM_step1 <- function(obj, s, e, t, obj.B, verbose = TRUE) {
  # Tensor Weighted Inner Product
  # Two identical copies A and B (e.g. odd and even indices)
  # tilde{A}^{s,e}(t) = sum_{u=s+1}^{e} omega_{s,e}^t(u) A(u), same for B
  # omega_{s,e}^t(u) = ifelse(u in (s+1):t, sqrt((e-t)/((e-s)*(t-s))), -sqrt((t-s)/((e-s)*(e-t)))
  # hat{D}^{s,e}(t) = \sum_{i,j,l=1}^{n,n,L} tilde{A}^{s,e}(t)_{i,j,l} tilde{B}^{s,e}(t)_{i,j,l} 
  if (verbose) {print(paste0("s = ", s, ", e = ", e, ", t = ", t, "."))}

  if (!all(dim(obj) == dim(obj.B))) {
    stop("Error: obj (A) and obj.B (B) must have the same dimensions.")
  }
  
  # Compute weights
  omega <- numeric(e - s)
  for (u in (s+1):t) {
    omega[u - s] <- sqrt((e - t) / ((e - s) * (t - s)))
  }
  for (u in (t+1):e) {
    omega[u - s] <- -sqrt((t - s) / ((e - s) * (e - t)))
  }

  weighted_A <- array(0, dim = dim(obj)[2:4])
  weighted_B <- array(0, dim = dim(obj)[2:4])
  
  for (u in (s+1):e) {
    weighted_A <- weighted_A + omega[u - s] * obj[u, , , ]
    weighted_B <- weighted_B + omega[u - s] * obj.B[u, , , ]
  }
  
  # Compute the CUSUM statistic as an inner product sum
  D_hat <- abs(sum(weighted_A * weighted_B))
  return(D_hat)
}

####################
# Local Refinement #
####################


CUSUM_refinement <- function(obj, s, e, t, obj.B, nu, rank) {
  # CUSUM local refinement (Def)
  # Two identical copies A and B (e.g. odd and even indices)
  # tilde{A}^{s,e}(t) = sum_{u=s+1}^{e} omega_{s,e}^t(u) A(u), same for B
  # omega_{s,e}^t(u) = ifelse(u in (s+1):t, sqrt((e-t)/((e-s)*(t-s))), -sqrt((t-s)/((e-s)*(e-t)))
  # hat{P}^{s,e}(nu) = estimate_thpca(tilde{A}^{s,e}(nu), rank, ub, lb)
  # hat{D}^{s,e}(t) = ||hat{P}^{s,e}(nu)||^{-1} \sum_{i,j,l=1}^{n,n,L} tilde{A}^{s,e}(t)_{i,j,l} hat{P}^{s,e}(nu)_{i,j,l} 
  
  if (!all(dim(obj) == dim(obj.B))) {
    stop("Error: obj (A) and obj.B (B) must have the same dimensions.")
  }
  
  # tilde A
  omega_A <- numeric(e - s)
  for (u in (s+1):t) {
    omega_A[u - s] <- sqrt((e - t) / ((e - s) * (t - s)))
  }
  for (u in (t+1):e) {
    omega_A[u - s] <- -sqrt((t - s) / ((e - s) * (e - t)))
  }
  
  omega_B <- numeric(e - s)
  for (u in (s+1):nu) {
    omega_B[u - s] <- sqrt((e - nu) / ((e - s) * (nu - s)))
  }
  for (u in (nu+1):e) {
    omega_B[u - s] <- -sqrt((nu - s) / ((e - s) * (e - nu)))
  }
  
  weighted_A <- array(0, dim = dim(obj)[2:4])
  weighted_B <- array(0, dim = dim(obj.B)[2:4])
  for (u in (s+1):e) {
    weighted_A <- weighted_A + omega_A[u - s] * obj[u, , , ]
    weighted_B <- weighted_B + omega_B[u - s] * obj[u, , , ]
  }
  
  P_hat <- estimate_thpca(as.tensor(weighted_B), rank, 2*sqrt((e-nu)*(nu-s)/(e-s)), -2*sqrt((e-nu)*(nu-s)/(e-s)))
  D_hat <- 1/(sum((P_hat)^2)^0.5)*abs(sum(P_hat * weighted_A))
  return(D_hat)
}

refinement1 <- function(detected_CP, A, B, rank, verbose = FALSE) {
  K <- length(detected_CP)
  if(K == 0) {return(detected_CP)}
  
  nu <- c(0, detected_CP, dim(A)[1])
  b <- nu
  
  for (k in 2:(K+1)) {
    sk <- floor((nu[k-1]+nu[k])/2)
    ek <- ceiling((nu[k+1]+nu[k])/2)
    if (verbose) {
      cat("k = ", k, ", nu[c(k-1, k, k+1)] = ", nu[c(k-1, k, k+1)], ", sk = ", sk, ", ek = ", ek, ".\n")
    }
    
    max_gain <- 0
    
    if ((ek-sk) > 1) {
      for (t in (sk+1):(ek-1)) {
        gain <- CUSUM_refinement(A, sk, ek, t, B, nu[k], rank)
        if (verbose) {
          cat("\tsk = ", sk, ", ek = ", ek, ", t = ", t, ", gain = ", gain, ".\n")
        }
        if (gain > max_gain) {
          max_gain <- gain
          b[k] <- t
        }
      }
    } 
  }
  
  return(b[2:(K+1)])
}

refinement2 <- function(detected_CP, A, B, rank, verbose = FALSE) {
  K <- length(detected_CP)
  if(K == 0) {return(detected_CP)}
  
  nu <- c(0, detected_CP, dim(A)[1])
  eta_bar <- nu
  for (k in 2:(K+1)) {
    sk <- floor(9*nu[k-1]/10 + nu[k]/10)
    ek <- ceiling(nu[k]/10 + 9*nu[k+1]/10)
    if (verbose) {
      cat("k = ", k, ", nu[c(k-1, k, k+1)] = ", nu[c(k-1, k, k+1)], ", sk = ", sk, ", ek = ", ek, ".\n")
    }
    
    min_Qk <- Inf
    
    if ((ek-sk) > 1) {
      for (eta in (sk+1):(ek-1)) {
        B_nu_prev <- (1/(nu[k]-nu[k-1])) * as.tensor( apply(B[(nu[k-1]+1):nu[k], , , , drop = FALSE], c(2, 3, 4), sum) )
        B_nu_next <- (1/(nu[k+1]-nu[k])) * as.tensor( apply(B[(nu[k]+1):nu[k+1], , , , drop = FALSE], c(2, 3, 4), sum) )
        
        P_nu_prev <- estimate_thpca(B_nu_prev, rank, tmax = 20)
        P_nu_next <- estimate_thpca(B_nu_next, rank, tmax = 20)
        
        Qk <- sum(c(vapply( (sk+1):eta, function(t) (diff_frobenius(A[t, , , ], P_nu_prev))^2, numeric(1) ), 
                    vapply( (eta+1):ek, function(t) (diff_frobenius(A[t, , , ], P_nu_next))^2, numeric(1) )))
        
        if (verbose) {
          cat("\tsk = ", sk, ", ek = ", ek, ", eta = ", eta, ", Qk = ", Qk, ".\n")
        }
        if (Qk < min_Qk) {
          min_Qk <- Qk
          eta_bar[k] <- eta
        }
      }
    } 
  }
  return(eta_bar[2:(K+1)])
  
}
    
