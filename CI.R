source("utility.R")


refinement2 <- function(detected_CP, A, B, rank, verbose = FALSE) {
  K <- length(detected_CP)
  if(K == 0) {return(detected_CP)}
  
  nu <- c(0, detected_CP, dim(A)[1])
  eta_bar <- nu
  for (k in 2:(K+1)) {
    sk <- floor((nu[k-1]+nu[k])/2)
    ek <- ceiling((nu[k+1]+nu[k])/2)
    # sk <- floor(9*nu[k-1]/10 + nu[k]/10)
    # ek <- ceiling(nu[k]/10 + 9*nu[k+1]/10)
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

construct_CI <- function(alpha, detected_CP, A, B, rank, verbose = FALSE) {
  # detected CP is after first refinement
  # second refinement done within this function 
  
  K <- length(detected_CP)
  T <- dim(A)[1]
  intervals <- matrix(nrow = K, ncol = 4) # estimate, refinement, lb, ub 

  b <- refinement2(detected_CP, A, B, rank, verbose)
  nu <- c(0, detected_CP, T) 
  
  for (k in 2:(K+1)) {
    if (verbose) {
      cat("k = ", k, ", nu[c(k-1, k, k+1)] = ", nu[c(k-1, k, k+1)], ".\n")
    }
    # Step 1: 
    B_nu_prev <- (1/(nu[k]-nu[k-1])) * as.tensor( apply(B[(nu[k-1]+1):nu[k], , , , drop = FALSE], c(2, 3, 4), sum) )
    B_nu_next <- (1/(nu[k+1]-nu[k])) * as.tensor( apply(B[(nu[k]+1):nu[k+1], , , , drop = FALSE], c(2, 3, 4), sum) )
    
    P_nu_prev <- estimate_thpca(B_nu_prev, rank, tmax = 20)
    P_nu_next <- estimate_thpca(B_nu_next, rank, tmax = 20)
    
    k_k <- sum((P_nu_next - P_nu_prev)^2)^0.5
    Psi_k <- 1/k_k * (P_nu_next - P_nu_prev)
    
    # Step 2: 
    var_k <- c(nu)
    # Question, should this be K+1 or K+2
    # I think K+2, since this corresponds to K+1 in the paper. nu[k_prime] would be T
    for (k_prime in 2:(K+2)) {
      # if (verbose) {
      #   cat("\tk' = ", k_prime, ", nu[c(k'-1, k')] = ", nu[c(k_prime-1, k_prime)], ".\n")
      # }
      B_prime <- (1/(nu[k_prime]-nu[k_prime-1])) * as.tensor( apply(B[(nu[k_prime-1]+1):nu[k_prime], , , , drop = FALSE], c(2, 3, 4), sum) )
      P_prime <- estimate_thpca(B_prime, rank, tmax = 20)
      
      var_k[k_prime] <- (1 / (nu[k_prime] - nu[k_prime - 1] - 1)) * sum(
        apply(A[(nu[k_prime - 1] + 1):nu[k_prime], , , , drop = FALSE], 1, function(x) {
          sum(Psi_k * (x - P_prime))^2
        })
      )
    }
    
    if (verbose) {
      cat("k = ", k, ", sigma^2_k[c(k, k+1)] = ", var_k[c(k, k+1)], ".\n")
    }
    
    # Step 3: 
    # Question, what is B, M and N in practice? 
    big_B <- 500 # Number of bootstrap replicates
    M <- T       # Support for r
    N <- T       # Length of interval used for simulation (can customize)
    grid_r <- seq(-M, M, length.out = 201)
    u_hat_k <- numeric(big_B)
    
    for (index in 1:big_B) {
      z_b <- rnorm(length(-floor(N * M):floor(N * M)))
      P_vals <- sapply(grid_r, function(r) {
        if (r < 0) {
          return(-r + (2 * sqrt(var_k[k]) / sqrt(N)) * sum(z_b[ceiling(N * r):(-1) + (floor(N * M) + 1)]))
        } else if (r > 0) {
          return(r + (2 * sqrt(var_k[k+1]) / sqrt(N)) * sum(z_b[1:floor(N * r) + (floor(N * M) + 1)]))
        } else {
          return(0)
        }
      })
      u_hat_k[index] <- grid_r[which.min(P_vals)]
    }
    
    # hist(u_hat_k, breaks = 50, main = "Distribution of u_hat_k", xlab = "u_hat_k")
    # Step 4: 
    # Question: why is this b_k in the paper? 
    ci_upper <- b[k-1] - ifelse(k_k^2 == 0, 0, quantile(u_hat_k, probs = alpha / 2) / (k_k^2))
    ci_lower <- b[k-1] - ifelse(k_k^2 == 0, 0, quantile(u_hat_k, probs = 1 - alpha / 2) / (k_k^2))
    
    if (verbose) {
      # cat("Range of u_k = ", range(u_hat_k), ". Quantiles : ", c(quantile(u_hat_k, probs = alpha / 2), quantile(u_hat_k, probs = 1 - alpha / 2)), ".\n")
      # cat("Range of P_vals = ", range(P_vals), ".\n")
      cat("b[k] = ", b[k-1], ", Interval = ", c(ci_lower, ci_upper), ".\n\n")
    }
    
    intervals[k-1, ] <- c(nu[k], b[k-1], ci_lower, ci_upper)
  }  
  intervals
}

# coverage <- function(true_CP, locs, starts, ends) {
#   if (length(true_CP) != length(starts)) {
#     warning("Incorrect number of intervals supplied.")
#     # return(NA)
#   }
#   
#   covered <- logical(length(true_CP))
#   
#   for (i in seq_along(true_CP)) {
#     cp <- true_CP[i]
#     covered[i] <- any(starts <= cp & cp <= ends)
#   }
#   
#   match to middle
#   covered \in {-1, 0, 1}
#   lengths (again negative )
#   
#   return(covered, lengths)
# }


coverage <- function(true_CP, locs, starts, ends) {
  covered <- rep(-1L, length(true_CP))  # default to -1 (no match)
  lengths <- rep(-1L, length(true_CP))  # default to -1 (no match)
  
  # Case 1: fewer locs than true CPs => only evaluate closest true CPs
  if (length(locs) < length(true_CP)) {
    available <- rep(TRUE, n_true)
    
    for (i in 1:length(locs)) {
      # Match each loc to closest available true_CP
      dists <- abs(true_CP - locs[i])
      dists[!available] <- Inf  # already matched
      min_index <- which.min(dists)
      
      if (is.finite(dists[min_index])) {
        available[min_index] <- FALSE
        if (starts[i] <= true_CP[min_index] && true_CP[min_index] <= ends[i]) {
          covered[min_index] <- 1L
          lengths[min_index] <- ends[i] - starts[i]
        } else {
          covered[min_index] <- 0L
          lengths[min_index] <- ends[i] - starts[i]
        }
      }
    }
  } # Case 2: enough locs to evaluate every true CP
  else {
    for (i in 1:length(true_CP)) {
      # Find closest loc
      dists <- abs(locs - true_CP[i])
      j <- which.min(dists)
      
      if (starts[j] <= true_CP[i] && true_CP[i] <= ends[j]) {
        covered[i] <- 1L
        lengths[i] <- ends[j] - starts[j]
      } else {
        covered[i] <- 0L
        lengths[i] <- ends[j] - starts[j]
      }
    }
  }
  
  return(list(covered = covered, lengths = lengths))
}
