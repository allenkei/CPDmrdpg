#library(devtools) # install.packages("devtools")
#install_github("etam4260/kneedle") # install the package "kneedle" via "devtools"
library(kneedle)

cal_log_lik <- function(Y, detected_CP, hat.rank){
  # Y: data with size (T, n, n, l)
  # detected_CP: vector of change points (cp1, cp2, ..., cpk) from a particular threshold # i.e. c(50,100)
  # hat.rank: c(15,15,15) for the function 'estimate_thpca()'
  
  num_T <- dim(Y)[1]     # T
  n_node <- dim(Y)[2]    # n
  num_layer <- dim(Y)[4] # l
  log_lik_full <- 0 # BIC_first_term = -2 * log_lik_full

  CP <- c(detected_CP, num_T) # add the last time point # i.e. c(50, 100, 150)
  
  start_t <- 0
  for (i in 1:length(CP)) {
    end_t <- CP[i] 
    
    A_within <- A.tensor[(start_t+1):end_t, , , ] # (end_t-start_t, n, n, l)
    dim(A_within)
    
    tensor_sum <- (1/(end_t-start_t)) * as.tensor( apply(A_within, c(2, 3, 4), sum) ) # (n, n, l)
    P_matrix  <- estimate_thpca(tensor_sum, hat.rank, tmax = 20) # (n, n, l) # fixed within interval
    # Probably getting bad estimates for P when not enough data 
    
    eps <- 1e-6
    
    for(t_within in 1:dim(A_within)[1]){
      # log_lik_t <- sum(A_within[t_within,,,]*log(P_matrix) + (1-A_within[t_within,,,]) * log(1-P_matrix)) # here has NaN
      log_lik_t <- sum(A_within[t_within,,,] * log(pmax(P_matrix, eps)) +
                         (1 - A_within[t_within,,,]) * log(pmax(1 - P_matrix, eps)))
      log_lik_full <- log_lik_full + log_lik_t
    }
    
    start_t <- end_t # make sure to update start_t
  }
  
  # cat("Candidates: ", paste(detected_CP, collapse = ", "), ". log-Likelihood = ", log_lik_full, "\n", sep = "")
  
  return(log_lik_full)
}

BIC <- function(Y, detected_CP, verbose, hat.rank){
  
  # BIC = −2log_lik + log(T*N_net) × num_par_within × (K+1)
  # Y: data with size (T, n, n, l)
  # detected_CP: vector of change points (cp1, cp2, ..., cpk) from a particular threshold # i.e. c(50,100)
  # hat.rank: c(15,15,15) for the function 'estimate_thpca()'
  
  num_T <- dim(Y)[1]     # T
  n_node <- dim(Y)[2]    # n
  num_layer <- dim(Y)[4] # l
  log_lik_full <- cal_log_lik(Y, detected_CP, hat.rank)
  
  # Calculate the second term in BIC
  num_par_within <- n_node*hat.rank[1]*hat.rank[2] + n_node*hat.rank[2]*hat.rank[3] + n_node*hat.rank[1]*hat.rank[3]
  num_data <- num_T * n_node * n_node * num_layer # Why not just num_T? 
  
  BIC <- (-2) * log_lik_full + log(num_data) * num_par_within * (length(detected_CP)+1)
  if(verbose) {
    cat("Candidates: ", paste(detected_CP, collapse = ", "), ". BIC = ", BIC, ". log-Likelihood = ", log_lik_full, "\n", sep = "")
  }
  
  return(BIC) # choose the threshold (and corresponding results) with lowest BIC
}

MDL <- function(Y, detected_CP, verbose, hat.rank){
  
  # 1. The penalty of a real-valued parameter estimated by n data points is log2 n;
  # 2. The penalty of an unbounded integer parameter K is 2 log2 K;
  # 3. The penalty of an integer bounded by a known integer N is 2 log2 N.
  # https://www.mdpi.com/1099-4300/26/1/50
  # Y: data with size (T, n, n, l)
  # detected_CP: vector of change points (cp1, cp2, ..., cpk) from a particular threshold # i.e. c(50,100)
  # hat.rank: c(15,15,15) for the function 'estimate_thpca()'
  
  num_T <- dim(Y)[1]     # T
  n_node <- dim(Y)[2]    # n
  num_layer <- dim(Y)[4] # l
  log_lik_full <- cal_log_lik(Y, detected_CP, hat.rank)
  
  K <- length(detected_CP)
  p <- n_node*hat.rank[1]*hat.rank[2] + n_node*hat.rank[2]*hat.rank[3] + n_node*hat.rank[1]*hat.rank[3]
  
  penalty <- 2 * log(max(K, 1)) + 2*K*log(num_T * n_node * n_node * num_layer) 
  CP <- c(detected_CP, num_T) # add the last time point # i.e. c(50, 100, 150)
  
  start_t <- 0
  for (i in 1:length(CP)) {
    
    end_t <- CP[i] 
    penalty <- penalty + p*log(n_node * n_node * num_layer * (end_t - start_t))
    start_t <- end_t # make sure to update start_t
  }
  
  MDL <- (-2) * log_lik_full + penalty
  if(verbose) {
    cat("Candidates: ", paste(detected_CP, collapse = ", "), ". MDL = ", MDL, ". log-Likelihood = ", log_lik_full, "\n", sep = "")
  }

  return(MDL) # choose the threshold (and corresponding results) with lowest MDL
}

AIC <- function(Y, detected_CP, verbose, hat.rank){
  
  # AIC = −2log_lik + 2 x (T*N_net) × num_par_within × (K+1)
  # Y: data with size (T, n, n, l)
  # detected_CP: vector of change points (cp1, cp2, ..., cpk) from a particular threshold # i.e. c(50,100)
  # hat.rank: c(15,15,15) for the function 'estimate_thpca()'
  
  num_T <- dim(Y)[1]     # T
  n_node <- dim(Y)[2]    # n
  num_layer <- dim(Y)[4] # l
  log_lik_full <- cal_log_lik(Y, detected_CP, hat.rank)
  
  num_par_within <- n_node*hat.rank[1]*hat.rank[2] + n_node*hat.rank[2]*hat.rank[3] + n_node*hat.rank[1]*hat.rank[3]
  
  AIC <- (-2) * log_lik_full + 2 * num_par_within * (length(detected_CP)+1)
  if(verbose) {
    cat("Candidates: ", paste(detected_CP, collapse = ", "), ". AIC = ", AIC, ". log-Likelihood = ", log_lik_full, "\n", sep = "")
  }
  
  return(AIC) # choose the threshold (and corresponding results) with lowest AIC
}

l1 <- function(Y, detected_CP, verbose, hat.rank, beta = 2){
  
  # BIC = −2log_lik + 2 sum || hat{P}^{s,t}, hat{P}^{s+1,t} ||_1
  # Y: data with size (T, n, n, l)
  # detected_CP: vector of change points (cp1, cp2, ..., cpk) from a particular threshold # i.e. c(50,100)
  # hat.rank: c(15,15,15) for the function 'estimate_thpca()'
  
  num_T <- dim(Y)[1]     # T
  n_node <- dim(Y)[2]    # n
  num_layer <- dim(Y)[4] # l
  log_lik_full <- cal_log_lik(Y, detected_CP, hat.rank)

  
  penalty <- 0 
  if (length(detected_CP) > 0) {
    CP <- c(detected_CP, num_T) # add the last time point # i.e. c(50, 100, 150)
    start_t <- 0
    cp_t <- CP[1]
    
    A_bottom <- A.tensor[(start_t+1):cp_t, , , ] # (end_t-start_t, n, n, l)
    sum_bottom <- (1/(cp_t-start_t)) * as.tensor( apply(A_bottom, c(2, 3, 4), sum) )
    P_bottom  <- estimate_thpca(sum_bottom, hat.rank, tmax = 20)
    
    for (i in 2:length(CP)) {
      end_t <- CP[i] 
      
      A_top <- A.tensor[(cp_t+1):end_t, , , ] # (end_t-start_t, n, n, l)
      sum_top <- (1/(end_t-cp_t)) * as.tensor( apply(A_top, c(2, 3, 4), sum) )
      P_top  <- estimate_thpca(sum_top, hat.rank, tmax = 20)
      # Probably getting bad estimates for P when not enough data 
      penalty <- penalty + sum(abs(P_top - P_bottom))
      
      # Updates
      start_t <- cp_t
      cp_t <- end_t
      A_bottom <- A_top
      sum_bottom <- sum_top
      P_bottom <- P_top
    }
  }
  if(verbose) {  
    cat("Candidates: ", paste(detected_CP, collapse = ", "), 
        ". Penalty = ", penalty, ". log-Likelihood = ", log_lik_full, "\n", sep = "")
  }

  return(-log_lik_full + beta * penalty) # choose the threshold (and corresponding results) with lowest BIC
}

elbow <- function(obj, cps, verbose, hat.rank) {
  # https://raghavan.usc.edu/papers/kneedle-simplex11.pdf
  log_lik <- cal_log_lik(obj, cps, hat.rank)
  if (verbose) {cat("Candidates: ", paste(cps, collapse = ", "), ". log-Likelihood = ", log_lik, "\n", sep = "")}
  
  return(log_lik)
}

model_selection <- function(results, obj, method = "BIC", verbose = "TRUE", ...) {
  # Using results constructed from Seeded Binary Segmentation, 
  # calculates model selection statistic, and selects best model.
  # The best model is the highest threshold minimizer of the statistic.
  # @param results Output of seeded binary segmentation
  # @param obj Data or model input needed for the method
  # @param method Function to calculate model selection statistic (default: cal_BIC)
  # @param ... Additional arguments for the selection method
  # FIXED TO RUN WITH STEP 1: MULTIPLYING CANDIDATES BY 2
  # TODO: add elbow as an option
  fun <- match.fun(method)

  best_cps <- NULL
  best_stat <- Inf  # Set initial stat to a high value
  best_index = -1
  stats <- numeric(length(results) - 1)
  num_cps <- numeric(length(results) - 1)
  
  for (i in 2:length(results)) {
    # Print Statement Handled Inside fun()
    cps <- 2 * sort(results[[i]]$results[, 1])
    stats[i-1] <- fun(obj, cps, verbose, ...)[1] # This allows other methods to return more arguments 
    num_cps[i-1] <- length(cps)
    
    # Store best model based on minimum Stat
    if (stats[i-1] < best_stat) {
      best_stat <- stats[i-1]
      best_cps <- cps
      best_index <- i 
    }
  }
  
  if (method == "elbow") {
    # Fancy calculation of "best"
    # https://raghavan.usc.edu/papers/kneedle-simplex11.pdf
    
    # unique(stats[num_cps < 15]), unique(num_cps[num_cps < 15])
    knee <- kneedle(1:length(unique(stats)), unique(stats), concave = TRUE, decreasing = FALSE)
    best_stat <- knee[2]
    best_index <- which(stats == best_stat)[1] + 1
    best_cps <- 2 * sort(results[[best_index]]$results[, 1])
  } 
  
  return(list("candidates" = as.vector(best_cps), "stat" = best_stat, "threshold" = results[[best_index]]$threshold, 
              "range" = stats, "ncps" = num_cps))
}
