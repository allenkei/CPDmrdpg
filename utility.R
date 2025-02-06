cal_BIC <- function(Y, detected_CP, hat.rank){
  
  # BIC = −2log_lik + log(T*N_net) × num_par_within × (K+1)
  # Y: data with size (T, n, n, l)
  # detected_CP: vector of change points (cp1, cp2, ..., cpk) from a particular threshold # i.e. c(50,100)
  # hat.rank: c(15,15,15) for the function 'estimate_thpca()'
  # no need to enter this function if detected_CP is empty
  
  num_T <- dim(Y)[1]     # T
  n_node <- dim(Y)[2]    # n
  num_layer <- dim(Y)[4] # l
  log_lik_full <- 0 # BIC_first_term = -2 * log_lik_full
  
  # Calculate the second term in BIC
  num_par_within <- n_node*hat.rank[1]*hat.rank[2] + n_node*hat.rank[2]*hat.rank[3] + n_node*hat.rank[1]*hat.rank[3]
  num_data <- log(num_T * n_node * n_node * num_layer)
  BIC_second_term <- num_data * num_par_within * (length(detected_CP)+1)
  
  # Calculate the first term in BIC
  CP <- c(detected_CP, num_T) # add the last time point # i.e. c(50, 100, 150)

  start_t <- 1
  for (i in 1:length(CP)) {
    end_t <- CP[i] 
    
    A_within <- A.tensor[start_t:end_t, , , ] # (end_t-start_t, n, n, l)
    dim(A_within)
    
    tensor_sum <- (1/(end_t-start_t)) * as.tensor( apply(A_within, c(2, 3, 4), sum) ) # (n, n, l)
    P_matrix  <- estimate_thpca(tensor_sum, hat.rank, tmax = 20) # (n, n, l) # fixed within interval
    
    for(t_within in 1:dim(A_within)[1]){
      log_lik_t <- sum(A_within[t_within,,,]*log(P_matrix) + (1-A_within[t_within,,,]) * log(1-P_matrix))
      log_lik_full <- log_lik_full + log_lik_t
    }
    
    start_t <- CP[i] + 1 # make sure to update start_t
  }
  
  BIC <- (-2) * log_lik_full + BIC_second_term
  cat(log_lik_full, "\t", BIC_second_term, "\t")
  
  return(BIC) # choose the threshold (and corresponding results) with lowest BIC
}


model_selection <- function(results, obj, method = cal_BIC, ...) {
  # Using results constructed from Seeded Binary Segmentation, 
  # calculates model selection statistic, and selects best model.
    # The best model is the highest threshold minimizer of the statistic.
  # @param results Output of seeded binary segmentation
  # @param obj Data or model input needed for the method
  # @param method Function to calculate model selection statistic (default: cal_BIC)
  # @param ... Additional arguments for the selection method
  # FIXED TO RUN WITH STEP 1: MULTIPLYING CANDIDATES BY 2
  method <- match.fun(method)
  
  best_cps <- NULL
  best_BIC <- Inf  # Set initial BIC to a high value
  best_index = -1
  
  for (i in 2:length(results)) {
    cps <- 2 * sort(results[[i]]$results[, 1])
    BIC <- method(obj, cps, ...)
    
    cat("Candidates: ", paste(cps, collapse = ", "), ". BIC = ", BIC, "\n", sep = "")
    
    # Store best model based on minimum BIC
    if (!is.na(BIC) && BIC < best_BIC) {
      best_BIC <- BIC
      best_cps <- cps
      best_index <- i 
    }
  }
  
  return(list(Candidates = as.vector(best_cps), BIC = best_BIC, threshold = results[[best_index]]$threshold))
}
