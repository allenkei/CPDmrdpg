
library(rTensor)




get_blockwise_const_mat <- function(n, n_c, p_1, p_2){
  # n: num node, n_c: num block, p1: between prob, p2: inter prob
  # output n by n matrix with same block size
  P = matrix(p_1, n, n)
  size_c = floor(n / n_c)
  
  for (k in 1:n_c){
    if (k < n_c){
      P[(1 + size_c*(k-1)):(size_c * k), (1 + size_c*(k-1)):(size_c * k)] = p_2  
    } else {
      P[(1 + size_c*(n_c-1)):n, (1 + size_c*(n_c-1)):n] = p_2  
    }
  }
  
  return(P)
}




get_blockwise_var_size_mat <- function(n, block_sizes, p_1, p_2){
  # output n by n matrix with block size varied
  # i.e. n = 10, block_sizes = c(3,7) where 3 + 7 = 10
  P = matrix(p_1, n, n)
  start_idx = 1
  
  for (k in 1:length(block_sizes)){
    size_k = block_sizes[k]
    
    end_idx = start_idx + size_k - 1
    if (end_idx > n) end_idx = n # ensure not exceed the matrix size
    
    P[start_idx:end_idx, start_idx:end_idx] = p_2  
    start_idx = end_idx + 1
  }
  
  return(P)
}





get_sbm_params <- function(n, L, n_c=c(4, 4), flip_layer=TRUE){
  
  # output probability matrix before and after, each with (n,n L)
  # block size fixed
  
  probability_1 = array(NA, c(n, n, L))
  probability_2 = array(NA, c(n, n, L))
  
  prob = seq(0,1,  1/(4*L))
  for (layer in 1: L)
  {
    p_1 = runif(1, prob[2*L+layer], prob[2*L+layer+1])
    p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
    
    P = get_blockwise_const_mat(n, n_c[1], p_1, p_2)
    probability_1[, , layer] = P
    
    P = get_blockwise_const_mat(n, n_c[2], p_1, p_2)
    
    if (flip_layer){
      probability_2[, , L - layer + 1] = P
    } else {
      probability_2[, , layer] = P
    }
    
  }
  
  return(list(probability_1, probability_2))
}





get_sbm_var_size_params <- function(n, L, block_sizes1, block_sizes2){
  
  # output probability matrix before and after, each with (n,n L)
  # block size varied
  
  probability_1 = array(NA, c(n, n, L))
  probability_2 = array(NA, c(n, n, L))
  
  prob = seq(0,1,  1/(4*L))
  for (layer in 1: L)
  {
    p_1 = runif(1, prob[2*L+layer], prob[2*L+layer+1])
    p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
    
    P = get_blockwise_var_size_mat(n, block_sizes1, p_1, p_2)
    probability_1[, , layer] = P
    
    P = get_blockwise_var_size_mat(n, block_sizes2, p_1, p_2)
    probability_2[, , layer] = P
  }
  
  return(list(probability_1, probability_2))
}






get_sbm_VS_FL_params <- function(n, L, block_sizes1, block_sizes2, n_c=c(4, 4)){
  # VS: varied block sizes
  # FL: first layer differs only # block_sizes1 and block_sizes2 apply to first layer
  
  # output probability matrix before and after, each with (n,n L)
  
  probability_1 = array(NA, c(n, n, L))
  probability_2 = array(NA, c(n, n, L))
  
  prob = seq(0,1,  1/(4*L))
  
  # iterate over layers
  for (layer in 1: L) 
  {
    
    if(layer == 1){ # first layer
      
      p_1 = runif(1, prob[2*L+layer], prob[2*L+layer+1])
      p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
      
      probability_1[, , layer] = get_blockwise_var_size_mat(n, block_sizes1, p_1, p_2)
      probability_2[, , layer] = get_blockwise_var_size_mat(n, block_sizes2, p_1, p_2)
      
    } else { # rest of the layers 
      
      p_1 = runif(1, prob[2*L+layer], prob[2*L+layer+1])
      p_2 = runif(1, prob[3*L+layer], prob[3*L+layer+1])
      
      probability_1[, , layer] = get_blockwise_const_mat(n, n_c[1], p_1, p_2)
      probability_2[, , layer] = get_blockwise_const_mat(n, n_c[2], p_1, p_2)
      
    }
    
  }
  
  return(list(probability_1, probability_2))
}

get_sbm_params_spa_inc <- function(n, L, n_c=c(4, 4, 4), epsilon = 0.05){
  # output probability matrix before and after, each with (n,n L)
  # block size fixed
  
  probability_1 = array(NA, c(n, n, L))
  probability_2 = array(NA, c(n, n, L))
  probability_3 = array(NA, c(n, n, L))
  
  # Define increasing probability ranges for each window
  prob_ranges <- list(
    c(0.21, 0.25),   # sparse
    c(0.21 + epsilon, 0.25 + epsilon),    # moderate
    c(0.21 + 2*epsilon, 0.25 + 2*epsilon)    # dense
  )
  
  for (layer in 1: L) {
    
    p1_vals <- sapply(prob_ranges, function(rng) runif(1, rng[1]/2, rng[2]/2)) # inter is lower
    p2_vals <- sapply(prob_ranges, function(rng) runif(1, rng[1], rng[2]))
    
    # Generate the blockwise adjacency matrices for each set of probabilities
    P = get_blockwise_const_mat(n, n_c[1], p1_vals[1], p2_vals[1])
    probability_1[, , layer] = P
    
    P = get_blockwise_const_mat(n, n_c[2], p1_vals[2], p2_vals[2])
    probability_2[, , layer] = P
    
    P = get_blockwise_const_mat(n, n_c[3], p1_vals[3], p2_vals[3])
    probability_3[, , layer] = P
  }
  
  return(list(probability_1, probability_2, probability_3))
}





generate_tensor_probability_directed <- function(n_1, n_2, L, probability){
  # generate adjacency matrix in {0,1}
  
  dim_ = c(n_1, n_2, L)
  A = array(NA,dim_)
  
  for (layer in 1: L)
  {
    A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
  }
  
  return(A)
}



sim_SBM_array <- function(num_seq=1, n=50, rho=0.5, L=4, true_CP=c(25,50,75), num_time=100){
  
  # This function directly output the adjacency matrix with size (num_seq, num_time, n, n, L)
  # rho controls the temporal dependency
  
  A.all_seq <- array(NA, c(num_seq, num_time, n, n, L)) # i.e. 10 sequences empty
  
  for(seq_iter in 1:num_seq){
    
    
    A.tensor <- array(NA, c(num_time, n, n, L)) # 1 sequence
    
    K = 3
    v = true_CP
    
    
    for (layer in 1: L){ ### FOR EACH LAYER, THEN FOR EACH TIME T
      
      for(t in 1:num_time){
        
        if( t==1 || t==(v[2]+1) ){ # t=1 or t=51
          
          P =  matrix(0.3,n,n)
          P[1:floor(n/K), 1:floor(n/K)] = 0.5
          P[(1+floor(n/K)):(2*floor(n/K)),(1+floor(n/K)):(2*floor(n/K)) ] = 0.5
          P[(1+2*floor(n/K)):n,(1+2*floor(n/K)):n ] = 0.5
          diag(P) = 0
          
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
          
        }
        
        if(t == (v[1]+1) || t == (v[3]+1)){ # t=26 or t=76
          
          Q =  matrix(0.2,n,n)
          Q[1:floor(n/K), 1:floor(n/K)] = 0.45
          Q[(1+floor(n/K)):(2*floor(n/K)),(1+floor(n/K)):(2*floor(n/K)) ] = 0.45
          Q[(1+2*floor(n/K)):n,(1+2*floor(n/K)):n ] = 0.45
          diag(Q) = 0
          
          A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),Q),n,n)
          
        }
        
        if( (t > 1 && t <= v[1])  ||  (t > v[2] && t <= v[3]) ){ # t=2 to t=25 or t=51 to t=75
          
          aux1 = (1-P)*rho + P # (1-E(t+1))*rho + E(t+1) if A(t) = 1
          aux2 = P*(1-rho) # (1-E(t+1))*rho + E(t+1) if A(t) = 0
          
          aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
          aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
          A =  aux1*A + aux2*(1-A)
          
        }
        
        if( (t > v[1] && t <= v[2]) || ((t > v[3] && t <= num_time)) ){ # t=27 to t=50 or t=77 to t=100
          
          aux1 = (1-Q)*rho + Q
          aux2 = Q*(1-rho)
          
          aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
          aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
          A =  aux1*A + aux2*(1-A)
          
        }
        
        diag(A) <- 0
        A.tensor[t,,,layer] = A
        
      }
      
    }
    
    
    A.all_seq[seq_iter,,,,] <- A.tensor
    
  }
  
  return(A.all_seq)
}



