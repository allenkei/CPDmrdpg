
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


