source("simulate_data.R")

set.seed(123)


scenario <- "s3" 
num_node <- 50


#######################
# FIXED CONFIGURATION #
#######################

if(scenario == "s1"){
  
  # LAYER FLIPPED
  num_seq <- 10
  num_time <- 150
  num_layer <- 4
  num_block_before <- 4 # same
  num_block_after <- 4 # same 
  FL = TRUE # FLIP
}else if(scenario == "s2"){
  
  # BLOCK NUMBER CHANGED
  num_seq <- 10
  num_time <- 150
  num_layer <- 4
  num_block_before <- 3 # different
  num_block_after <- 4 # different
  FL = FALSE # NO FLIP
}else if(scenario == "s3"){
  
  # COMBINE s1 and s2 WITH LONGER TIME SPAN
  num_seq <- 10
  num_time <- 400 # LONGER
  num_layer <- 4
  
  # BLOCK NUMBER K CHANGED
  num_block_before_K <- 4 # different
  num_block_after_K <- 3 # different
  FL_K = FALSE # same
  
  # LAYER L FLIPPED
  num_block_before_L <- 4 # same
  num_block_after_L <- 4 # same 
  FL_L = TRUE # different
}

#############
# SAVE DATA #
#############

if(scenario == "s1" | scenario == "s2"){
  
  sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before, num_block_after), flip_layer=FL)
  probability_1 = sbm_params[[1]]
  probability_2 = sbm_params[[2]]
  
  A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
  
  # begin simulate data
  for(seq_iter in 1:num_seq){
    
    A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
    
    # T from 1 to 150 (otherwise change the for loop)
    for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
    for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
    for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
  
    A.all_seq[seq_iter,,,,] <- A.tensor
  }; rm(seq_iter, A.tensor)
  
  
  dim(A.all_seq) 
  save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
  
}else if(scenario == "s3"){
  
  # BLOCK
  sbm_params_K <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_K, num_block_after_K), FL_K)
  probability_K_1 = sbm_params_K[[1]]
  probability_K_2 = sbm_params_K[[2]]
  
  # LAYER
  sbm_params_L <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_L, num_block_after_L), FL_L)
  probability_L_1 = sbm_params_L[[1]]
  probability_L_2 = sbm_params_L[[2]]
  
  A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
  
  # begin simulate data
  for(seq_iter in 1:num_seq){
    
    A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
    
    # T from 1 to 150 (otherwise change the for loop)
    for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
    for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
    for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
    for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
    for(t_iter in 201:250) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_1)
    for(t_iter in 251:300) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
    for(t_iter in 301:350) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
    for(t_iter in 351:400) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
    
    A.all_seq[seq_iter,,,,] <- A.tensor
  }; rm(seq_iter, A.tensor)
  
  
  dim(A.all_seq) 
  save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
  
}

