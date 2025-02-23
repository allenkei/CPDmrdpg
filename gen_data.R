source("simulate_data.R")

set.seed(123)


scenario <- "s2" 
num_node <- 50





# FIXED
if(scenario == "s1"){
  
  # LAYER FLIPPED
  num_seq <- 10
  num_time <- 150
  num_layer <- 4
  num_block_before <- 4 # same
  num_block_after <- 4 # same 
  FL = TRUE
}else if(scenario == "s2"){
  
  # BLOCK NUMBER CHANGED
  num_seq <- 10
  num_time <- 150
  num_layer <- 4
  num_block_before <- 3 # different
  num_block_after <- 4 # different
  FL = FALSE
}


A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty


sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before, num_block_after), flip_layer=FL)
probability_1 = sbm_params[[1]]
probability_2 = sbm_params[[2]]

# begin simulate data
for(seq_iter in 1:num_seq){
  
  A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
  
  # T from 1 to 150 (otherwise change the for loop)
  for(t_iter in 1:50) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
  for(t_iter in 51:100) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
  for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
  
  A.all_seq[seq_iter,,,,] <- A.tensor
}; rm(seq_iter, A.tensor)


dim(A.all_seq) 
save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists





