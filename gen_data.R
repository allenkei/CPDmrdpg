source("simulate_data.R")

set.seed(123)


scenario <- "s1" 
num_node <- 100

# FIXED
if(scenario == "s1"){
  num_seq <- 10
  num_time <- 150
  num_layer <- 4
  num_block <- 4
}


A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty


sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(4, 4), flip_layer=TRUE)
probability_1 = sbm_params[[1]]
probability_2 = sbm_params[[2]]

# begin simulate data
for(seq_iter in 1:num_seq){
  
  A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
  for(t_iter in 1:50) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
  for(t_iter in 51:100) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_2)
  for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
  
  A.all_seq[seq_iter,,,,] <- A.tensor
}; rm(seq_iter, A.tensor)

save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,"s1.RData")) # data folder exists





