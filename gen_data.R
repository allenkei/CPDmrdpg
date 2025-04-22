source("utility.R")



# Kyles Generate function
# f1: No change, T = 200
# f2: Flip layers and change block numbers, T = 200 (5 change points) 
# f3: Block size change, first layer only, change to T = 200 (3 change points 1-2-3-1)
# f4: Sparsity fluctuates, T = 200 (5 change points 1-2-3-2-1-2)
# f5: Sparsity, weak difference, T = 200 (3 change points 1-2-3-1)
# f6: Directed Dirichlet, T = 150 (3 change points)


generate <- function(scenario, num_node = 50, num_seq = 1, save = FALSE) {
  
  num_time <- 200
  num_layer <- 4
  
  if(scenario == "f1"){
    # Old Scenario 4, extended
    num_block_before <- 4 # same
    num_block_after <- 4 # same
    FL = FALSE # NO FLIP
    
    sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before, num_block_after), flip_layer=FL)
    probability_1 = sbm_params[[1]] # ONLY USE THIS
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      # T from 1 to 150 (otherwise change the for loop)
      for(t_iter in 1:200) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  } else if(scenario == "f2"){
    # Old Scenario 3b
    # COMBINE s1 and s2 WITH LONGER TIME SPAN
    
    # BLOCK NUMBER K CHANGED
    num_block_before_K <- 4 # different
    num_block_after_K <- 3 # different
    FL_K = FALSE # same
    
    # LAYER L FLIPPED
    num_block_before_L <- 4 # same
    num_block_after_L <- 4 # same 
    FL_L = TRUE # different
    
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
      
      for(t_iter in 1:20) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
      for(t_iter in 21:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
      for(t_iter in 61:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
      for(t_iter in 81:160) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_1)
      for(t_iter in 161:180) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
      for(t_iter in 181:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
      
      # c(20, 60, 80, 160, 180)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  }else if(scenario == "f3"){
    # Old Scenario 6, extended
    block_size1 <- floor(c(3, 4, 3) / 10 * num_node) # fixed ratio
    block_size2 <- floor(c(4, 3, 3) / 10 * num_node) # fixed ratio
    block_size3 <- floor(c(5, 3, 2) / 10 * num_node) # fixed ratio
    
    sbm_params1 <- get_sbm_VS_FL_params(n=num_node, L=num_layer, block_size1, block_size2)
    sbm_params2 <- get_sbm_VS_FL_params(n=num_node, L=num_layer, block_size2, block_size3)
    probability_1 = sbm_params1[[1]]
    probability_2 = sbm_params1[[2]]
    probability_3 = sbm_params2[[2]] 
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      # T from 1 to 150 (otherwise change the for loop)
      for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
      for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  } else if(scenario == "f4"){
    # Old Scenario 8b
    epsilon <- 0.1
    
    sbm_params <- get_sbm_params_spa_inc(n=num_node, L=num_layer, epsilon = epsilon)
    probability_1 = sbm_params[[1]]
    probability_2 = sbm_params[[2]]
    probability_3 = sbm_params[[3]]
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      # T from 1 to 150 (otherwise change the for loop)
      for(t_iter in 1:20) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 21:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 61:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
      for(t_iter in 81:160) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 161:180) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 181:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      
      # c(20, 60, 80, 160, 180)
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
    
  } else if(scenario == "f5") {
    # Old Scenario 8, extended
    epsilon <- 0.05
    
    sbm_params <- get_sbm_params_spa_inc(n=num_node, L=num_layer, epsilon = epsilon)
    probability_1 = sbm_params[[1]]
    probability_2 = sbm_params[[2]]
    probability_3 = sbm_params[[3]]
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      
      A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
      
      # T from 1 to 150 (otherwise change the for loop)
      for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
      for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
      for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
      
      A.all_seq[seq_iter,,,,] <- A.tensor
    }; rm(seq_iter, A.tensor)
  } else if(scenario == "f6"){
    
    cp_truth <- c(70, 140)
    d <- 3
    
    A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
    
    # begin simulate data
    for(seq_iter in 1:num_seq){
      params <- get_dirichlet_params(num_node, num_node, num_layer, d)
      dirichlet_xy_1 <- params[[1]]; dirichlet_xy_2 <- params[[2]]
      W_1 <- params[[3]]; W_2 <- params[[4]]
      
      A.all_seq[seq_iter,,,,] <- get_data_cp_dirichlet(num_time, cp_truth, num_node, num_node, 
                                                       num_layer, dirichlet_xy_1, W_1, dirichlet_xy_2, W_2, directed = TRUE)
    }
  }
  
  if (save == TRUE) {save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData"))}
  return(A.all_seq)
}



















#set.seed(123)
#scenario <- "s8b" ### "s1","s2","s3b","s4","s5","s6b","s7"
#num_node <- 50 ### 50, 100
#num_seq <- 10 ### 50, 100 # 10 is for testing the code



#######################
# FIXED CONFIGURATION #
#######################

# if(scenario == "s1"){
#   
#   # LAYER FLIPPED
#   num_time <- 150
#   num_layer <- 4
#   num_block_before <- 4 # same
#   num_block_after <- 4 # same 
#   FL = TRUE # FLIP
#   
# }else if(scenario == "s2"){
#   
#   # BLOCK NUMBER CHANGED
#   num_time <- 150
#   num_layer <- 4
#   num_block_before <- 3 # different
#   num_block_after <- 4 # different
#   FL = FALSE # NO FLIP
#   
# }else if(scenario == "s3"){
#   
#   # COMBINE s1 and s2 WITH LONGER TIME SPAN
#   num_time <- 300 # longer time span
#   num_layer <- 4
#   
#   # BLOCK NUMBER K CHANGED
#   num_block_before_K <- 4 # different
#   num_block_after_K <- 3 # different
#   FL_K = FALSE # same
#   
#   # LAYER L FLIPPED
#   num_block_before_L <- 4 # same
#   num_block_after_L <- 4 # same 
#   FL_L = TRUE # different
#   
# }else if(scenario == "s3b"){
#   
#   # COMBINE s1 and s2 WITH LONGER TIME SPAN
#   num_time <- 200 # longer time span
#   num_layer <- 4
#   
#   # BLOCK NUMBER K CHANGED
#   num_block_before_K <- 4 # different
#   num_block_after_K <- 3 # different
#   FL_K = FALSE # same
#   
#   # LAYER L FLIPPED
#   num_block_before_L <- 4 # same
#   num_block_after_L <- 4 # same 
#   FL_L = TRUE # different
#   
# }else if(scenario == "s4"){
#   
#   # NO CHANGED
#   num_time <- 150
#   num_layer <- 4
#   num_block_before <- 4 # same
#   num_block_after <- 4 # same
#   FL = FALSE # NO FLIP
#   
# }else if(scenario == "s5"){
#   
#   # VARIED BLOCK SIZES
#   num_time <- 150
#   num_layer <- 4
#   block_size1 <- floor(c(3, 4, 3) / 10 * num_node) # fixed ratio
#   block_size2 <- floor(c(2, 3, 5) / 10 * num_node) # fixed ratio
#   block_size3 <- floor(c(4, 3, 3) / 10 * num_node) # fixed ratio
#   
# }else if(scenario == "s6"){
#   
#   # FIRST LAYER DIFFERS
#   num_time <- 150
#   num_layer <- 4
#   block_size1 <- floor(c(3, 4, 3) / 10 * num_node) # fixed ratio
#   block_size2 <- floor(c(4, 3, 3) / 10 * num_node) # fixed ratio
#   
# }else if(scenario == "s6b"){
#   
#   # FIRST LAYER DIFFERS
#   num_time <- 150
#   num_layer <- 4
#   block_size1 <- floor(c(3, 4, 3) / 10 * num_node) # fixed ratio
#   block_size2 <- floor(c(5, 2, 3) / 10 * num_node) # fixed ratio
#   
# }else if(scenario == "s6.06"){
#   # .02, .04, .06, .08
#   # FIRST LAYER DIFFERS
#   num_time <- 150
#   num_layer <- 4
#   epsilon <- 0.6
#   block_size1 <- floor(c(3, 3+epsilon, 3) / 10 * num_node) # fixed ratio
#   block_size2 <- floor(c(3+epsilon, 3, 3) / 10 * num_node) # fixed ratio
#   
# }else if(scenario == "s7"){
#   
#   # SEE 'utility.R' FOR DETAILS
#   # THIS SCENARIO HAS TEMPORAL DEPENDENCY
#   num_time <- 100
#   num_layer <- 4
#   true_CP <- c(20,50,80) # RIGHT NOW, need exactly 3 CPs; no need to be evenly-spaced; see codes
#   rho <- 0.1 # TEMPORAL DEPENDENCY; can try (0.1, 0.5, 0.9)
#   
# }else if(scenario == "s7b"){
#   
#   # SEE 'utility.R' FOR DETAILS
#   # THIS SCENARIO HAS TEMPORAL DEPENDENCY
#   num_time <- 100
#   num_layer <- 4
#   true_CP <- c(20,50,80) # RIGHT NOW, need exactly 3 CPs; no need to be evenly-spaced; see codes
#   rho <- 0.5
#   
# }else if(scenario == "s8"){
#   num_time <- 150
#   num_layer <- 4
#   epsilon <- 0.05 # Default 0.1
# }else if(scenario == "s8.05"){
#   num_time <- 150
#   num_layer <- 4
#   epsilon <- 0.05 
# }else if(scenario == "s8b"){
#   num_time <- 200
#   num_layer <- 4
#   epsilon <- 0.1 
# }
# 
# 






#############
# SAVE DATA #
#############

# if(scenario == "s1" | scenario == "s2"){
# if(scenario == "s1"){
#   
#   sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before, num_block_after), flip_layer=FL)
#   probability_1 = sbm_params[[1]]
#   probability_2 = sbm_params[[2]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:40) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     for(t_iter in 41:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 61:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#   
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# }else if(scenario == "s2"){
#   
#   #sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before, num_block_after), flip_layer=FL)
#   probability_1 = sbm_params[[1]]
#   probability_2 = sbm_params[[2]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# }else if(scenario == "s3"){
#   
#   # BLOCK
#   sbm_params_K <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_K, num_block_after_K), FL_K)
#   probability_K_1 = sbm_params_K[[1]]
#   probability_K_2 = sbm_params_K[[2]]
#   
#   # LAYER
#   sbm_params_L <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_L, num_block_after_L), FL_L)
#   probability_L_1 = sbm_params_L[[1]]
#   probability_L_2 = sbm_params_L[[2]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 300 (otherwise change the for loop)
#     for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
#     for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
#     for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
#     for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_1)
#     for(t_iter in 201:250) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
#     for(t_iter in 251:300) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
# 
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# }else if(scenario == "s3b"){
#   
#   # BLOCK
#   sbm_params_K <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_K, num_block_after_K), FL_K)
#   probability_K_1 = sbm_params_K[[1]]
#   probability_K_2 = sbm_params_K[[2]]
#   
#   # LAYER
#   sbm_params_L <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before_L, num_block_after_L), FL_L)
#   probability_L_1 = sbm_params_L[[1]]
#   probability_L_2 = sbm_params_L[[2]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 300 (otherwise change the for loop)
#     # for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
#     # for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
#     # for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
#     # for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_1)
#     # for(t_iter in 201:250) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
#     # for(t_iter in 251:300) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
#     
#     for(t_iter in 1:20) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
#     for(t_iter in 21:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
#     for(t_iter in 61:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_2)
#     for(t_iter in 81:160) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_L_1)
#     for(t_iter in 161:180) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_2)
#     for(t_iter in 181:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_K_1)
# 
#     # c(20, 60, 80, 160, 180)
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# }else if(scenario == "s4"){
#   
#   sbm_params <- get_sbm_params(n=num_node, L=num_layer, n_c=c(num_block_before, num_block_after), flip_layer=FL)
#   probability_1 = sbm_params[[1]] # ONLY USE THIS
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:150) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# }else if(scenario == "s5"){
#   
#   sbm_params1 <- get_sbm_var_size_params(n=num_node, L=num_layer, block_size1, block_size2)
#   sbm_params2 <- get_sbm_var_size_params(n=num_node, L=num_layer, block_size2, block_size3)
#   probability_1 = sbm_params1[[1]]
#   probability_2 = sbm_params1[[2]]
#   probability_3 = sbm_params2[[2]] # the second prob mat from sbm_params2
#   
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:40) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     for(t_iter in 41:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 61:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# } else if(startsWith(scenario, "s6")){
#   
#   sbm_params <- get_sbm_VS_FL_params(n=num_node, L=num_layer, block_size1, block_size2)
#   probability_1 = sbm_params[[1]]
#   probability_2 = sbm_params[[2]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# } else if(scenario == "s7" || scenario == "s7b"){
#   
#   # this function directly output A.all_seq with (num_seq, num_time, n, n, L)
#   A.all_seq <- sim_SBM_array(num_seq, num_node, rho, num_layer, true_CP, num_time)
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# } else if(scenario == "s8b"){
#   sbm_params <- get_sbm_params_spa_inc(n=num_node, L=num_layer, epsilon = epsilon)
#   probability_1 = sbm_params[[1]]
#   probability_2 = sbm_params[[2]]
#   probability_3 = sbm_params[[3]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:20) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     for(t_iter in 21:60) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 61:120) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
#     for(t_iter in 121:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 151:200) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
#   
# } else if(startsWith(scenario, "s8")) {
#   
#   sbm_params <- get_sbm_params_spa_inc(n=num_node, L=num_layer, epsilon = epsilon)
#   probability_1 = sbm_params[[1]]
#   probability_2 = sbm_params[[2]]
#   probability_3 = sbm_params[[3]]
#   
#   A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
#   
#   # begin simulate data
#   for(seq_iter in 1:num_seq){
#     
#     A.tensor <- array(NA, c(num_time, num_node, num_node, num_layer)) # 1 sequence
#     
#     # T from 1 to 150 (otherwise change the for loop)
#     for(t_iter in 1:50) A.tensor[t_iter,,,]    <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_1)
#     for(t_iter in 51:100) A.tensor[t_iter,,,]  <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_2)
#     for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=num_node, n_2=num_node, L=num_layer, probability_3)
#     
#     A.all_seq[seq_iter,,,,] <- A.tensor
#   }; rm(seq_iter, A.tensor)
#   
#   dim(A.all_seq) 
#   save(A.all_seq, file = paste0("data/seq",num_seq,"n",num_node,scenario,".RData")) # data folder exists
# }
  
  

 

