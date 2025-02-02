
library(rTensor)
source("SBS.R")
source("CUSUM.R")

##################
# simulation SBM #
##################

get_blockwise_const_mat <- function(n, n_c, p_1, p_2){
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


get_sbm_params <- function(n, L, n_c=c(4, 4), flip_layer=TRUE){
  
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




#### sbm directed
generate_tensor_probability_directed <- function(n_1, n_2, L, probability){
  
  dim_ = c(n_1, n_2, L)
  A = array(NA,dim_)
  
  for (layer in 1: L)
  {
    A[, , layer] = matrix(rbinom(matrix(1,n_1,n_2),matrix(1,n_1,n_2), probability[ , , layer]),n_1,n_2)
  }
  
  return(A)
}



###################
# generating data #
###################

A.tensor <- array(NA, c(150, 50, 50, 4))
sbm_params <- get_sbm_params(n=50, L=4, n_c=c(4, 4), flip_layer=TRUE)

probability_1 = sbm_params[[1]]
probability_2 = sbm_params[[2]]

for(t_iter in 1:50) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
for(t_iter in 51:100) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_2)
for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)

# for(t_iter in 1:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
# for(t_iter in 81:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_2)



##################################
# Example Code (Not CUSUM_step1) # 
##################################

hat.rank <- c(15, 15, 15)
s <- 0
e <- 150
frobenius_holder <- numeric(149)
source("CUSUM.R")

for(t in 2:148){
  frobenius_holder[t] <- CUSUM_frob_SBS(A.tensor, s, e, t, hat.rank)
}

plot(1:149, frobenius_holder, type='l')

intervals <- construct_intervals(150, 1/2, 4)
results_one <- cusum_on_intervals(CUSUM_frobenius, A.tensor, c(30, 38), rank = hat.rank)
results_all_Frobenius <- cusum_on_intervals(CUSUM_frobenius, A.tensor, intervals, rank = hat.rank)
results_all_Layer <- cusum_on_intervals(CUSUM_layer, A.tensor, intervals, rank = hat.rank)

# Calculate intervals & CUSUM within
# results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, 150, alpha = 1/2, m = 4, 
#                              threshold = 18, method = "Narrowest", rank = hat.rank)

# Pass in CUSUM results (as a matrix) to save computation
results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, 150, CUSUM_res = results_all_Layer, 
                             threshold = c(10, 9, 8), method = "Greedy", rank = hat.rank)

results[[2]]$results

##############################
# Example Code (CUSUM_step1) # 
##############################

A.tensor.even <- A.tensor[seq(2, 150, by = 2), , , ]
B.tensor.odd  <- A.tensor[seq(1, 149, by = 2), , , ]
if (dim(A.tensor.even)[1] != dim(B.tensor.odd)[1]) {
  stop("Make sure even and odd have same length")
}

hat.rank <- c(15, 15, 15)
s <- 0
e <- 75
frobenius_holder <- numeric(74)

for(t in 2:73){
  frobenius_holder[t] <- CUSUM_step1(A.tensor.even, s, e, t, obj.B = B.tensor.odd)
}

plot(1:74, frobenius_holder, type='l')

intervals <- construct_intervals(75, sqrt(1/2), 2)
results_one <- cusum_on_intervals(CUSUM_step1, A.tensor.even, c(30, 38), obj.B = B.tensor.odd)
results_all_step1 <- cusum_on_intervals(CUSUM_step1, A.tensor.even, intervals, obj.B = B.tensor.odd)


# Calculate intervals & CUSUM within
# results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, alpha = 1/2, m = 4,
#                              threshold = 200, method = "Narrowest", obj.B = B.tensor.odd)

# Pass in CUSUM results (as a matrix) to save computation
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = c(1000, 500, 250, 50), method = "Narrowest", obj.B = B.tensor.odd)
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = c(1000, 500, 250, 50), method = "Greedy", obj.B = B.tensor.odd)

results[[2]]$results



# Candidate Selection: 
# Narrowest is statistically (slightly) preferred, 
# but is computationally more expensive for multiple thresholds
# Greedy is most efficient for multiple thresholds

# We are sensitive to the endpoints for Frobenius CUSUMs
# Fan's weighting seems to fix endpoint issue

# TODO Model selection