
library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("utility.R")



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
                             threshold = c(1800, 600, 200, 50), method = "Narrowest", obj.B = B.tensor.odd)
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = c(1000, 500, 250, 50), method = "Greedy", obj.B = B.tensor.odd)

results[[2]]


###########################
# Example Model Selection # 
###########################
source("SBS.R")
init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1)
threshold_list <- seq(1, max(init[[1]]$results[, 2]), length.out=25)
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)

# Manual calculation
for (i in seq_along(threshold_list)) {
  candidates <- 2*sort(results[[i+1]]$results[, 1])
  BIC <- cal_BIC(A.tensor, candidates, hat.rank)
  cat("Candidates: ", paste(candidates, collapse = ", "), ". BIC = ", BIC, "\n", sep = "")
}

# using utilities/model selection
model_selection(results, A.tensor, hat.rank = hat.rank)

# We can pass in other function to determine selections statistics 
# model_selection(results, A.tensor, method = cal_BIC, hat.rank = hat.rank)



