source("simulate_data.R")
source("CUSUM.R")
source("utility.R")


A.tensor <- array(NA, c(150, 50, 50, 4))
sbm_params <- get_sbm_params(n=50, L=4, n_c=c(4, 4), flip_layer=TRUE)

probability_1 = sbm_params[[1]]
probability_2 = sbm_params[[2]]

for(t_iter in 1:50) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
for(t_iter in 51:100) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_2)
for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)


cal_BIC(A.tensor, c(50,100), hat.rank = c(15,15,15))
cal_BIC(A.tensor, c(25, 50, 100), hat.rank = c(15,15,15))
cal_BIC(A.tensor, c(25, 50, 75, 100,125), hat.rank = c(15,15,15))
cal_BIC(A.tensor, c(100), hat.rank = c(15,15,15)) 

# the log_lik for 1 CP is similar to the log_lik for 2 CP
# but the num_par for 1 CP is substantially less than the num_par for 2 CP
# so the BIC for 1 CP is small than the BIC for 2 CP (2 CP is ground truth)
