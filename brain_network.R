
library(rTensor)
source("SBS.R")

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


##########################
# Change Point Detection #
##########################





diff_frobenius <- function(A, B) sum((A - B)^2)^0.5



Hetero_PCA_test <- function(Y, r, tmax = 20, vartol = 1e-6){
  
  N_t = Y
  r = min(c(r, dim(N_t)))
  diag(N_t) = 0
  U_t = matrix(NA, nrow = dim(Y)[1], r)
  t = 1
  approx = -1
  
  while(t <= tmax){ # Stop criterion: convergence or maximum number of iteration reached
    temp = svd(N_t)
    U_t = temp$u[,1:r]
    V_t = temp$v[,1:r]
    if (r > 1){
      tilde_N_test_t = U_t %*% diag(temp$d[1:r]) %*% t(V_t)
    }
    else{
      tilde_N_test_t = temp$d[1] * U_t %*% t(V_t)
    }
    
    N_test_new = N_t
    diag(N_test_new) = diag(tilde_N_test_t)
    N_t = N_test_new
    svector = diag(tilde_N_test_t)
    if (abs(sum(svector^2) - approx) > vartol){
      t = t+1
      approx = sum(svector^2)
    }
    else {
      break
    }
  }
  return(U_t)
}



Tensor_Hetero_PCA_test <- function(Y, r, tmax = 20){
  
  p = dim(Y)
  d = length(p)
  U_0 = list()
  for (i in 1:d){
    MY = k_unfold(Y, i)@data
    MY_Hetero_PCA = Hetero_PCA_test(MY %*% t(MY), r[i], tmax)
    U_0 = c(U_0, list(MY_Hetero_PCA))
  }
  return(U_0)
}


estimate_thpca <- function(Y.tensor, hat.rank, tmax = 20){
  U.hat = Tensor_Hetero_PCA_test(Y.tensor, hat.rank, tmax)
  P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
  P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
  P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
  Y.hat = ttm(ttm(ttm(Y.tensor, P.U1, 1), P.U2, 2), P.U3, 3) 
  
  P_hat  = Y.hat@data
  
  P_hat[P_hat > 1]  = 1
  P_hat[P_hat < 0]  = 0
  
  return(P_hat)
}


CUSUM_frobenius <- function(obj, s, e, t, rank) {
  print(paste0("s = ", s, ", e = ", e, ", t = ", t, "."))
  if ((t - s) == 1) {
    sum_s_t <- as.tensor(A.tensor[(s+1):t, , , ])
  } else {
    sum_s_t  <- (1/(t - s)) * as.tensor( apply(A.tensor[(s+1):t, , , ], c(2, 3, 4), sum) )
  }
  
  if ((e - t )== 1) {
    sum_t_e <- as.tensor(A.tensor[(t+1):e, , , ])
  } else {
    sum_t_e <- (1/(e - (t))) * as.tensor( apply(A.tensor[(t+1):e, , , ], c(2, 3, 4), sum) )
  }
  
  P_s_t  <- estimate_thpca(sum_s_t, rank, tmax = 20)
  P_t_e <- estimate_thpca(sum_t_e, rank, tmax = 20)
  
  return(diff_frobenius(P_s_t, P_t_e))
}

CUSUM_layer <- function(obj, s, e, t, rank) {
  # Implement CUSUM layer
}


hat.rank <- c(15, 15, 15)
s <- 0
e <- 150
frobenius_holder <- numeric(149)

for(t in 2:148){
  # print(paste0("s = ", s, ", e = ", e, ", t = ", t, "."))
  # if ((t - s) == 1) {
  #   sum_s_t <- as.tensor(A.tensor[(s+1):t, , , ])
  # } else {
  #   sum_s_t  <- (1/(t - s)) * as.tensor( apply(A.tensor[(s+1):t, , , ], c(2, 3, 4), sum) )
  # }
  # 
  # if ((e - t )== 1) {
  #   sum_t_e <- as.tensor(A.tensor[(t+1):e, , , ])
  # } else {
  #   sum_t_e <- (1/(e - (t))) * as.tensor( apply(A.tensor[(t+1):e, , , ], c(2, 3, 4), sum) )
  # }
  # 
  # P_s_t  <- estimate_thpca(sum_s_t, hat.rank, tmax = 20)
  # P_t_e <- estimate_thpca(sum_t_e, hat.rank, tmax = 20)
  # 
  # frobenius_holder[t] <- (diff_frobenius(P_s_t, P_t_e))
  
  frobenius_holder[t] <- CUSUM_frobenius(A.tensor, s, e, t, hat.rank)
}

plot(1:149, frobenius_holder, type='l')

intervals <- construct_intervals(150, 1/2, 4)
results_one <- cusum_on_intervals(CUSUM_frobenius, A.tensor, c(30, 38), rank = hat.rank)
results_all <- cusum_on_intervals(CUSUM_frobenius, A.tensor, intervals, rank = hat.rank)

# Calculate intervals & CUSUM within
# results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, 150, alpha = 1/2, m = 4, 
#                              threshold = 18, method = "Narrowest", rank = hat.rank)

# Pass in CUSUM results (as a matrix) to save computation
results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, 150, CUSUM_res = results_all, 
                             threshold = c(22, 18, 16), method = "Greedy", rank = hat.rank)

results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, 150, CUSUM_res = results_all, 
                             threshold = c(22, 18, 16), method = "Narrowest", rank = hat.rank)

results[[1]]$results

# Candidate Selection: 
# Narrowest is statistically (slightly) preferred, 
# but is computationally more expensive for multiple thresholds
# Greedy is most efficient for multiple thresholds

# We are sensitive to the endpoints
# Do we increase m?
# Or try other statistics first

# Check on m doing anticipated behavior 

# Check on weights for CUSUM / SBS version of CUSUM