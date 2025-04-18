#### local ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("../generate_data.R")
source("../simulation_wrap.R")
Rcpp::sourceCpp("../cpd_hpca.cpp")
Rcpp::sourceCpp("../cpd_uase.cpp")
source("../other_methods.R")


#### main ####

library(mvtnorm)
library(dirmult)
library(rTensor)
library(Rcpp)
library(tictoc)

library(multiness)
library(gStream)

library(egg)
library(ggplot2)


load("trade_tensor.RData")


#### binarize tensors
cut_off = 0.5
time = dim(data_year_product_tensor)[4]


for (t in 1:time){
  data_year = data_year_product_tensor[, , , t]
  t_cut_off = quantile(as.numeric(data_year), cut_off)
  data_year[data_year <= t_cut_off] = 0
  data_year[data_year > t_cut_off] = 1
  data_year_product_tensor[, , , t] = data_year
}


### cpd

B = 100

K_max = 500

n_1 = dim(data_year_product_tensor)[1]

n_2 = dim(data_year_product_tensor)[2]

L = dim(data_year_product_tensor)[3]



T_burn = 10

TT = time - T_burn

h_kernel = (K_max*log(TT*n_1* n_2)/(TT*n_1*n_2))^{1/L}
h_kernel

hat.rank = c(10, 10, 10)
max_D_rescale_B = rep(0, B)
max_D_rescale_B_thpca = rep(0, B)
max_D_rescale_B_uase = rep(0, B)
max_D_rescale_B_multi = rep(0, B)
alpha = 0.05

tic()




set.seed(0) 

for(b in 1:B){
  b_ix = sample(1:T_burn, replace = FALSE)
  A_b =  data_year_product_tensor[, , , b_ix]
  
  A_b_list = list()
  for (t in 1:T_burn){
    A_b_list[[t]] = A_b[, , , t]
  }
  
  max_D_rescale_B[b] = max(max_D_s_t_rescale_fixed_cpp(A_b_list, hat.rank, verbose = FALSE))
  max_D_rescale_B_thpca[b] = max(max_D_s_t_rescale_fixed_thpca_cpp(A_b_list, hat.rank, verbose = FALSE))
  max_D_rescale_B_uase[b] = max(max_D_s_t_rescale_fixed_uase_cpp(A_b_list, hat.rank[1], verbose = FALSE))
  max_D_rescale_B_multi[b] = max(max_D_s_t_rescale_fixed_multi(A_b, hat.rank[1], verbose = FALSE))
  print(paste0("b = ", b))
}
toc()


tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
tau_factor_thpca = quantile(max_D_rescale_B_thpca, 1 - alpha, type = 1)
tau_factor_uase = quantile(max_D_rescale_B_uase, 1 - alpha, type = 1)
tau_factor_multi = quantile(max_D_rescale_B_multi, 1 - alpha, type = 1)


A_list  =  data_year_product_tensor[, , , (T_burn+1) : time]

A_list_list = list()
for (t in 1:TT){
  A_list_list[[t]] = A_list[, , , t]
}



result_online_cpd = online_cpd_fixed_cpp(A_list_list, tau_factor, hat.rank, verbose = FALSE)
result_online_cpd_thpca = online_cpd_fixed_thpca_cpp(A_list_list, tau_factor_thpca, hat.rank, verbose = FALSE)
result_online_cpd_uase = online_cpd_fixed_uase_cpp(A_list_list, tau_factor_uase, hat.rank[1], verbose = FALSE)
result_online_cpd_multi = online_cpd_fixed_multi(A_list, tau_factor_multi, hat.rank[1], verbose = FALSE)



tau_factor


k_nn = 3
n0 = floor(0.3 * T_burn)
n1 = floor(0.7 * T_burn)
ARL = 2000

A_all = data_year_product_tensor
n_all = dim(A_all)[4]

dist_mat_all = get_dist_mat(A_all, distance_Y)
diag(dist_mat_all) = max(dist_mat_all) + 100

t_hat = tryCatch(
  expr = {
    res = gstream(dist_mat_all, L = T_burn, N0 = T_burn, k_nn, statistics = "m",
                  n0, n1, ARL, alpha, skew.corr = TRUE, asymp = FALSE)
    tauhat = res$tauhat$max.type
    tauhat
  },
  error = function(cond){return (Inf)}
)
t_hat
length(t_hat)
min(t_hat)

data_names[T_burn + result_online_cpd$t]
data_names[T_burn + result_online_cpd_thpca$t]
data_names[T_burn + result_online_cpd_uase$t]
data_names[T_burn + result_online_cpd_multi$t]
data_names[T_burn +min(t_hat)]

