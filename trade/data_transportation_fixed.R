
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




#### data preprocess ####


path = "air_transportation/"

data_whole_22 = read.csv(paste(c(path, "Air Transportation_2022.csv"), collapse =  ""), header = T) 
data_whole_22$MONTH  = 202200 + data_whole_22$MONTH 
data_whole_21 = read.csv(paste(c(path, "Air Transportation_2021.csv"), collapse =  ""), header = T) 
data_whole_21$MONTH  = 202100 + data_whole_21$MONTH
data_whole_20 = read.csv(paste(c(path, "Air Transportation_2020.csv"), collapse =  ""), header = T) 
data_whole_20$MONTH  = 202000 + data_whole_20$MONTH
data_whole_19 = read.csv(paste(c(path, "Air Transportation_2019.csv"), collapse =  ""), header = T) 
data_whole_19$MONTH  = 201900 + data_whole_19$MONTH
data_whole_18 = read.csv(paste(c(path, "Air Transportation_2018.csv"), collapse =  ""), header = T) 
data_whole_18$MONTH  = 201800 + data_whole_18$MONTH
data_whole_17 = read.csv(paste(c(path, "Air Transportation_2017.csv"), collapse =  ""), header = T) 
data_whole_17$MONTH  = 201700 + data_whole_17$MONTH
data_whole_16 = read.csv(paste(c(path, "Air Transportation_2016.csv"), collapse =  ""), header = T) 
data_whole_16$MONTH  = 201600 + data_whole_16$MONTH
data_whole_15 = read.csv(paste(c(path, "Air Transportation_2015.csv"), collapse =  ""), header = T) 
data_whole_15$MONTH  = 201500 + data_whole_15$MONTH

data_whole = rbind(data_whole_22,data_whole_21,  data_whole_20, data_whole_19, data_whole_18, data_whole_17, data_whole_16, data_whole_15)


colnames(data_whole)

AIRLINE_sort = rle(sort(data_whole[,1]))
AIRLINE = AIRLINE_sort$value[order(AIRLINE_sort $lengths, decreasing = T)[1:4]]

#19393 "Southwest Airlines Co.: WN"
#20304  "SkyWest Airlines Inc.: OO"
#19977  "United Air Lines Inc.: UA"
#19790  "Delta Air Lines Inc.: DL"
#20107  "Federal Express Corporation: FX"

data_airport_id = read.csv(paste(c(path, "L_AIRPORT_ID.csv"), collapse =  ""), header = T) 


ORIGIN_AIRPORT = rle(sort(data_whole[,2]))
og_id_100 = ORIGIN_AIRPORT$value[order(ORIGIN_AIRPORT$lengths, decreasing = T)[1:50]]

ORIGIN_AIRPORT_names = data_airport_id$Description[data_airport_id$Code %in% og_id_100]


DEST_AIRPORT = rle(sort(data_whole[,4]))
de_id_100 = DEST_AIRPORT$value[order(DEST_AIRPORT$lengths, decreasing = T)[1:50]]
DEST_AIRPORT_names = data_airport_id$Description[data_airport_id$Code %in% de_id_100]

og_id_100_nan = data_airport_id$Description[data_airport_id$Code %in%  og_id_100[which(! og_id_100 %in% de_id_100)]]

de_id_100_nan = data_airport_id$Description[data_airport_id$Code %in%  de_id_100[which(! de_id_100 %in% og_id_100)]]


month_time = sort(unique(data_whole$MONTH))

#node og_id_100
#node de_id_100
#layer AIRLINE
#time month_time

n_1 = length(og_id_100)
n_2 = length(de_id_100)
L =  length(AIRLINE)
TT = length(month_time)
dim = c(n_1, n_2, L, TT)
data_month_product_tensor = array(rep(NA, n_1*n_2*L*TT), dim)



for (u in 1:TT){
  
  month = month_time[u]
  
  data_month = data_whole[which(data_whole$MONTH == month),]
  
  for(l in 1:L){
    
    data_month_airline = data_month[which(data_month$AIRLINE_ID == AIRLINE[l]), ]
    
  
    for(i in 1: n_1){
      for(j in 1: n_2){
        x = which( data_month_airline$ORIGIN_AIRPORT_ID == og_id_100[i] & 
                     data_month_airline$DEST_AIRPORT_ID == de_id_100[j] )
        if ( length(x)!=0 ){
          data_month_product_tensor[i, j, l, u] = 1
        }
        else{
          data_month_product_tensor[i, j, l, u] = 0
        }
     }
   }

 }
}






#### burn-in ####

B = 100

K_max = 500

T_burn = 30

TT_list = TT  - T_burn


h_kernel = (K_max*log(TT_list*n_1 * n_2)/(TT_list*n_1*n_2))^{1/L}
h_kernel 

hat.rank = rep(10, 3)


max_D_rescale_B = rep(0, B)
max_D_rescale_B_thpca = rep(0, B)
max_D_rescale_B_uase = rep(0, B)
max_D_rescale_B_multi = rep(0, B)

tic()

set.seed(0) 

for(b in 1:B){
  b_ix = sample(1:T_burn, replace = FALSE)
  A_b = data_month_product_tensor[, , , b_ix]
  
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




alpha = 0.05

tau_factor = quantile(max_D_rescale_B, 1 - alpha, type = 1)
tau_factor_thpca = quantile(max_D_rescale_B_thpca, 1 - alpha, type = 1)
tau_factor_uase = quantile(max_D_rescale_B_uase, 1 - alpha, type = 1)
tau_factor_multi = quantile(max_D_rescale_B_multi, 1 - alpha, type = 1)



#### online cpd ####


A_list  =  data_month_product_tensor[, , , (T_burn+1) : TT]

A_list_list = list()
for (t in 1:TT_list){
  A_list_list[[t]] = A_list[, , , t]
}


result_online_cpd = online_cpd_fixed_cpp(A_list_list, tau_factor, hat.rank, verbose = FALSE)
result_online_cpd_thpca = online_cpd_fixed_thpca_cpp(A_list_list, tau_factor_thpca, hat.rank, verbose = FALSE)
result_online_cpd_uase = online_cpd_fixed_uase_cpp(A_list_list, tau_factor_uase, hat.rank[1], verbose = FALSE)
result_online_cpd_multi = online_cpd_fixed_multi(A_list, tau_factor_multi, hat.rank[1], verbose = FALSE)


k_nn = 3
n0 = floor(0.3 * T_burn)
n1 = floor(0.7 * T_burn)
ARL = 5*10^7

A_all = data_month_product_tensor
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

month_time[T_burn + result_online_cpd$t]
month_time[T_burn + result_online_cpd_thpca$t]
month_time[T_burn + result_online_cpd_uase$t]
month_time[T_burn + result_online_cpd_multi$t]
month_time[T_burn +min(t_hat)]







