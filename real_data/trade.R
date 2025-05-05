load("real_data/trade/trade_tensor.RData")

library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("competitor.R")
source("CI.R")



# data preprocessing
cut_off = 0.5
num_node = dim(data_year_product_tensor)[1]
num_layer = dim(data_year_product_tensor)[3]
num_time = dim(data_year_product_tensor)[4]

data <- array(NA,dim=c(num_time, num_node, num_node, num_layer))

for (t in 1:num_time){
  data_year = data_year_product_tensor[, , , t] # 75 75  4
  t_cut_off = quantile(as.numeric(data_year), cut_off)
  data_year[data_year <= t_cut_off] = 0
  data_year[data_year > t_cut_off] = 1
  data[t, , ,] = data_year
}

dim(data) # 35 75 75  4

A.tensor <- data[1:34,,,] # even number time span
dim(A.tensor)




# proposed method
hat.rank <- c(10, 10, num_layer)
threshold_ratio <- c(0.2, 0.3, 0.4, 0.5)
threshold_list <- rev(threshold_ratio * num_node*sqrt(num_layer)*(log(num_time/2))^(3/2))
intervals <- construct_intervals(num_time/2, sqrt(1/2), 4)




A.tensor.even <- A.tensor[seq(2, num_time, by = 2), , , ]
B.tensor.odd  <- A.tensor[seq(1, num_time-1, by = 2), , , ] 

gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = TRUE, intervals, obj.B = B.tensor.odd)
results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_time/2, CUSUM_res = gains, verbose = FALSE,
                               threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)

output <- list()

for (i in 1:(length(results_g)-1)) {
  detected_CP_g <- sort(results_g[[i+1]]$results[, 1])
  detected_CP_gl1 <- refinement1(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
  CI <- construct_CI(alpha=0.05, detected_CP_gl1, A.tensor.even, B.tensor.odd, hat.rank)
  
  output[[i]] <- list()
  output[[i]]$threshold <- results_g[[i+1]]$threshold
  output[[i]]$thres_ratio <- rev(threshold_ratio)[i]
  output[[i]]$detected_CP <- data_names[2* detected_CP_gl1]
  output[[i]]$CI <- CI*2
  
}




# competitor methods
result_kerSeg_net <- Evaluation_kerSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "kerSeg_net")
result_kerSeg_fro <- Evaluation_kerSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "kerSeg_fro")
result_gSeg_net <- Evaluation_gSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "gSeg_net")
result_gSeg_fro <- Evaluation_gSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "gSeg_fro")


# save the result
kerSeg_net <- data_names[result_kerSeg_net]
kerSeg_fro <- data_names[result_kerSeg_fro]
gSeg_net <- data_names[result_gSeg_net]
gSeg_fro <- data_names[result_gSeg_fro]


# printing
output
kerSeg_net
kerSeg_fro
gSeg_net
gSeg_fro


#save(output, file = paste0("real_data/results/trade_proposed.RData"))
#save(kerSeg_net, file = paste0("real_data/results/trade_kerSeg_net.RData"))
#save(kerSeg_fro, file = paste0("real_data/results/trade_kerSeg_fro.RData"))
#save(gSeg_net, file = paste0("real_data/results/trade_gSeg_net.RData"))
#save(gSeg_fro, file = paste0("real_data/results/trade_gSeg_fro.RData"))


