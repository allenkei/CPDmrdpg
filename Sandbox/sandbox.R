library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("model_selection.R")
source("eval.R")
#library(devtools) # install.packages("devtools")
#install_github("etam4260/kneedle") # install the package "kneedle" via "devtools"
library(kneedle)

#############
# LOAD DATA #
#############
# MANUALLY UN-COMMENT THE LINE TO LOAD DATA "A.all_seq"

scenario <- "s5" # CHOOSE "s1","s2","s3","s4","s5"

#load("data/seq10n50s1.RData")  # Scenario 1 with node 50
#load("data/seq10n100s1.RData") # Scenario 1 with node 100
#load("data/seq10n50s2.RData")  # Scenario 2 with node 50
#load("data/seq10n100s2.RData") # Scenario 2 with node 100
#load("data/seq10n50s3.RData")  # Scenario 3 with node 50
#load("data/seq10n100s3.RData") # Scenario 3 with node 100
#load("data/seq10n50s4.RData")  # Scenario 4 with node 50
#load("data/seq10n100s4.RData") # Scenario 4 with node 100
load("data/seq10n50s5.RData")  # Scenario 5 with node 50
#load("data/seq10n100s5.RData") # Scenario 5 with node 100


#################
# CONFIGURATION #
#################





if(scenario %in% c("s1", "s2", "s5")){
  
  true_CP <- c(50,100) # "s1", "s2", "s5" have 2 CP
  
}else if(scenario == "s3") {
  
  true_CP <- c(50,100,150,200,250) # "s3" has 5 CP
  
}else if(scenario == "s4") {
  
  true_CP <- c() # "s4" has no CP
  
}





dim(A.all_seq) # 10 150  50  50   4

num_seq <- dim(A.all_seq)[1] # 10 sequences
num_T <- dim(A.all_seq)[2] # 150 or 300 time points
num_node <- dim(A.all_seq)[3]  # 50 or 100 nodes
num_layer <- dim(A.all_seq)[5] # 4 layers
hat.rank <- c(15, 15, 15) # needed for model selection (Question: should be used as input to some FUNC)


#threshold <- num_node*num_layer*sqrt(log(num_T))

# construct intervals (FIXED for all sequences)
intervals <- construct_intervals(num_T/2, sqrt(1/2), 2) # half of full time span






####################
# SIMULATION STUDY #
####################


seq_iter <- 1 # used to test INSIDE the for-loop


output_holder <- matrix(NA, nrow = num_seq, ncol = 4) # 4 metrics


# report mean of metric for all simulated sequences
# can suppress print statements with verbose = FALSE (default TRUE)
for(seq_iter in 1:num_seq){
  
  # if(seq_iter == 2) break # break at seq_iter=1 for checking
  
  A.tensor <- A.all_seq[seq_iter,,,,] # a particular sequence with dim 150  50  50   4
  
  # splitting data in half
  A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
  B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
  
  # CUSUM_step1 is a FUNC from CUSUM.R
  # obtain CP candidate for each interval
  results_all_step1 <- cusum_on_intervals(CUSUM_step1, A.tensor.even, intervals, verbose = TRUE, obj.B = B.tensor.odd)
  
  # CUSUM_step1 is a FUNC from CUSUM.R
  # obtain initial result for threshold
  init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                            threshold = 0, method = "Greedy", verbose = TRUE, obj.B = B.tensor.odd)
  
  # construct list of threshold
  max <- max(init[[1]]$results[, 2])
  threshold_list = c(max + 1, init[[2]]$results[, 2][1:10]) # Maximum number of changepoints 
  
  
  results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                  threshold = threshold_list, method = "Greedy", verbose = TRUE, obj.B = B.tensor.odd)
  
  
  # Model selection
  # Currently, all methods need hat.rank
  out <- model_selection(results_ms, A.tensor, method = "elbow", verbose = TRUE, hat.rank = hat.rank)
  
  
  # Visualization
  vals <- out[[4]]
  ncps <- out[[5]]
  plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
  text(ncps, vals, ncps, pos = 1, cex = 0.8)
  symbols(length(out[[1]]), out[[2]], circles = 1, add = TRUE, inches = 0.08, fg = "red")
  text(length(out[[1]]), out[[2]], "Knee", pos = 3, cex = 0.8, col = "red")
  
  
  metric_list <- eval_CP(true_CP, detected_CP = out[[1]], num_T)
  
  cat("detected CP after model selection:", out[[1]], "\n")
  cat("metrics: ",metric_list[[1]], metric_list[[2]], metric_list[[3]], metric_list[[4]])
  
  output_holder[seq_iter ,] <- c(metric_list[[1]], metric_list[[2]], metric_list[[3]], metric_list[[4]])
  break
}

out[[1]]
results_ms[[12]]

