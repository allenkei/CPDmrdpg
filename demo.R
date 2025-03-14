library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("model_selection.R")
source("eval.R")
#library(devtools) # install.packages("devtools")
#install_github("etam4260/kneedle") # install the package "kneedle" via "devtools"
library(kneedle)


########
# DEMO #
########

load("data/seq10n100s3.RData") # Scenario 1 with node 50
# load("data/seq10n100s1.RData") # Scenario 1 with node 100
# load("data/seq10n50s2.RData") # Scenario 2 with node 50
# load("data/seq10n100s2.RData") # Scenario 2 with node 100

dim(A.all_seq) # 10 150  50  50   4

num_seq <- dim(A.all_seq)[1] # 10 sequences
num_T <- dim(A.all_seq)[2] # 150 time points
num_node <- dim(A.all_seq)[3] 
num_layer <- dim(A.all_seq)[5] 
hat.rank <- c(15, 15, 15) # needed for model selection (Question: should be used as input to some FUNC)
# true_CP <- c(50,100) # Sce 1, 2
true_CP <- c(50,100,150,200,250,300,350) # Sce 3
# true_CP <- c() # Sce 4

threshold <- num_node*num_layer*sqrt(log(num_T))

# construct intervals (FIXED for all sequences)
intervals <- construct_intervals(num_T/2, sqrt(1/2), 2) # half of full time span

seq_iter <- 1 # used to test INSIDE the for-loop

output_holder <- matrix(NA, nrow = num_seq, ncol = 4) # 4 metrics


# report mean of metric for all simulated sequences
# can suppress print statements with verbose = FALSE (default TRUE)
for(seq_iter in 1:num_seq){
  
  # if(seq_iter == 6) break 
  
  A.tensor <- A.all_seq[seq_iter,,,,] # a particular sequence with dim 150  50  50   4
  
  # splitting data in half
  A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
  B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
  
  # threshold = threshold*c(1, 0.5, 0.1) # If I want to test more
  results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, threshold = threshold,
                               method = "Narrowest", verbose = FALSE, obj.B = B.tensor.odd)
  
  detected_CP <- 2*sort(results[[2]]$results[, 1])
  
  metric_list <- eval_CP(true_CP, detected_CP, num_T)
  
  cat("\nIteration", seq_iter, "detected CP:", detected_CP, "\n")
  cat("\tthreshold: ",threshold, "\n")
  cat("\tmetrics: ",as.numeric(metric_list), "\n")
  print(results[[2]]$results)
  
  output_holder[seq_iter ,] <- c(metric_list[[1]], metric_list[[2]], metric_list[[3]], metric_list[[4]])
}

output_holder


  ##################################
# Example Code (Not CUSUM_step1) # 
##################################

A.tensor <- A.all_seq[1,,,,]

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
results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, num_T, CUSUM_res = results_all_Layer, 
                             threshold = c(10, 9, 8), method = "Greedy", rank = hat.rank)

results[[2]]$results

##############################
# Example Code (CUSUM_step1) # 
##############################

A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ]
if (dim(A.tensor.even)[1] != dim(B.tensor.odd)[1]) {
  stop("Make sure even and odd have same length")
}

hat.rank <- c(15, 15, 15)
s <- 0
e <- 75
frobenius_holder <- numeric(74)

for(t in 2:(num_T/2-2)){
  frobenius_holder[t] <- CUSUM_step1(A.tensor.even, s, e, t, obj.B = B.tensor.odd)
}

plot(1:(num_T/2-1), frobenius_holder, type='l')

intervals <- construct_intervals(num_T/2, sqrt(1/2), 2)
results_one <- cusum_on_intervals(CUSUM_step1, A.tensor.even, c(30, 38), obj.B = B.tensor.odd)
results_all_step1 <- cusum_on_intervals(CUSUM_step1, A.tensor.even, intervals, obj.B = B.tensor.odd)


# Calculate intervals & CUSUM within
# results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, alpha = 1/2, m = 4,
#                              threshold = 200, method = "Narrowest", obj.B = B.tensor.odd)

# Pass in CUSUM results (as a matrix) to save computation
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                             threshold = c(1800, 900, 500, 200, 1), method = "Narrowest", obj.B = B.tensor.odd)
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                             threshold = c(1000, 500, 250, 70, 1), method = "Greedy", obj.B = B.tensor.odd)

mean(results[[6]]$results)
results[[5]]

###################
# Model Selection # 
###################

source("model_selection.R")
# Kyle's guess for best method
# Thresholding with 10 changepoints
# Elbow method via kneedle: https://ieeexplore.ieee.org/document/5961514
init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                          threshold = 0, method = "Greedy", obj.B = B.tensor.odd)
max <-max(init[[1]]$results[, 2])
init[[2]]$results[, 2]

threshold_list = c(max + 1, init[[2]]$results[, 2][1:10])
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- out[[4]]
ncps <- out[[5]]
plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
text(ncps, vals, ncps, pos = 1, cex = 0.8)
symbols(length(out[[1]]), out[[2]], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(length(out[[1]]), out[[2]], "Knee", pos = 3, cex = 0.8, col = "red")



################################
# More Model Selection Options # 
################################
source("utility.R")
# using utilities/model selection
# Both AIC and BIC prefer fewest (0) changepoints
# New l1 penalty - which requires tuning parameter
init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                          threshold = 0, method = "Greedy", obj.B = B.tensor.odd)
max <-max(init[[1]]$results[, 2])
threshold_list <- exp(seq(1, log(max) + 1, length.out=50))
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)

model_selection(results_ms, A.tensor, method = "l1", hat.rank = hat.rank, beta = 1)[1:3]  # many 
model_selection(results_ms, A.tensor, method = "l1", hat.rank = hat.rank, beta = 2)[1:3]  # (default) many 
model_selection(results_ms, A.tensor, method = "l1", hat.rank = hat.rank, beta = 5)[1:3]  # {50, 100}
model_selection(results_ms, A.tensor, method = "l1", hat.rank = hat.rank, beta = 10)[1:3] # {50, 100}
model_selection(results_ms, A.tensor, method = "l1", hat.rank = hat.rank, beta = 20)[1:3] # {0}

out <- model_selection(results_ms, A.tensor, method = "MDL", hat.rank = hat.rank)
plot(unique(out[[4]]))
out[[1]]

# Other elbow options
# Thresholding using Max Gain
threshold_list <- seq(1, max + 1, length.out=50)
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- unique(out[[4]])
ncps <- unique(out[[5]])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Equally Spaced")
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
symbols(which(length(out[[1]]) == ncps)[1], out[[2]], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(which(length(out[[1]]) == ncps)[1], out[[2]], "Knee", pos = 3, cex = 0.8, col = "red")


# Thresholding using Max Gain, log spacing
threshold_list <- exp(seq(1, log(max) + 1, length.out=50))
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- unique(out[[4]])
ncps <- unique(out[[5]])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Log-Spaced")
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
symbols(which(length(out[[1]]) == ncps)[1], out[[2]], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(which(length(out[[1]]) == ncps)[1], out[[2]], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Changepoints from Greedy, top 20 Gains
threshold_list = c(max + 1, init[[2]]$results[, 2][1:20])
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- out[[4]]
ncps <- out[[5]]
plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
text(ncps, vals, ncps, pos = 1, cex = 0.8)
symbols(length(out[[1]]), out[[2]], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(length(out[[1]]), out[[2]], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Changepoints from Greedy, all
threshold_list = c(max + 1, init[[2]]$results[, 2])
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- out[[4]]
ncps <- out[[5]]
plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
text(ncps, vals, ncps, pos = 1, cex = 0.8)
symbols(length(out[[1]]), out[[2]], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(length(out[[1]]), out[[2]], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Max Gain, < 15 Changepoints
threshold_list <- seq(1, max + 1, length.out=50)
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- unique(out[[4]][out[[5]] < 15])
ncps <- unique(out[[5]][out[[5]] < 15])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Equally Spaced")
mtext("Less than 15 Changepoints", side = 3, line = 0.5, cex = 0.9)
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
# Need to recalculate knee
knee <- kneedle(1:length(ncps), vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Max Gain, log spacing, < 15 Changepoints
threshold_list <- exp(seq(1, log(max) + 1, length.out=50))
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)
vals <- unique(out[[4]][out[[5]] < 15])
ncps <- unique(out[[5]][out[[5]] < 15])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Log-Spaced")
mtext("Less than 15 Changepoints", side = 3, line = 0.5, cex = 0.9)
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
# Need to recalculate knee
knee <- kneedle(1:length(ncps), vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")



