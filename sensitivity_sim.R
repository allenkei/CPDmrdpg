library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("eval.R")


########
# DEMO #
########

true_CP <- c(40, 60) # Sce 1
# true_CP <- c(50,100) # Sce 2, 5
# true_CP <- c(50,100,150,200,250) # Sce 3
# true_CP <- c() # Sce 4

load("data/seq10n50s1.RData") # Scenario 1 with node 50

num_seq <- dim(A.all_seq)[1] # 10 sequences
num_T <- dim(A.all_seq)[2] # 150 time points
num_node <- dim(A.all_seq)[3] 
num_layer <- dim(A.all_seq)[5] 
hat.rank <- c(15, 15, num_layer) # needed for model selection (Question: should be used as input to some FUNC)

threshold_list <- rev(c(0.2, 0.4, 0.5, 0.6, 0.8) * num_node*sqrt(num_layer)*(log(num_T/2))^(3/2))

seq_iter <- 1 # used to test INSIDE the for-loop
# output_holder <- matrix(NA, nrow = num_seq, ncol = 4) # 4 metrics
intervals <- construct_intervals(num_T/2, sqrt(1/2), 4)

output_holder_g <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_gl1 <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_gl2 <- array(NA, dim = c(num_seq, length(threshold_list), 4))

output_holder_n <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_nl1 <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_nl2 <- array(NA, dim = c(num_seq, length(threshold_list), 4))

# report mean of metric for all simulated sequences
# can suppress print statements with verbose = FALSE (default TRUE)
for(seq_iter in 1:num_seq){
  cat("\nIteration", seq_iter, "begin.\n")
  
  A.tensor <- A.all_seq[seq_iter,,,,] # a particular sequence with dim 150  50  50   4
  
  # splitting data in half
  A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
  B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
  
  gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = FALSE, intervals, obj.B = B.tensor.odd)
  results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                 threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
  
  for (i in 1:length(threshold_list)) {
    detected_CP_g <- sort(results_g[[i+1]]$results[, 1])
    
    detected_CP_gl1 <- refinement1(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
    detected_CP_gl2 <- refinement2(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
    
    output_holder_g[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_g, num_T))
    output_holder_gl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_gl1, num_T))
    output_holder_gl2[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_gl2, num_T))
    
    cat("Threshold: ", threshold_list[i], "\n")
    cat("\tDetected Greedy CP  :", 2*detected_CP_g, ". Metrics: ", output_holder_g[seq_iter, i, ], "\n")
    cat("\tRefinement Greedy   :", 2*detected_CP_gl1, ". Metrics: ", output_holder_gl1[seq_iter, i, ], "\n")
    cat("\tRefinement 2 Greedy :", 2*detected_CP_gl2, ". Metrics: ", output_holder_gl2[seq_iter, i, ], "\n")
  }
  
  results_n <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                 threshold = threshold_list, method = "Narrowest", obj.B = B.tensor.odd)

  for (i in 1:length(threshold_list)) {
    detected_CP_n <- sort(results_n[[i+1]]$results[, 1])

    detected_CP_nl1 <- refinement1(detected_CP_n, A.tensor.even, B.tensor.odd, hat.rank)
    detected_CP_nl2 <- refinement2(detected_CP_n, A.tensor.even, B.tensor.odd, hat.rank)

    output_holder_n[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_n, num_T))
    output_holder_nl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_nl1, num_T))
    output_holder_nl2[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_nl2, num_T))

    cat("Threshold: ", threshold_list[i], "\n")
    cat("\tDetected Narrowest CP  :", 2*detected_CP_n, ". Metrics: ", output_holder_n[seq_iter, i, ], "\n")
    cat("\tRefinement Narrowest   :", 2*detected_CP_nl1, ". Metrics: ", output_holder_nl1[seq_iter, i, ], "\n")
    cat("\tRefinement 2 Narrowest :", 2*detected_CP_nl2, ". Metrics: ", output_holder_nl2[seq_iter, i, ], "\n")
  }
}


sce1_50 <- list(
  greedy = output_holder_g,
  greedyl1 = output_holder_gl1,
  greedyl2 = output_holder_gl2,
  narrowest = output_holder_n,
  narrowestl1 = output_holder_nl1,
  narrowestl2 = output_holder_nl2
)
save(sce1_50, file = "results/sce1_50.RData")




load("data/seq10n100s1.RData") # Scenario 1 with node 100

num_seq <- dim(A.all_seq)[1] # 10 sequences
num_T <- dim(A.all_seq)[2] # 150 time points
num_node <- dim(A.all_seq)[3] 
num_layer <- dim(A.all_seq)[5] 
hat.rank <- c(15, 15, num_layer) # needed for model selection (Question: should be used as input to some FUNC)

threshold_list <- rev(c(0.2, 0.4, 0.5, 0.6, 0.8) * num_node*sqrt(num_layer)*(log(num_T/2))^(3/2))

seq_iter <- 1 # used to test INSIDE the for-loop
# output_holder <- matrix(NA, nrow = num_seq, ncol = 4) # 4 metrics
intervals <- construct_intervals(num_T/2, sqrt(1/2), 4)

output_holder_g <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_gl1 <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_gl2 <- array(NA, dim = c(num_seq, length(threshold_list), 4))

output_holder_n <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_nl1 <- array(NA, dim = c(num_seq, length(threshold_list), 4))
output_holder_nl2 <- array(NA, dim = c(num_seq, length(threshold_list), 4))

# report mean of metric for all simulated sequences
# can suppress print statements with verbose = FALSE (default TRUE)
for(seq_iter in 1:num_seq){
  cat("\nIteration", seq_iter, "begin.\n")
  
  A.tensor <- A.all_seq[seq_iter,,,,] # a particular sequence with dim 150  50  50   4
  
  # splitting data in half
  A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
  B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
  
  gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = FALSE, intervals, obj.B = B.tensor.odd)
  results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                 threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
  
  for (i in 1:length(threshold_list)) {
    detected_CP_g <- sort(results_g[[i+1]]$results[, 1])
    
    detected_CP_gl1 <- refinement1(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
    detected_CP_gl2 <- refinement2(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
    
    output_holder_g[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_g, num_T))
    output_holder_gl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_gl1, num_T))
    output_holder_gl2[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_gl2, num_T))
    
    cat("Threshold: ", threshold_list[i], "\n")
    cat("\tDetected Greedy CP  :", 2*detected_CP_g, ". Metrics: ", output_holder_g[seq_iter, i, ], "\n")
    cat("\tRefinement Greedy   :", 2*detected_CP_gl1, ". Metrics: ", output_holder_gl1[seq_iter, i, ], "\n")
    cat("\tRefinement 2 Greedy :", 2*detected_CP_gl2, ". Metrics: ", output_holder_gl2[seq_iter, i, ], "\n")
  }
  
  results_n <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                 threshold = threshold_list, method = "Narrowest", obj.B = B.tensor.odd)
  
  for (i in 1:length(threshold_list)) {
    detected_CP_n <- sort(results_n[[i+1]]$results[, 1])
    
    detected_CP_nl1 <- refinement1(detected_CP_n, A.tensor.even, B.tensor.odd, hat.rank)
    detected_CP_nl2 <- refinement2(detected_CP_n, A.tensor.even, B.tensor.odd, hat.rank)
    
    output_holder_n[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_n, num_T))
    output_holder_nl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_nl1, num_T))
    output_holder_nl2[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_nl2, num_T))
    
    cat("Threshold: ", threshold_list[i], "\n")
    cat("\tDetected Narrowest CP  :", 2*detected_CP_n, ". Metrics: ", output_holder_n[seq_iter, i, ], "\n")
    cat("\tRefinement Narrowest   :", 2*detected_CP_nl1, ". Metrics: ", output_holder_nl1[seq_iter, i, ], "\n")
    cat("\tRefinement 2 Narrowest :", 2*detected_CP_nl2, ". Metrics: ", output_holder_nl2[seq_iter, i, ], "\n")
  }
  
}

sce1_100 <- list(
  greedy = output_holder_g,
  greedyl1 = output_holder_gl1,
  greedyl2 = output_holder_gl2,
  narrowest = output_holder_n,
  narrowestl1 = output_holder_nl1,
  narrowestl2 = output_holder_nl2
)
save(sce1_100, file = "results/sce1_100.RData")

