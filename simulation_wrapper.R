library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("eval.R")
source("gen_data.R")
source("CI.R")


simulate_sensitivity <- function(scenario, true_cp, num_node = 50, num_seq = 10) {
  # Simulate_sensitivity computes 4 metrics
  # For a range of thresholds
  # Returning a 3-way tensor
  A.all_seq <- generate(scenario, num_node, 1, FALSE)
  
  num_T <- dim(A.all_seq)[2] 
  num_node <- dim(A.all_seq)[3] 
  num_layer <- dim(A.all_seq)[5] 
  hat.rank <- c(15, 15, num_layer) # needed for model selection (Question: should be used as input to some FUNC)
  
  threshold_list <- rev(c(0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25) * num_node*sqrt(num_layer)*(log(num_T/2))^(3/2))
  
  intervals <- construct_intervals(num_T/2, sqrt(1/2), 4)
  
  output_holder_g <- array(NA, dim = c(num_seq, length(threshold_list), 4))
  output_holder_gl1 <- array(NA, dim = c(num_seq, length(threshold_list), 4))
  
  # report mean of metric for all simulated sequences
  # can suppress print statements with verbose = FALSE (default TRUE)
  for(seq_iter in 1:num_seq) {
    cat("\nIteration", seq_iter, "begin.\n")
    set.seed(seq_iter)
    # Generate Data 1-by-1
    A.all_seq <- generate(scenario, num_node, 1, FALSE)
    A.tensor <- A.all_seq[1,,,,] # a particular sequence with dim 150  50  50   4
    
    # splitting data in half
    A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
    B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
    
    gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = FALSE, intervals, obj.B = B.tensor.odd)
    results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                   threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
    
    for (i in 1:(length(results_g)-1)) {
      detected_CP_g <- sort(results_g[[i+1]]$results[, 1])
      detected_CP_gl1 <- refinement1(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
      
      output_holder_g[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_g, num_T))
      output_holder_gl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_gl1, num_T))
      
      cat("Threshold: ", threshold_list[i], "\n")
      cat("\tDetected Greedy CP  :", 2*detected_CP_g, "Metrics: ", output_holder_g[seq_iter, i, ], "\n")
      cat("\tRefinement Greedy   :", 2*detected_CP_gl1, "Metrics: ", output_holder_gl1[seq_iter, i, ], "\n")
    }
  }
  
  results <- list(
    greedy = output_holder_g,
    greedyl1 = output_holder_gl1
  )
  #save(results, file = paste0("results/", scenario, "_", num_node, ".RData"))
  
  return(results)
}


simulate_coverage <- function(scenario, true_cp, num_node = 50, num_seq = 10, 
                              threshold_C = 0.1, alpha = 0.05) {
  # Simulate_coverage returns 4 metrics & interval coverage
  # For a single threshold
  # Returning a matrix
  A.all_seq <- generate(scenario, num_node, 1, FALSE)
  
  num_T <- dim(A.all_seq)[2] 
  num_node <- dim(A.all_seq)[3] 
  num_layer <- dim(A.all_seq)[5] 
  hat.rank <- c(15, 15, num_layer) # needed for model selection (Question: should be used as input to some FUNC)
  
  intervals <- construct_intervals(num_T/2, sqrt(1/2), 4)
  
  output_holder <- array(NA, dim = c(num_seq, 4))
  coverage_holder <- array(NA, dim = c(num_seq, length(true_cp)))
  lengths_holder <- array(NA, dim = c(num_seq, length(true_cp)))
  
  # report mean of metric for all simulated sequences
  # can suppress print statements with verbose = FALSE (default TRUE)
  for(seq_iter in 1:num_seq) {
    cat("\nIteration", seq_iter, "begin.\n")
    set.seed(seq_iter)
    # Generate Data 1-by-1
    A.all_seq <- generate(scenario, num_node, 1, FALSE)
    A.tensor <- A.all_seq[1,,,,] # a particular sequence with dim 150  50  50   4
    
    # splitting data in half
    A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
    B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
    
    gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = FALSE, intervals, obj.B = B.tensor.odd)
    results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                   threshold = threshold_C * num_node*sqrt(num_layer)*(log(num_T/2))^(3/2), 
                                   method = "Greedy", obj.B = B.tensor.odd)
    detected_CP <- refinement1(sort(results_g[[2]]$results[, 1]), A.tensor.even, B.tensor.odd, hat.rank)
    output_holder[seq_iter, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP, num_T))
    cat("\tDetected Greedy CP  :", 2*detected_CP, "Metrics: ", output_holder[seq_iter, ], "\n")
    
    
    CI <- construct_CI(alpha, detected_CP, A.tensor.even, B.tensor.odd, hat.rank)
    coverage_out <- coverage(true_cp, CI[, 2]*2, CI[, 3]*2, CI[, 4]*2)
    coverage_holder[seq_iter, ] <- coverage_out$covered
    lengths_holder[seq_iter, ] <- coverage_out$lengths
    cat("\tRefined (Coverage) CP :", 2*CI[, 2], "Coverage: ", coverage_holder[seq_iter, ], 
        "Lengths: ", lengths_holder[seq_iter, ],"\n")
  }
  
  results <- output_holder
  results <- list(
    output_holder,
    coverage_holder,
    lengths_holder
  )
  save(results, file = paste0("results/coverage_", scenario, "_", num_node, ".RData"))
  
  return(results)
}
 
###########
# Run one #
###########

scenario <- "f1" # f1, f2, f3, f4, f5, f6
if (scenario == "f1") {
  true_CP <- c(70, 140)
} else if (scenario == "f2") {
  true_CP <- c(20, 60, 80, 160, 180)
} else if (scenario == "f3") {
  true_CP <- c(50, 100, 150)
} else if (scenario == "f4") {
  true_CP <- c(20, 60, 80, 160, 180)
} else if (scenario == "f5") {
  true_CP <- c(50, 100, 150)
} else if (scenario == "f6") {
  true_CP <- c()
} else {
  stop("Invalid scenario!")
}

num_node <- 100
num_seq <- 100

#results <- simulate_sensitivity(scenario, true_CP, num_node, num_seq)
results <- simulate_coverage(scenario, true_CP, num_node, num_seq)

###########
# Run all #
###########


{
  timing_summary <- data.frame(
    scenario = character(),
    elapsed_minutes = numeric(),
    stringsAsFactors = FALSE
  )
  
  start_time <- Sys.time()
  cat("Start time:", format(start_time), "\n")
  for (scenario in c("f1","f2", "f3", "f4", "f5")) {
    
    num_node <- 50
    num_seq <- 100
    
    if (scenario == "f1") {
      true_CP <- c(70, 140)
    } else if (scenario == "f2") {
      true_CP <- c(20, 60, 80, 160, 180)
    } else if (scenario == "f3") {
      true_CP <- c(50, 100, 150)
    } else if (scenario == "f4") {
      true_CP <- c(20, 60, 80, 160, 180)
    } else if (scenario == "f5") {
      true_CP <- c(50, 100, 150)
    } else if (scenario == "f6") {
      true_CP <- c()
    } else {
      stop("Invalid scenario!")
    }
    
    loc_start <- Sys.time()
    cat("\n==== Running scenario", scenario, "====\n")
    cat("Start time:", format(loc_start), "\n")
    
    # results <- simulate_sensitivity(scenario, true_CP, num_node, num_seq)
    results <- simulate_coverage(scenario, true_CP, num_node, num_seq)
    
    end_time <- Sys.time()
    elapsed <- difftime(end_time, loc_start, units = "mins")
    cat("\nEnd time:", format(end_time), "\n")
    cat("Elapsed time:", round(elapsed, 2), "minutes\n")
    
    # Log timing
    timing_summary <- rbind(
      timing_summary,
      data.frame(scenario = scenario, elapsed_minutes = round(as.numeric(elapsed), 2))
    )
  }
  
  timing_summary
}

############
# Analysis #
############

# For multiple-threshold results
load("results/f6_50.RData")

# One metric, refinement and original 
num_thresholds <- dim(results[[1]])[2]
summary_matrix <- matrix(NA, nrow = length(results), ncol = num_thresholds)

for (i in 1:length(results)) {
  # ADJUST FOR DESIRED STAT
  stat_matrix <- results[[i]][, , 4] 
  summary_matrix[i, ] <- colMeans(stat_matrix, na.rm = TRUE)
  # summary_matrix[i, ] <- apply(stat_matrix, 2, median, na.rm = TRUE)
  # summary_matrix[i, ] <- apply(stat_matrix, 2, function(x) sum(x != 0, na.rm = TRUE))
  # summary_matrix[i, ] <- apply(stat_matrix, 2, function(x) sum(is.finite(x)))
}

rownames(summary_matrix) <- c("Greedy", "Greedy1")
colnames(summary_matrix) <- paste0(rev(c(0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)))
summary_matrix

# All metrics, refinement only 
load("results/f6_100.RData")
num_thresholds <- dim(results[[1]])[2]
summary_matrix <- matrix(NA, nrow = 4, ncol = num_thresholds)

for (i in 1:4) {
  # ADJUST FOR DESIRED STAT
  stat_matrix <- results[[2]][, , i] 
  summary_matrix[i, ] <- colMeans(stat_matrix, na.rm = TRUE)
  # summary_matrix[i, ] <- apply(stat_matrix, 2, median, na.rm = TRUE)
  # summary_matrix[i, ] <- apply(stat_matrix, 2, function(x) sum(x != 0, na.rm = TRUE))
  # summary_matrix[i, ] <- apply(stat_matrix, 2, function(x) sum(is.finite(x)))
}

colnames(summary_matrix) <- paste0(rev(c(0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)))
summary_matrix


############
# Analysis #
############
# For single-threshold results (with coverage)
load("results/coverage_f5_50.RData")

colMeans(results[[1]]) # Original 4 metrics
results[[2]] # Coverage indicators
results[[3]] # Lengths, (-1 if not estimated)

# Number of non-evaluated
apply(results[[2]], 2, function(col) {sum(col == -1, na.rm = TRUE)})

# Coverage per change point, filtering -1
apply(results[[2]], 2, function(col) {mean(col[col != -1], na.rm = TRUE)})
# Average coverage, filtering -1
mean(results[[2]][results[[2]] != -1], na.rm = TRUE)

# Lengths per change point, filtering -1
apply(results[[3]], 2, function(col) {mean(col[col != -1], na.rm = TRUE)})
# Average lenghts, filtering -1
mean(results[[3]][results[[3]] != -1], na.rm = TRUE)

# Median lengths per change point, filtering -1
apply(results[[3]], 2, function(col) {median(col[col != -1], na.rm = TRUE)})



# f1 100% overall, 
  # 100% each, 
  # 0 not evaluated, 
  # length 0.01357

# f2 79.4% overall, 
  # 41%, 100%, 100%, 92%, 64% 
  # 0 not evaluated
  # length 1.517

# f3 56.2% overall, 
  # 17%, 64%, 81%
  # 19 not evaluated (17, 2, 0)
  # length 7.57

# f4 100% overall, 
  # 100% each, 
  # 0 not evaluated, 
  # length 1.962

# f5 100% overall, 
  # 100% each, 
  # 80 not evaluated (35, 45, 0)
  # length 1.961


