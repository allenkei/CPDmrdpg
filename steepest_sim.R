library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("eval.R")
source("model_selection.R")


########
# DEMO #
########

run_sensitivity <- function(both = TRUE) {
  confirm <- function(label) {
    repeat {
      ans <- readline(paste0("Did you check ", label, "? (y/n): "))
      if (tolower(ans) == "y") break
      if (tolower(ans) == "n") stop(paste("Please fix ", label, " before proceeding."))
      cat("Please answer 'y' or 'n'.\n")
    }
  }
  
  confirm("true_cp")
  confirm("load files for 50 and 100")
  confirm("save files for 50 and 100")
  
  cat("All settings confirmed. Proceeding...\n")
  
  # true_CP <- c(40, 60) # Sce 1, 5
  true_CP <- c(50,100) # Sce 2, 6, 6b
  # true_CP <- c(50,100,150,200,250) # Sce 3
  # true_CP <- c(20, 60, 80, 160, 180) #Sce 3b
  # true_CP <- c() # Sce 4
  # true_CP <- c(20,50,80) # Sce 7, 7b
  
  load("data/seq10n50s2.RData") # Scenario 1 with node 50
  
  num_seq <- dim(A.all_seq)[1] # 10 sequences
  num_T <- dim(A.all_seq)[2] # 150 time points
  num_node <- dim(A.all_seq)[3] 
  num_layer <- dim(A.all_seq)[5] 
  hat.rank <- c(15, 15, num_layer) # needed for model selection (Question: should be used as input to some FUNC)
  
  # c(0.1, 0.2, 0.25, 0.3, 0.4)
  level_list <- rev(c(0.05, 0.1, 0.2, 0.5) * num_node*sqrt(num_layer)*(log(num_T/2))^(3/2))
  
  seq_iter <- 1 # used to test INSIDE the for-loop
  intervals <- construct_intervals(num_T/2, sqrt(1/2), 4)
  
  output_holder_s <- array(NA, dim = c(num_seq, length(level_list), 4))
  output_holder_sl1 <- array(NA, dim = c(num_seq, length(level_list), 4))
  output_holder_r <- array(NA, dim = c(num_seq, length(level_list), 4))
  output_holder_rl1 <- array(NA, dim = c(num_seq, length(level_list), 4))
  
  # report mean of metric for all simulated sequences
  # can suppress print statements with verbose = FALSE (default TRUE)
  for(seq_iter in 1:num_seq) {
    cat("\nIteration", seq_iter, "begin.\n")
    
    A.tensor <- A.all_seq[seq_iter,,,,] # a particular sequence with dim 150  50  50   4
    
    # splitting data in half
    A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
    B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
    
    gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = FALSE, intervals, obj.B = B.tensor.odd)
    init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                                   threshold = 1, method = "Greedy", obj.B = B.tensor.odd)
    
    for (i in 1:length(level_list)) {
      detected_CP_s <- sort(steepest_drop(init[[2]]$results, level_list[i], verbose = FALSE)$candidates)
      detected_CP_sl1 <- refinement1(detected_CP_s, A.tensor.even, B.tensor.odd, hat.rank)
      
      detected_CP_r <- sort(steepest_drop_relative(init[[2]]$results, level_list[i], verbose = FALSE)$candidates)
      detected_CP_rl1 <- refinement1(detected_CP_r, A.tensor.even, B.tensor.odd, hat.rank)
      
      output_holder_s[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_s, num_T))
      output_holder_sl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_sl1, num_T))
      output_holder_r[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_r, num_T))
      output_holder_rl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_rl1, num_T))
      
      cat("Threshold: ", level_list[i], "\n")
      cat("\tDetected SDLL CP :", 2*detected_CP_s, ". Metrics: ", output_holder_s[seq_iter, i, ], "\n")
      cat("\tRefinement SDLL  :", 2*detected_CP_sl1, ". Metrics: ", output_holder_sl1[seq_iter, i, ], "\n")
      cat("\tDetected Rel CP  :", 2*detected_CP_r, ". Metrics: ", output_holder_r[seq_iter, i, ], "\n")
      cat("\tRefinement Rel   :", 2*detected_CP_rl1, ". Metrics: ", output_holder_rl1[seq_iter, i, ], "\n")
    }
  }
  
  sce_50 <- list(
    SDLL = output_holder_s,
    SDLL1 = output_holder_sl1, 
    relSDLL = output_holder_r, 
    relSDLL1 = output_holder_rl1
  )
  save(sce_50, file = "results/S_sce2_50.RData")
  
  
  
  if(both == FALSE) {return ()}
  load("data/seq10n100s2.RData") 
  
  num_seq <- dim(A.all_seq)[1] # 10 sequences
  num_T <- dim(A.all_seq)[2] # 150 time points
  num_node <- dim(A.all_seq)[3] 
  num_layer <- dim(A.all_seq)[5] 
  hat.rank <- c(15, 15, num_layer) # needed for model selection (Question: should be used as input to some FUNC)
  
  # c(0.1, 0.2, 0.25, 0.3, 0.4)
  level_list <- rev(c(0.05, 0.1, 0.2, 0.5) * num_node*sqrt(num_layer)*(log(num_T/2))^(3/2))
  
  seq_iter <- 1 # used to test INSIDE the for-loop
  intervals <- construct_intervals(num_T/2, sqrt(1/2), 4)
  
  output_holder_s <- array(NA, dim = c(num_seq, length(level_list), 4))
  output_holder_sl1 <- array(NA, dim = c(num_seq, length(level_list), 4))
  output_holder_r <- array(NA, dim = c(num_seq, length(level_list), 4))
  output_holder_rl1 <- array(NA, dim = c(num_seq, length(level_list), 4))
  
  # report mean of metric for all simulated sequences
  # can suppress print statements with verbose = FALSE (default TRUE)
  for(seq_iter in 1:num_seq) {
    cat("\nIteration", seq_iter, "begin.\n")
    
    A.tensor <- A.all_seq[seq_iter,,,,] # a particular sequence with dim 150  50  50   4
    
    # splitting data in half
    A.tensor.even <- A.tensor[seq(2, num_T, by = 2), , , ]
    B.tensor.odd  <- A.tensor[seq(1, num_T-1, by = 2), , , ] # named as B.tensor
    
    gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = FALSE, intervals, obj.B = B.tensor.odd)
    init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_T/2, CUSUM_res = gains, verbose = FALSE,
                              threshold = 1, method = "Greedy", obj.B = B.tensor.odd)
    
    for (i in 1:length(level_list)) {
      detected_CP_s <- sort(steepest_drop(init[[2]]$results, level_list[i], verbose = FALSE)$candidates)
      detected_CP_sl1 <- refinement1(detected_CP_s, A.tensor.even, B.tensor.odd, hat.rank)
      
      detected_CP_r <- sort(steepest_drop_relative(init[[2]]$results, level_list[i], verbose = FALSE)$candidates)
      detected_CP_rl1 <- refinement1(detected_CP_r, A.tensor.even, B.tensor.odd, hat.rank)
      
      output_holder_s[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_s, num_T))
      output_holder_sl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_sl1, num_T))
      output_holder_r[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_r, num_T))
      output_holder_rl1[seq_iter, i, ] <- as.numeric(eval_CP(true_CP, 2*detected_CP_rl1, num_T))
      
      cat("Threshold: ", level_list[i], "\n")
      cat("\tDetected SDLL CP :", 2*detected_CP_s, ". Metrics: ", output_holder_s[seq_iter, i, ], "\n")
      cat("\tRefinement SDLL  :", 2*detected_CP_sl1, ". Metrics: ", output_holder_sl1[seq_iter, i, ], "\n")
      cat("\tDetected Rel CP  :", 2*detected_CP_r, ". Metrics: ", output_holder_r[seq_iter, i, ], "\n")
      cat("\tRefinement Rel   :", 2*detected_CP_rl1, ". Metrics: ", output_holder_rl1[seq_iter, i, ], "\n")
    }
  }
  
  sce_100 <- list(
    SDLL = output_holder_s,
    SDLL1 = output_holder_sl1, 
    relSDLL = output_holder_r, 
    relSDLL1 = output_holder_rl1
  )
  save(sce_100, file = "results/S_sce2_100.RData")

}

run_sensitivity()

