

################################
# Evaluation Function (kerSeg) #
################################


library(kerSeg)

Evaluation_kerSeg <- function(y_array, p_threshold, is_experiment=FALSE, true_CP, pre_process){
  
  # is_experiment = TRUE  # for real data experiment
  # is_experiment = FALSE # for simulated study
  # pre_process: kerSeg_net, kerSeg_fro
  
  binary_search <- function(temp, holder, lower_index, upper_index, p_threshold){
    
    if(dim(temp)[1] > 5){
      
      K = gaussiankernel(temp) # Gaussian kernel matrix
      a = kerseg1(dim(temp)[1], K, pval.perm=TRUE, B=1000)
      cur_index <- a$tauhat 
      cat('lower: ', lower_index, '\n')
      cat('upper: ', upper_index, '\n')
      cat('dim: ', dim(temp)[1], '\n')
      cat('est_CP: ', cur_index, '\n')
      cat('actual index: ', lower_index + cur_index - 1, '\n')
      cat('p-value: ',a$appr$fGKCP1_bon, '\n')
      
      if(a$appr$fGKCP1_bon < p_threshold){
        holder <- c(holder, lower_index+cur_index-1)
        cat('holder: ', holder, '\n')
        
        cat('left\n')
        temp1 <- temp[1:cur_index,]
        holder <- binary_search(temp1, holder, lower_index, lower_index+cur_index-1, p_threshold)
        cat('right\n')
        temp2 <- temp[(cur_index+1):dim(temp)[1],]
        holder <- binary_search(temp2, holder, lower_index+cur_index, upper_index, p_threshold)
      }
      
    }else{cat('less than 5\n')}
    
    return(holder)
    
  }
  
  if(!is_experiment){
    
    repetition_iteration <- dim(y_array)[1] 
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_array[rep_iter,,,,] # 150  50  50   4
      num_time <- dim(y_data)[1]
      n <- dim(y_data[1,,,])[1] # 50
      L <- dim(y_data[1,,,])[3] # 4
      
      
      if(pre_process == "kerSeg_net"){
        
        # vectorize network
        temp <- matrix(NA, nrow=num_time, ncol=n*n*L)
        for(iter in 1:num_time){temp[iter,] <- c(y_data[iter,,,])}
        
      }else if(pre_process == "kerSeg_fro"){
        
        # Frobenius norm by layer
        temp <- matrix(NA, nrow=num_time, ncol=L)
        for(iter in 1:num_time){
          for (l in 1:L) {
            temp[iter, l] <- norm(as.matrix(y_data[iter, , , l]), type = "F")  
          }
        }
        
      }
      
      est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
      est_CP <- sort(est_CP)
      
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(0, true_CP, num_time)
      est_CP_corrected <- c(0, est_CP, num_time)
      
      gt_list <- est_list <- list()
      for (i in 2:length(gt_CP_corrected)) {
        gt_list[[i-1]] <- (gt_CP_corrected[i-1] + 1):gt_CP_corrected[i]
      }
      for (i in 2:length(est_CP_corrected)) {
        est_list[[i-1]] <- (est_CP_corrected[i-1] + 1):est_CP_corrected[i]
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/num_time
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
    }
    
    return(result)
    
  }
  
}



##############################
# Evaluation Function (gSeg) #
##############################


library(gSeg)
library(ade4)

Evaluation_gSeg <- function(y_array, p_threshold, is_experiment=FALSE, true_CP, pre_process){
  
  # is_experiment = TRUE  # for real data experiment
  # is_experiment = FALSE # for simulated study
  # pre_process: gSeg_net, gSeg_fro
  
  binary_search <- function(temp, holder, lower_index, upper_index, p_threshold){
    
    if(dim(temp)[1] > 5){
      
      temp_dist <- dist(temp)
      temp_mst <- mstree(temp_dist)
      temp_result <- gseg1(dim(temp)[1], temp_mst, statistics="o") # "m" "w" "g"
      cur_index <- temp_result$scanZ$ori$tauhat # temp_result$scanZ$max$tauhat # temp_result$scanZ$weighted$tauhat # temp_result$scanZ$generalized$tauhat
      p_value <- temp_result$pval.appr$ori # temp_result$pval.appr$max # temp_result$pval.appr$weighted # temp_result$pval.appr$generalized
      
      cat('lower: ', lower_index, '\n')
      cat('upper: ', upper_index, '\n')
      cat('dim: ', dim(temp)[1], '\n')
      cat('est_CP: ', cur_index, '\n')
      cat('actual index: ', lower_index + cur_index - 1, '\n')
      cat('p-value: ',p_value, '\n')
      
      if(p_value < p_threshold){
        holder <- c(holder, lower_index+cur_index-1)
        cat('holder: ', holder, '\n')
        
        cat('left\n')
        temp1 <- temp[1:cur_index,]
        holder <- binary_search(temp1, holder, lower_index, lower_index+cur_index-1, p_threshold)
        cat('right\n')
        temp2 <- temp[(cur_index+1):dim(temp)[1],]
        holder <- binary_search(temp2, holder, lower_index+cur_index, upper_index, p_threshold)
      }
      
    }else{cat('less than 5\n')}
    
    return(holder)
    
  }
  
  if(!is_experiment){
    
    repetition_iteration <- dim(y_array)[1] 
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_array[rep_iter,,,,] # 150  50  50   4
      num_time <- dim(y_data)[1]
      n <- dim(y_data[1,,,])[1] # 50
      L <- dim(y_data[1,,,])[3] # 4
      
      if(pre_process == "gSeg_net"){
        
        # vectorize network
        temp <- matrix(NA, nrow=num_time, ncol=n*n*L)
        for(iter in 1:num_time){temp[iter,] <- c(y_data[iter,,,])}
        
      }else if(pre_process == "gSeg_fro"){
        
        # Frobenius norm by layer
        temp <- matrix(NA, nrow=num_time, ncol=L)
        for(iter in 1:num_time){
          for (l in 1:L) {
            temp[iter, l] <- norm(as.matrix(y_data[iter, , , l]), type = "F")  
          }
        }
        
      }
      
      est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
      est_CP <- sort(est_CP)
      
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(0, true_CP, num_time)
      est_CP_corrected <- c(0, est_CP, num_time)
      
      gt_list <- est_list <- list()
      for (i in 2:length(gt_CP_corrected)) {
        gt_list[[i-1]] <- (gt_CP_corrected[i-1] + 1):gt_CP_corrected[i]
      }
      for (i in 2:length(est_CP_corrected)) {
        est_list[[i-1]] <- (est_CP_corrected[i-1] + 1):est_CP_corrected[i]
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/num_time
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
    }
    
    return(result)
    
  }
  
}





#########
# USAGE #
#########
#load(file.choose()) # choose RData
#true_CP <- c(40,60)
#Evaluation_kerSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_CP, "kerSeg_net")
#Evaluation_gSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_CP, "gSeg_net")


