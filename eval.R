eval_CP <- function(true_CP, detected_CP, num_T){

  detected_CP <- detected_CP + 1 # Question: is 50 CP or 51 CP (which definition of CP)?
  num_CP <- length(detected_CP)
  
  gt_CP_corrected <- c(1, true_CP, num_T) # c(1,50,100,150)
  detected_CP_corrected <- c(1, detected_CP, num_T)
  
  # intervals
  gt_list <- detected_list <- list();
  for(i in 2:length(gt_CP_corrected)){
    gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
  }
  for(i in 2:length(detected_CP_corrected)){
    detected_list[[i-1]] <- detected_CP_corrected[i-1]:(detected_CP_corrected[i]-1)
  }
  
  # one-sided Hausdorff distance
  if(num_CP == 0){
    dist_detected_gt <- Inf
    dist_gt_detected <- -Inf
    covering_metric <- 0
  }else{
    
    holder <- c()
    for(i in true_CP){
      dist_diff <- c()
      for(j in detected_CP){dist_diff <- c(dist_diff, abs(j-i))}
      holder <- c(holder, min(dist_diff))
    }
    dist_detected_gt <- max(holder)
    
    holder <- c()
    for(i in detected_CP){
      dist_diff <- c()
      for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
      holder <- c(holder, min(dist_diff))
    }
    dist_gt_detected <- max(holder)
    
    # Coverage
    covering_metric <- 0
    for(i in 1:length(gt_list)){
      A <- gt_list[[i]]
      jaccard <- c()
      for(j in 1:length(detected_list)){
        A_prime <- detected_list[[j]]
        jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
      }
      covering_metric <- covering_metric + length(A)*max(jaccard)
    }
    covering_metric <- covering_metric/num_T
    
  }
  
  # difference in number of CP
  abs_error <- abs(num_CP - length(true_CP))
  
  return(list(abs_error = abs_error,
              dist_detected_gt = dist_detected_gt,
              dist_est_detected = dist_gt_detected,
              covering_metric = covering_metric))
  
}


#eval_CP(true_CP = c(50, 100), detected_CP = c(50, 100), num_T = 150)
