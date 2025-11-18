
#download.file(url = "https://github.com/allenkei/CPDmrdpg/archive/master.zip", destfile = "master.zip")
#unzip(zipfile = "master.zip")
#setwd("~/RDPGM/CPDmrdpg-main")


library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("eval.R")
source("gen_data.R")
source("competitor.R")


simulate_scenario_competitor <- function(scenario, true_cp, num_node, num_seq, competitor) {
  # competitor: kerSeg_net, kerSeg_fro, gSeg_net, gSeg_fro
  true_cp_store <- true_cp
  if (is.matrix(true_cp_store)) { true_cp <- as.numeric(true_cp_store[1, ])} 
  temp <- generate(scenario, true_cp, num_node, 1, FALSE)
  
  num_time <- dim(temp)[2]; num_layer <- dim(temp)[5]; rm(temp)
  
  results <- matrix(NA, nrow = num_seq, ncol = 4) # 4 metrics
  
  
  for(seq_iter in 1:num_seq) {
    # Generate Data 1-by-1
    set.seed(seq_iter)
    cat("\nIteration", seq_iter, "begin.\n")
    

    if (is.matrix(true_cp_store)) { true_cp <- as.numeric(true_cp_store[seq_iter, ])} 
    A.all_seq <- generate(scenario, true_cp, num_node, 1, FALSE) # (1,200,50,50,4) # NOTE the first dim = 1
    
    if(competitor == "kerSeg_net"){
      results[seq_iter,] <- Evaluation_kerSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_cp, competitor)
      
    }else if(competitor == "kerSeg_fro"){
      results[seq_iter,] <- Evaluation_kerSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_cp, competitor)
      
    }else if(competitor == "gSeg_net"){
      results[seq_iter,] <- Evaluation_gSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_cp, competitor)
      
    }else if(competitor == "gSeg_fro"){
      results[seq_iter,] <- Evaluation_gSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_cp, competitor)
      
    }
    
    
  }
  
  save(results, file = paste0("results/comp/",scenario,"_n",num_node,"_",competitor,".RData"))
  
  return(results)
}





###########
# Run all #
###########


scenario <- "f1"

for (scenario in c("f1")) { ##### c("f1", "f2", "f3", "f4", "f5")
  
  num_node <- 100
  num_seq <- 100
  
  if (scenario == "f1") {
    true_CP <- true_CP <- c(70, 140)
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
  }
  
  
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "kerSeg_net")
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "kerSeg_fro")
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "gSeg_net")
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "gSeg_fro")
  
}


