
#download.file(url = "https://github.com/allenkei/CPDmrdpg/archive/master.zip", destfile = "master.zip")
#unzip(zipfile = "master.zip")
#setwd("~/RDPGM/CPDmrdpg-main")


library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("eval.R")
source("gen_data.R")
source("competitor.R")


simulate_scenario_competitor <- function(scenario, true_cp, num_node = 50, num_seq = 10, competitor = "kerSeg_net") {
  # competitor: kerSeg_net, kerSeg_fro, gSeg_net, gSeg_fro
  
  temp <- generate(scenario, num_node, 1, FALSE)
  num_time <- dim(temp)[2]; num_layer <- dim(temp)[5]; rm(temp)
  
  
  A.all_seq <- array(NA, c(num_seq, num_time, num_node, num_node, num_layer)) # i.e. 10 sequences empty
  
  for(seq_iter in 1:num_seq) {
    # Generate Data 1-by-1
    set.seed(seq_iter)
    A.all_seq[seq_iter,,,,] <- generate(scenario, num_node, 1, FALSE) # (1,200,50,50,4) # NOTE the first dim = 1
  }
  
  
  if(competitor == "kerSeg_net"){
    results <- Evaluation_kerSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_CP, competitor)
    
  }else if(competitor == "kerSeg_fro"){
    results <- Evaluation_kerSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_CP, competitor)
    
  }else if(competitor == "gSeg_net"){
    results <- Evaluation_gSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_CP, competitor)
    
  }else if(competitor == "gSeg_fro"){
    results <- Evaluation_gSeg(A.all_seq, p_threshold=0.05, is_experiment=FALSE, true_CP, competitor)
    
  }
  
  save(results, file = paste0("results/comp/",scenario,"_n",num_node,"_",competitor,".RData"))
  
  return(results)
}





###########
# Run all #
###########




for (scenario in c("f1", "f2", "f3", "f4", "f5")) {
  
  num_node <- 50
  num_seq <- 10
  
  if (scenario == "f1") {
    true_CP <- c()
  } else if (scenario == "f2") {
    true_CP <- c(20, 60, 80, 160, 180)
  } else if (scenario == "f3") {
    true_CP <- c(50, 100, 150)
  } else if (scenario == "f4") {
    true_CP <- c(20, 60, 80, 160, 180)
  } else if (scenario == "f5") {
    true_CP <- c(50, 100, 150)
  } 
  

  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "kerSeg_net")
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "kerSeg_fro")
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "gSeg_net")
  results <- simulate_scenario_competitor(scenario, true_CP, num_node, num_seq, "gSeg_fro")
  
}
  

