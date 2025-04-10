###########
# Reading for node=50 

load("results/sce6b_50.RData")

num_thresholds <- dim(sce_50[[1]])[2]
summary_matrix <- matrix(NA, nrow = length(sce_50), ncol = num_thresholds)

for (i in 1:length(sce_50)) {
  # ADJUST FOR DESIRED STAT
  stat_matrix <- sce_50[[i]][, , 4] 
  summary_matrix[i, ] <- colMeans(stat_matrix, na.rm = TRUE)
}

rownames(summary_matrix) <- c("Greedy", "Greedy1", "Greedy2", 
                              "Narrowest", "Narrowest1", "Narrowest2")
colnames(summary_matrix) <- paste0("Threshold_", 1:num_thresholds)
summary_matrix


##########
# Reading for node=100 

load("results/sce6_100.RData")

num_thresholds <-  dim(sce_100[[1]])[2]
summary_matrix <- matrix(NA, nrow = length(sce_100), ncol = num_thresholds)

for (i in 1:length(sce_100)) {
  # ADJUST FOR DESIRED STAT
  stat_matrix <- sce_100[[i]][, , 4] 
  summary_matrix[i, ] <- colMeans(stat_matrix, na.rm = TRUE)
}

rownames(summary_matrix) <- c("Greedy", "Greedy1", "Greedy2", 
                              "Narrowest", "Narrowest1", "Narrowest2")
colnames(summary_matrix) <- paste0("Threshold_", 1:num_thresholds)
summary_matrix




