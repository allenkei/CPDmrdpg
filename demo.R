
library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("utility.R")
source("simulate_data.R")
library(kneedle)

###################
# generating data #
###################

A.tensor <- array(NA, c(150, 50, 50, 4))
sbm_params <- get_sbm_params(n=50, L=4, n_c=c(4, 4), flip_layer=TRUE)

probability_1 = sbm_params[[1]]
probability_2 = sbm_params[[2]]

for(t_iter in 1:50) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
for(t_iter in 51:100) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_2)
for(t_iter in 101:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)

# for(t_iter in 1:80) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_1)
# for(t_iter in 81:150) A.tensor[t_iter,,,] <- generate_tensor_probability_directed(n_1=50, n_2=50, L=4, probability_2)



##################################
# Example Code (Not CUSUM_step1) # 
##################################

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
results <- seeded_binary_seg(CUSUM_frobenius, A.tensor, 150, CUSUM_res = results_all_Layer, 
                             threshold = c(10, 9, 8), method = "Greedy", rank = hat.rank)

results[[2]]$results

##############################
# Example Code (CUSUM_step1) # 
##############################

A.tensor.even <- A.tensor[seq(2, 150, by = 2), , , ]
B.tensor.odd  <- A.tensor[seq(1, 149, by = 2), , , ]
if (dim(A.tensor.even)[1] != dim(B.tensor.odd)[1]) {
  stop("Make sure even and odd have same length")
}

hat.rank <- c(15, 15, 15)
s <- 0
e <- 75
frobenius_holder <- numeric(74)

for(t in 2:73){
  frobenius_holder[t] <- CUSUM_step1(A.tensor.even, s, e, t, obj.B = B.tensor.odd)
}

plot(1:74, frobenius_holder, type='l')

intervals <- construct_intervals(75, sqrt(1/2), 2)
results_one <- cusum_on_intervals(CUSUM_step1, A.tensor.even, c(30, 38), obj.B = B.tensor.odd)
results_all_step1 <- cusum_on_intervals(CUSUM_step1, A.tensor.even, intervals, obj.B = B.tensor.odd)


# Calculate intervals & CUSUM within
# results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, alpha = 1/2, m = 4,
#                              threshold = 200, method = "Narrowest", obj.B = B.tensor.odd)

# Pass in CUSUM results (as a matrix) to save computation
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = c(1800, 900, 500, 200, 1), method = "Narrowest", obj.B = B.tensor.odd)
results <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = c(1000, 500, 250, 50, 1), method = "Greedy", obj.B = B.tensor.odd)

results[[6]]


###########################
# Example Model Selection # 
###########################

source("utility.R")
init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                          threshold = 0, method = "Greedy", obj.B = B.tensor.odd)
max <-max(init[[1]]$results[, 2])
init[[2]]$results[, 2]

# Thresholding using Max Gain
threshold_list <- seq(1, max + 1, length.out=50)
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                             threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- unique(out[[1]])
ncps <- unique(out[[2]])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Equally Spaced")
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(1:length(ncps), vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Max Gain, log spacing
threshold_list <- exp(seq(1, log(max) + 1, length.out=50))
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- unique(out[[1]])
ncps <- unique(out[[2]])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Log-Spaced")
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(1:length(ncps), vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Changepoints from Greedy, top 10 Gains
threshold_list = c(max + 1, init[[2]]$results[, 2][1:10])
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- out[[1]]
ncps <- out[[2]]
plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
text(ncps, vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(ncps, vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Changepoints from Greedy, top 20 Gains
threshold_list = c(max + 1, init[[2]]$results[, 2][1:20])
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- out[[1]]
ncps <- out[[2]]
plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
text(ncps, vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(ncps, vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Changepoints from Greedy, all
threshold_list = c(max + 1, init[[2]]$results[, 2])
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- out[[1]]
ncps <- out[[2]]
plot(ncps, vals, xlab = "Number of Changepoints", ylab = "Log-Likelihood", main = "Thresholds via Changepoints")
text(ncps, vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(ncps, vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Max Gain, < 15 Changepoints
threshold_list <- seq(1, max + 1, length.out=50)
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- unique(out[[1]][out[[2]] < 15])
ncps <- unique(out[[2]][out[[2]] < 15])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Equally Spaced")
mtext("Less than 15 Changepoints", side = 3, line = 0.5, cex = 0.9)
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(1:length(ncps), vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# Thresholding using Max Gain, log spacing, < 15 Changepoints
threshold_list <- exp(seq(1, log(max) + 1, length.out=50))
results_ms <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
                                threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)
out <- model_selection(results_ms, A.tensor, method = "elbow", hat.rank = hat.rank)[c(4, 5)]
vals <- unique(out[[1]][out[[2]] < 15])
ncps <- unique(out[[2]][out[[2]] < 15])
plot(vals, ylab = "Log-Likelihood", main = "Thresholds Log-Spaced")
mtext("Less than 15 Changepoints", side = 3, line = 0.5, cex = 0.9)
text(1:length(ncps), vals, ncps, pos = 1, cex = 0.8)
knee <- kneedle(1:length(ncps), vals, concave = FALSE, decreasing = FALSE)
print(knee)
symbols(knee[1], knee[2], circles = 1, add = TRUE, inches = 0.08, fg = "red")
text(knee[1], knee[2], "Knee", pos = 3, cex = 0.8, col = "red")

# using utilities/model selection
# Both AIC and BIC prefer fewest (0) changepoints
source("utility.R")
init <- seeded_binary_seg(CUSUM_step1, A.tensor.even, 75, CUSUM_res = results_all_step1, 
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

model_selection(results_ms, A.tensor, method = "AIC", hat.rank = hat.rank)[[1:3]]
# Both AIC and BIC prefer fewest (0) changepoints




