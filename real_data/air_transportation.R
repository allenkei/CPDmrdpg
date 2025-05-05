
library(rTensor)
source("SBS.R")
source("CUSUM.R")
source("competitor.R")
source("CI.R")



path = "real_data/air_transportation/"


data_whole_22 = read.csv(paste(c(path, "Air Transportation_2022.csv"), collapse =  ""), header = T) 
data_whole_22$MONTH  = 202200 + data_whole_22$MONTH 
data_whole_21 = read.csv(paste(c(path, "Air Transportation_2021.csv"), collapse =  ""), header = T) 
data_whole_21$MONTH  = 202100 + data_whole_21$MONTH
data_whole_20 = read.csv(paste(c(path, "Air Transportation_2020.csv"), collapse =  ""), header = T) 
data_whole_20$MONTH  = 202000 + data_whole_20$MONTH
data_whole_19 = read.csv(paste(c(path, "Air Transportation_2019.csv"), collapse =  ""), header = T) 
data_whole_19$MONTH  = 201900 + data_whole_19$MONTH
data_whole_18 = read.csv(paste(c(path, "Air Transportation_2018.csv"), collapse =  ""), header = T) 
data_whole_18$MONTH  = 201800 + data_whole_18$MONTH
data_whole_17 = read.csv(paste(c(path, "Air Transportation_2017.csv"), collapse =  ""), header = T) 
data_whole_17$MONTH  = 201700 + data_whole_17$MONTH
data_whole_16 = read.csv(paste(c(path, "Air Transportation_2016.csv"), collapse =  ""), header = T) 
data_whole_16$MONTH  = 201600 + data_whole_16$MONTH
data_whole_15 = read.csv(paste(c(path, "Air Transportation_2015.csv"), collapse =  ""), header = T) 
data_whole_15$MONTH  = 201500 + data_whole_15$MONTH

data_whole = rbind(data_whole_22,data_whole_21,data_whole_20,data_whole_19,data_whole_18,data_whole_17,data_whole_16,data_whole_15)
#colnames(data_whole)

AIRLINE_sort = rle(sort(data_whole[,1]))
AIRLINE = AIRLINE_sort$value[order(AIRLINE_sort $lengths, decreasing = T)[1:4]]


#19393 "Southwest Airlines Co.: WN"
#20304  "SkyWest Airlines Inc.: OO"
#19977  "United Air Lines Inc.: UA"
#19790  "Delta Air Lines Inc.: DL"


data_airport_id = read.csv(paste(c(path, "L_AIRPORT_ID.csv"), collapse =  ""), header = T) 


ORIGIN_AIRPORT = rle(sort(data_whole[,2]))
og_id_100 = ORIGIN_AIRPORT$value[order(ORIGIN_AIRPORT$lengths, decreasing = T)[1:50]]

ORIGIN_AIRPORT_names = data_airport_id$Description[data_airport_id$Code %in% og_id_100]


DEST_AIRPORT = rle(sort(data_whole[,4]))
de_id_100 = DEST_AIRPORT$value[order(DEST_AIRPORT$lengths, decreasing = T)[1:50]]
DEST_AIRPORT_names = data_airport_id$Description[data_airport_id$Code %in% de_id_100]

og_id_100_nan = data_airport_id$Description[data_airport_id$Code %in%  og_id_100[which(! og_id_100 %in% de_id_100)]]
de_id_100_nan = data_airport_id$Description[data_airport_id$Code %in%  de_id_100[which(! de_id_100 %in% og_id_100)]]


month_time = sort(unique(data_whole$MONTH))

rm(data_whole_15,data_whole_16,data_whole_17,data_whole_18,data_whole_19,data_whole_20,
   data_whole_21,data_whole_22,DEST_AIRPORT,ORIGIN_AIRPORT,DEST_AIRPORT_names,ORIGIN_AIRPORT_names,
   AIRLINE_sort)


n_1 = length(og_id_100)
n_2 = length(de_id_100)
L =  length(AIRLINE)
TT = length(month_time)
dim = c(TT, n_1, n_2, L)
data_month_product_tensor = array(NA, dim)



for (u in 1:TT){
  
  month = month_time[u]
  
  data_month = data_whole[which(data_whole$MONTH == month),]
  
  for(l in 1:L){
    
    data_month_airline = data_month[which(data_month$AIRLINE_ID == AIRLINE[l]), ]
    
    
    for(i in 1: n_1){
      for(j in 1: n_2){
        x = which( data_month_airline$ORIGIN_AIRPORT_ID == og_id_100[i] & 
                     data_month_airline$DEST_AIRPORT_ID == de_id_100[j] )
        if ( length(x)!=0 ){
          data_month_product_tensor[u, i, j, l] = 1
        }
        else{
          data_month_product_tensor[u, i, j, l] = 0
        }
      }
    }
    
  }
}



dim(data_month_product_tensor) # 91 50 50  4

# save preprocessed data
#save(data_month_product_tensor, file = "real_data/AT_data.RData")




A.tensor <- data_month_product_tensor[1:90,,,]
dim(A.tensor)



# proposed method
num_time = dim(A.tensor)[1]
num_node = dim(A.tensor)[2]
num_layer = dim(A.tensor)[4]


hat.rank <- c(10, 10, num_layer)
threshold_ratio <- c( 0.2, 0.3, 0.4, 0.5)
threshold_list <- rev(threshold_ratio * num_node*sqrt(num_layer)*(log(num_time/2))^(3/2))
intervals <- construct_intervals(num_time/2, sqrt(1/2), 4)




A.tensor.even <- A.tensor[seq(2, num_time, by = 2), , , ]
B.tensor.odd  <- A.tensor[seq(1, num_time-1, by = 2), , , ] 

gains <- cusum_on_intervals(CUSUM_step1, A.tensor.even, verbose = TRUE, intervals, obj.B = B.tensor.odd)
results_g <- seeded_binary_seg(CUSUM_step1, A.tensor.even, num_time/2, CUSUM_res = gains, verbose = FALSE,
                               threshold = threshold_list, method = "Greedy", obj.B = B.tensor.odd)

output <- list()

for (i in 1:(length(results_g)-1)) {
  detected_CP_g <- sort(results_g[[i+1]]$results[, 1])
  detected_CP_gl1 <- refinement1(detected_CP_g, A.tensor.even, B.tensor.odd, hat.rank)
  CI <- construct_CI(alpha=0.05, detected_CP_gl1, A.tensor.even, B.tensor.odd, hat.rank)
  
  output[[i]] <- list()
  output[[i]]$threshold <- results_g[[i+1]]$threshold
  output[[i]]$thres_ratio <- rev(threshold_ratio)[i]
  output[[i]]$detected_CP <- 2*detected_CP_gl1
  output[[i]]$detected_month <- month_time[2*detected_CP_gl1]
  output[[i]]$CI <- CI*2
  
}


# competitor methods
result_kerSeg_net <- Evaluation_kerSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "kerSeg_net")
result_kerSeg_fro <- Evaluation_kerSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "kerSeg_fro")
result_gSeg_net <- Evaluation_gSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "gSeg_net")
result_gSeg_fro <- Evaluation_gSeg(A.tensor, p_threshold=0.01, is_experiment=TRUE, true_CP, "gSeg_fro")


# save the result
kerSeg_net <- month_time[result_kerSeg_net]
kerSeg_fro <- month_time[result_kerSeg_fro]
gSeg_net <- month_time[result_gSeg_net]
gSeg_fro <- month_time[result_gSeg_fro]



# printing
output
kerSeg_net
kerSeg_fro
gSeg_net
gSeg_fro

#save(output, file = paste0("real_data/results/AT_proposed.RData"))
#save(kerSeg_net, file = paste0("real_data/results/AT_kerSeg_net.RData"))
#save(kerSeg_fro, file = paste0("real_data/results/AT_kerSeg_fro.RData"))
#save(gSeg_net, file = paste0("real_data/results/AT_gSeg_net.RData"))
#save(gSeg_fro, file = paste0("real_data/results/AT_gSeg_fro.RData"))

