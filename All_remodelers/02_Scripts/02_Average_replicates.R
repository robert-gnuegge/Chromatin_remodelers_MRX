# info --------------------------------------------------------------------
# purpose: calculate mean, standard deviation, and sample size from experiment replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 01/19/24
# last modified: 01/19/24


# read data ---------------------------------------------------------------
cut <- read.table(file = "01_Raw_data/Cutting_replicates.txt", header = TRUE)
resect <- read.table(file = "01_Raw_data/Resection_replicates.txt", header = TRUE)



# calc and save cutting mean, sd, and sample size -------------------------
cut_modes <- aggregate(mean ~ strain + time, data = cut,
                       FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = 0))

# sort
cut_modes <- cut_modes[order(cut_modes$strain, cut_modes$time),]

# mean and sd are in a single column (as a matrix)
# convert into separate columns
cut_modes <- data.frame(cut_modes[, 1:(ncol(cut_modes)-1)], as.data.frame(cut_modes$mean))

# add sample size per strain and time point
for(i in 1:nrow(cut_modes)){
  cut_modes$n[i] <- length(cut$experiment_ID[cut$strain == cut_modes$strain[i] & cut$time == cut_modes$time[i]])
}

# save
write.table(x = cut_modes, file = "03_Processed_data/Cutting_averages.txt", row.names = FALSE, col.names = TRUE)



# calc and save resection mean, sd, and sample size -----------------------
resect_modes <- aggregate(mean ~ strain + primers + time, data = resect,
                          FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = 0))

# sort
resect_modes <- resect_modes[order(resect_modes$strain, resect_modes$primers,resect_modes$time),]

# mean and sd are in a single column (as a matrix)
# convert into separate columns
resect_modes <- data.frame(resect_modes[, 1:(ncol(resect_modes)-1)], as.data.frame(resect_modes$mean))

# add sample size per strain and time point
for(i in 1:nrow(resect_modes)){
  resect_modes$n[i] <- length(resect$experiment_ID[resect$strain == resect_modes$strain[i] & resect$primers == resect_modes$primers[i] & resect$time == resect_modes$time[i]])
}

# save
write.table(x = resect_modes, file = "03_Processed_data/Resection_averages.txt", row.names = FALSE, col.names = TRUE)
