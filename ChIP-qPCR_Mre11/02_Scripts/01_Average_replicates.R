# info --------------------------------------------------------------------
# purpose: calculate mean, standard deviation, and sample size from experiment replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 04/22/24
# last modified: 04/22/24


# read data ---------------------------------------------------------------
percent_input <- read.table(file = "01_Raw_data/percent_input_replicates.txt", header = TRUE)

# calc and save cutting mean, sd, and sample size -------------------------
percent_input_modes <- aggregate(mean ~ strain + sample + primers + time, data = percent_input,
                                 FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = 0))

# sort
percent_input_modes <- percent_input_modes[order(percent_input_modes$strain, percent_input_modes$sample, percent_input_modes$primers, percent_input_modes$time),]

# mean and sd are in a single column (as a matrix)
# convert into separate columns
percent_input_modes <- data.frame(percent_input_modes[, 1:(ncol(percent_input_modes)-1)], as.data.frame(percent_input_modes$mean))

# add sample size per strain and time point
for(i in 1:nrow(percent_input_modes)){
  percent_input_modes$n[i] <- length(percent_input$replicate[percent_input$strain == percent_input_modes$strain[i] 
                                                             & percent_input$sample == percent_input_modes$sample[i]
                                                             & percent_input$primers == percent_input_modes$primers[i]
                                                             & percent_input$time == percent_input_modes$time[i]])
}

# save
write.table(x = percent_input_modes, file = "03_Processed_data/percent_input_averages.txt", row.names = FALSE, col.names = TRUE)
