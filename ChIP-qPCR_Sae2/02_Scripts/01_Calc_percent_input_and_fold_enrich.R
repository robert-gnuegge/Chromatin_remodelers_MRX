# info --------------------------------------------------------------------
# purpose: derive % input from Cq values
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 02/25/25
# last modified: 05/11/25

library(errors)  # for error propagation


# read primer efficiencies ------------------------------------------------
efficiencies <- read.table(file = "../Src/primer_efficiencies.txt", header = TRUE)

# read data ---------------------------------------------------------------
Cq <- read.table(file = "01_Raw_data/Sae2_ChIP-qPCR_Cq_values.txt", header = TRUE)

# average data ------------------------------------------------------------
Cq_modes <- aggregate(Cq ~ strain + sample + primers + time + replicate, data = Cq,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE), 
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)))
# sort
Cq_modes <- Cq_modes[order(Cq_modes$strain, Cq_modes$sample, Cq_modes$primers, Cq_modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq_modes <- data.frame(Cq_modes[, 1:(ncol(Cq_modes)-1)], as.data.frame(Cq_modes$Cq))

# what is the largest cv?
max(Cq_modes$cv, na.rm = TRUE)
# 0.0217753

# which are the samples with large cv?
Cq_modes[Cq_modes$cv > 0.025, ]

# what is the ADH1 fold spread in Input samples?
tmp <- Cq_modes$mean[Cq_modes$primers == "oRG54_oRG55" & Cq_modes$sample == "Input"]
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 4.557841

# what is the ADH1 fold spread within each strain
aggregate(mean ~ strain + replicate, data = Cq_modes[Cq_modes$primers == "oRG54_oRG55" & Cq_modes$sample == "Input", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#    strain replicate     mean
# 1 LSY6097         1 1.668550
# 2 LSY6098         1 2.527822
# 3 LSY6097         2 1.431454
# 4 LSY6098         2 2.121790
# 5 LSY6098         3 3.479497


# calc percent input enrichment ===========================================

# create data frame with Cq means and assigned errors (as attribute)
Cq_err <- Cq_modes[, c("strain", "sample", "time", "primers", "mean", "replicate")]
errors(Cq_err$mean) <- Cq_modes$sd

# define funtion
Calc_percent_input <- function(Cq_Input, Cq_IP, Input_dil, E){
  percent_input <- E^(Cq_Input - log(1 / Input_dil, E) - Cq_IP) * 100
  return(percent_input)
}

# initialize data.frame to collect results
percent_input <- data.frame()

for (replicate in 1:2){
  for (strain in c("LSY6097", "LSY6098")){
    for (primers in c("oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55")){
      for (time in c(0, 1, 2, 4)){
        tmp <- Calc_percent_input(Cq_Input = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                  Cq_IP = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                  Input_dil = 0.01, 
                                  E = efficiencies$efficiency[efficiencies$primers == primers])
        percent_input <- rbind(percent_input,
                               data.frame(strain = strain, primers = primers, time = time, mean = as.numeric(tmp), sd = errors(tmp), replicate = replicate))
      }
    }
  }
}

for (replicate in 3){
  for (strain in c("LSY6098")){
    for (primers in c("oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55")){
      for (time in c(0, 1, 2, 4)){
        tmp <- Calc_percent_input(Cq_Input = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                  Cq_IP = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                  Input_dil = 0.01, 
                                  E = efficiencies$efficiency[efficiencies$primers == primers])
        percent_input <- rbind(percent_input,
                               data.frame(strain = strain, primers = primers, time = time, mean = as.numeric(tmp), sd = errors(tmp), replicate = replicate))
      }
    }
  }
}


# write data to file
write.table(x = percent_input, file = "03_Processed_data/percent_input_replicates.txt", row.names = FALSE, col.names = TRUE)


# calc fold enrichment over control =======================================

# create data frame with Cq means and assigned errors (as attribute)
Cq_err <- Cq_modes[, c("strain", "sample", "time", "primers", "mean", "replicate")]
errors(Cq_err$mean) <- Cq_modes$sd

# define funtion
Calc_fold_over_ctrl <- function(Cq_Input_target, Cq_IP_target, Cq_Input_ctrl, Cq_IP_ctrl, Input_dil, E_target, E_ctrl){
  fold_over_ctrl <- E_target^(Cq_Input_target - log(1 / Input_dil, E_target) - Cq_IP_target) / E_ctrl^(Cq_Input_ctrl - log(1 / Input_dil, E_ctrl) - Cq_IP_ctrl)
  return(fold_over_ctrl)
}

# initialize data.frame to collect results
fold_enrich <- data.frame()

for (replicate in 1:2){
  for (strain in c("LSY6097", "LSY6098")){
    for (primers in c("oRG50_oRG51", "oRG52_oRG53")){
      for (time in c(0, 1, 2, 4)){
        tmp <- Calc_fold_over_ctrl(Cq_Input_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                   Cq_Input_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time & Cq_err$replicate == replicate], 
                                   Cq_IP_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                   Cq_IP_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time & Cq_err$replicate == replicate],
                                   Input_dil = 0.01, 
                                   E_target = efficiencies$efficiency[efficiencies$primers == primers],
                                   E_ctrl = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"])
        fold_enrich <- rbind(fold_enrich,
                             data.frame(strain = strain, primers = primers, time = time, mean = as.numeric(tmp), sd = errors(tmp), replicate = replicate))
      }
    }
  }
}

for (replicate in 3){
  for (strain in c("LSY6098")){
    for (primers in c("oRG50_oRG51", "oRG52_oRG53")){
      for (time in c(0, 1, 2, 4)){
        tmp <- Calc_fold_over_ctrl(Cq_Input_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                   Cq_Input_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time & Cq_err$replicate == replicate], 
                                   Cq_IP_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == primers & Cq_err$time == time & Cq_err$replicate == replicate], 
                                   Cq_IP_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time & Cq_err$replicate == replicate],
                                   Input_dil = 0.01, 
                                   E_target = efficiencies$efficiency[efficiencies$primers == primers],
                                   E_ctrl = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"])
        fold_enrich <- rbind(fold_enrich,
                             data.frame(strain = strain, primers = primers, time = time, mean = as.numeric(tmp), sd = errors(tmp), replicate = replicate))
      }
    }
  }
}

# write data to file
write.table(x = fold_enrich, file = "03_Processed_data/fold_enrich_replicates.txt", row.names = FALSE, col.names = TRUE)
