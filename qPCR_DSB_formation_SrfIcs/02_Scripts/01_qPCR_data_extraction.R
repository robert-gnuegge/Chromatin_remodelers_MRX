# info --------------------------------------------------------------------
# purpose: extract and process qPCR data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/29/24
# last modified: 03/29/24

save_dir <- "03_Processed_data/"
efficiencies <- read.table(file = "../Src/primer_efficiencies.txt", header = TRUE)

# load libraries ----------------------------------------------------------
library(errors)  # for error propagation

# function definitions ====================================================
CalcCutFraction <- function(Cq_cut_0, Cq_cut_t, Cq_ref_0, Cq_ref_t, cut_primer_eff, ref_primer_eff){
  1 - cut_primer_eff^(Cq_cut_0 - Cq_cut_t) / ref_primer_eff^(Cq_ref_0 - Cq_ref_t)
}

# read, check, and summarize Cq values ========================================

# define plate layout ---------------------------------------------------------
qPCR <- data.frame(plate = c(rep(1, 384), rep(2, 336), rep(3, 336)), 
                   well = c(paste0(rep(LETTERS[1:16], rep(24, 16)), formatC(x = 1:24, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[1:16], rep(21, 16)), formatC(x = 2:22, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[1:16], rep(21, 16)), formatC(x = 2:22, width = 2, format = "d", flag = "0"))),
                   strain = c(rep(c("5934_1", "5935_1", "5934_2", "5935_2"), rep(96, 4)),
                              rep(c("5934_1", "5935_1", "5934_2", "5935_2"), rep(84, 4)),
                              rep(c("5934_1", "5935_1", "5934_2", "5935_2"), rep(84, 4))),
                   time = c(rep(rep(c(0, 1, 2, 4), rep(24, 4)), 4),
                            rep(rep(c(0, 1, 2, 4), rep(21, 4)), 4),
                            rep(rep(c(0, 1, 2, 4), rep(21, 4)), 4)),
                   primers = c(rep(rep(c("oRG54_oRG55", "oRG341_oRG342", "oRG345_oRG346", "oRG347_oRG348", "oRG353_oRG354", "oRG355_oRG356", "oRG361_oRG362", "oRG1075_oRG1076"), rep(3, 8)), 16),
                               rep(rep(c("oRG54_oRG55", "oRG1079_oRG1080", "oRG1081_oRG1082", "oRG1085_oRG1086", "oRG1091_oRG1092", "oRG1093_oRG1094", "oRG1099_oRG1100"), rep(3, 7)), 16),
                               rep(rep(c("oRG54_oRG55", "oRG1101_oRG1102", "oRG1107_oRG1108", "oRG1109_oRG1110", "oRG1115_oRG1116", "oRG1119_oRG1120", "oRG1121_oRG1122"), rep(3, 7)), 16)),
                   Cq = 0)

head(qPCR)
tail(qPCR)

# read Cq Results file --------------------------------------------------------
file_path <- grep(pattern = "Cq", x = list.files(path = "01_Raw_data", full.names = TRUE, recursive = TRUE), value = TRUE)

# save Cq values in qPCR data frame
for (plate in 1:length(file_path)){
  
  tmp <- read.table(file = file_path[plate], sep = "\t", header = TRUE)
  
  for (well in qPCR$well[qPCR$plate == plate]){
    qPCR$Cq[qPCR$well == well & qPCR$plate == plate] <- tmp$Cq[tmp$Well == well & tmp$Fluor == "SYBR"]
  }
  
}


head(qPCR)
tail(qPCR)

# calc Cq mean, sd, and cv ----------------------------------------------------
Cq_modes <- aggregate(Cq ~ strain + time + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                                          )
                      )

# sort
Cq_modes <- Cq_modes[order(Cq_modes$strain, Cq_modes$primers, Cq_modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq_modes <- data.frame(Cq_modes[, 1:(ncol(Cq_modes)-1)], as.data.frame(Cq_modes$Cq))

# what is the largest cv?
max(Cq_modes$cv, na.rm = TRUE)
# 0.05460837

# which are the samples with large cv?
Cq_modes[Cq_modes$cv > 0.025, ]
#     strain time         primers     mean        sd         cv
# 145 5934_1    0 oRG1109_oRG1110 21.21460 0.6193085 0.02919256
# 161 5934_1    0 oRG1115_oRG1116 20.69404 0.5622304 0.02716872
# 193 5934_1    0 oRG1121_oRG1122 21.34393 1.1655574 0.05460837

# remove outliers
qPCR[qPCR$strain == "5934_1" & qPCR$time == 0 & qPCR$primers == "oRG1109_oRG1110", ]
qPCR$Cq[qPCR$plate == 3 & qPCR$well == "A11"] <- NA

qPCR[qPCR$strain == "5934_1" & qPCR$time == 0 & qPCR$primers == "oRG1115_oRG1116", ]
qPCR$Cq[qPCR$plate == 3 & qPCR$well == "A15"] <- NA

qPCR[qPCR$strain == "5934_1" & qPCR$time == 0 & qPCR$primers == "oRG1121_oRG1122", ]
qPCR$Cq[qPCR$plate == 3 & qPCR$well == "A22"] <- NA


# re-calc Cq mean, sd, and cv ---------------------------------------------
Cq_modes <- aggregate(Cq ~ strain + time + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                      )
)

# sort
Cq_modes <- Cq_modes[order(Cq_modes$strain, Cq_modes$primers, Cq_modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq_modes <- data.frame(Cq_modes[, 1:(ncol(Cq_modes)-1)], as.data.frame(Cq_modes$Cq))

# what is the largest cv?
max(Cq_modes$cv, na.rm = TRUE)
# 0.02103476

# what is the ADH1 fold spread?
tmp <- Cq_modes$mean[Cq_modes$primers == "oRG54_oRG55"]
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 2.209467

# what is the ADH1 fold spread within each strain
aggregate(mean ~ strain, data = Cq_modes[Cq_modes$primers == "oRG54_oRG55", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#   strain     mean
# 1 5934_1 1.330586
# 2 5934_2 1.078481
# 3 5935_1 1.122662
# 4 5935_2 1.271164


# calc cut fraction over time =============================================

# create data frame with Cq means and assigned errors (as attribute)
Cq_mean_sd <- Cq_modes[, c("strain", "time", "primers", "mean")]
errors(Cq_mean_sd$mean) <- Cq_modes$sd

# initialize data.frame to collect results
Cut_fraction <- data.frame()
ref_primers <- "oRG54_oRG55"

# iterate through primers
for (cut_primers in c("oRG341_oRG342", "oRG345_oRG346", "oRG347_oRG348", "oRG353_oRG354", "oRG355_oRG356", 
                      "oRG361_oRG362", "oRG1075_oRG1076", "oRG1079_oRG1080", "oRG1081_oRG1082", "oRG1085_oRG1086", 
                      "oRG1091_oRG1092", "oRG1093_oRG1094", "oRG1099_oRG1100", "oRG1101_oRG1102", "oRG1107_oRG1108", 
                      "oRG1109_oRG1110", "oRG1115_oRG1116", "oRG1119_oRG1120", "oRG1121_oRG1122")){
  
  cat("\nProcessing data for", cut_primers, "...")
  
  # iterate through strains
  for (strain in c("5934_1", "5935_1", "5934_2", "5935_2")){
    
    # iterate through time points
    for (t in c(0, 1, 2, 4)){
      
      # calc cut fraction (with error propagation)
      tmp <- CalcCutFraction(Cq_cut_0 = Cq_mean_sd$mean[Cq_mean_sd$strain == strain & Cq_mean_sd$primers == cut_primers & Cq_mean_sd$time == 0],
                             Cq_cut_t = Cq_mean_sd$mean[Cq_mean_sd$strain == strain & Cq_mean_sd$primers == cut_primers & Cq_mean_sd$time == t],
                             Cq_ref_0 = Cq_mean_sd$mean[Cq_mean_sd$strain == strain & Cq_mean_sd$primers == ref_primers & Cq_mean_sd$time == 0],
                             Cq_ref_t = Cq_mean_sd$mean[Cq_mean_sd$strain == strain & Cq_mean_sd$primers == ref_primers & Cq_mean_sd$time == t],
                             cut_primer_eff = efficiencies$efficiency[efficiencies$primers == cut_primers],
                             ref_primer_eff = efficiencies$efficiency[efficiencies$primers == ref_primers])
      
      # save results
      Cut_fraction <- rbind(Cut_fraction,
                            data.frame(strain = strain, primers = cut_primers, time = t, mean = as.numeric(tmp), sd = errors(tmp)))
      
    }
  }
}

# add HOcs::SrfIcs data
tmp <- read.table(file = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin/FUN30_RAD9/24-02-08-S1-seq_MNase-seq_sample_collection/24-02-12-qPCR_resection_assay/Moments/Cutting.txt", header = TRUE)
tmp <- tmp[tmp$strain %in% c(5934, 5935) & tmp$primers == "oRG46_oRG47", ]
tmp$strain[tmp$strain == "5934"] <- "5934_1"
tmp$strain[tmp$strain == "5935"] <- "5935_1"
Cut_fraction <- rbind(Cut_fraction, tmp)

tmp <- read.table(file = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin/FUN30_RAD9/24-03-26-S1-seq_MNase-seq_sample_collection/24-03-28-qPCR_resection_assay/Moments/Cutting.txt", header = TRUE)
tmp <- tmp[tmp$strain %in% c(5934, 5935) & tmp$primers == "oRG46_oRG47", ]
tmp$strain[tmp$strain == "5934"] <- "5934_2"
tmp$strain[tmp$strain == "5935"] <- "5935_2"
Cut_fraction <- rbind(Cut_fraction, tmp)

# write data to file
write.table(x = Cut_fraction, file = "03_Processed_data/SrfIcs_cutting.txt", row.names = FALSE, col.names = TRUE)
