# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

library(errors)  # for error propagation

# create directory for processed data
dir.create(path = "Moments/", showWarnings = FALSE)

# read in helper files and functions
source(file = "/home/robert/Research/Software/R_scripts/Misc_helper_functions.R")

# read primer efficiencies
efficiencies <- read.table(file = "/home/robert/Research/Resources/qPCR_primer_efficiencies/primer_efficiencies.txt", header = TRUE)

# function definitions ====================================================
CalcCutFraction <- function(Cq_cut_0, Cq_cut_t, Cq_ref_0, Cq_ref_t, cut_primer_eff, ref_primer_eff){
  1 - cut_primer_eff^(Cq_cut_0 - Cq_cut_t) / ref_primer_eff^(Cq_ref_0 - Cq_ref_t)
}

CalcResectedFraction <- function(Cq_digest, Cq_mock, cut_fraction, resect_primer_eff){
  2 / ((resect_primer_eff^(Cq_digest - Cq_mock) + 1) * cut_fraction)
}

# read, check, and summarize Cq values ========================================

# define plate layout ---------------------------------------------------------
qPCR <- data.frame(well = c(paste0(rep(LETTERS[1:16], rep(18, 16)), formatC(x = 2:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[1:4], rep(3, 4)), formatC(x = 21:23, width = 2, format = "d", flag = "0"))),
                   sample = c(rep("6072", 72), rep("6073", 72), rep("6074", 72), rep("6075", 72), rep ("NTC", 12)),
                   time = c(rep(rep(c(0, 1, 2, 4), rep(18, 4)), 4), rep(0, 12)),
                   digest = c(rep(c(rep("none", 12), rep("RsaI", 6)), 16), rep("none", 12)),
                   primers = c(rep(rep(c("oRG54_oRG55", "oRG46_oRG47", "oRG50_oRG51", "oRG52_oRG53", "oRG50_oRG51", "oRG52_oRG53"), rep(3, 6)), 16),
                               rep(c("oRG54_oRG55", "oRG46_oRG47", "oRG50_oRG51", "oRG52_oRG53"), rep(3, 4))),
                   Cq = 0)

head(qPCR)
tail(qPCR)

# read Cq Results file --------------------------------------------------------
file_path <- grep(pattern = "Cq", x = list.files(path = "qPCR_data", full.names = TRUE, recursive = TRUE), value = TRUE)

# save Cq values in qPCR data frame
tmp <- read.table(file = file_path, sep = "\t", header = TRUE)
for (well in qPCR$well){
  qPCR$Cq[qPCR$well == well] <- tmp$Cq[tmp$Well == well & tmp$Fluor == "SYBR"]
}

head(qPCR)
tail(qPCR, n = 20)

# check NTCs ------------------------------------------------------------------

# min NTC Cq
min(qPCR$Cq[qPCR$sample == "NTC"], na.rm = TRUE)
# 32.67994
qPCR[qPCR$sample == "NTC", ]

# max non-NTC Cq
max(qPCR$Cq[qPCR$sample != "NTC"], na.rm = TRUE)
# 27.11356

# check non-NTC samples with high Cq values
qPCR[qPCR$Cq > 25 & qPCR$sample != "NTC",]
# RsaI digests

# calc Cq mean, sd, and cv ----------------------------------------------------
Cq_modes <- aggregate(Cq ~ sample + digest + time + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                                          )
                      )

# sort
Cq_modes <- Cq_modes[order(Cq_modes$sample, Cq_modes$primers, Cq_modes$digest, Cq_modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq_modes <- data.frame(Cq_modes[, 1:(ncol(Cq_modes)-1)], as.data.frame(Cq_modes$Cq))

# what is the largest cv?
max(Cq_modes$cv, na.rm = TRUE)
# 0.04034165

# which are the samples with large cv?
Cq_modes[Cq_modes$cv > 0.025, ]
#    sample digest time     primers     mean       sd         cv
# 5     NTC   none    0 oRG46_oRG47 36.10458 1.456519 0.04034165
# 55    NTC   none    0 oRG52_oRG53 35.29511 1.176960 0.03334627

# what is the ADH1 fold spread?
tmp <- Cq_modes$mean[Cq_modes$primers == "oRG54_oRG55" & Cq_modes$sample != "NTC"]
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 1.866728

# what is the ADH1 fold spread within each strain
aggregate(mean ~ sample, data = Cq_modes[Cq_modes$primers == "oRG54_oRG55" & Cq_modes$sample != "NTC", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#   sample     mean
# 1   6072 1.223795
# 2   6073 1.331714
# 3   6074 1.233178
# 4   6075 1.307893

# create data frame with Cq means and assigned errors (as attribute)
Cq_mean_sd <- Cq_modes[, c("sample", "time", "primers", "digest","mean")]
errors(Cq_mean_sd$mean) <- Cq_modes$sd


# calc cut fraction over time =============================================

# initialize data.frame to collect results
Cut_fraction <- data.frame()
ref_primers <- "oRG54_oRG55"
cut_primers <- "oRG46_oRG47"

strains <- unique(Cq_modes$sample)
strains <- strains[strains != "NTC"]

# iterate through strains
for (strain in strains){
  # iterate through time points
  for (t in unique(Cq_mean_sd$time[Cq_mean_sd$sample == strain])){
    
    # calc cut fraction (with error propagation)
    tmp <- CalcCutFraction(Cq_cut_0 = Cq_mean_sd$mean[Cq_mean_sd$sample == strain & Cq_mean_sd$primers == cut_primers & Cq_mean_sd$time == 0],
                           Cq_cut_t = Cq_mean_sd$mean[Cq_mean_sd$sample == strain & Cq_mean_sd$primers == cut_primers & Cq_mean_sd$time == t],
                           Cq_ref_0 = Cq_mean_sd$mean[Cq_mean_sd$sample == strain & Cq_mean_sd$primers == ref_primers & Cq_mean_sd$time == 0],
                           Cq_ref_t = Cq_mean_sd$mean[Cq_mean_sd$sample == strain & Cq_mean_sd$primers == ref_primers & Cq_mean_sd$time == t],
                           cut_primer_eff = efficiencies$efficiency[efficiencies$primers == cut_primers],
                           ref_primer_eff = efficiencies$efficiency[efficiencies$primers == ref_primers])
    
    # save results
    Cut_fraction <- rbind(Cut_fraction,
                          data.frame(strain = strain, primers = cut_primers, time = t, mean = as.numeric(tmp), sd = errors(tmp)))
    
  }
}

# write data to file
write.table(x = Cut_fraction, file = "Moments/Cutting.txt", row.names = FALSE, col.names = TRUE)

# calc resected fraction over time ============================================

# initialize data.frame to collect results
Resected_fraction <- data.frame()

# create data frame with cut fraction means and assigned errors (as attribute)
Cut_mean_sd <- Cut_fraction
errors(Cut_mean_sd$mean) <- Cut_fraction$sd[Cut_fraction$primers == "oRG46_oRG47"]

# iterate through strains
for (strain in strains){
  
  # iterate through resection evaluation amplicons
  for (resect_primers in c("oRG50_oRG51", "oRG52_oRG53")){
    
    # iterate through time points
    for (t in unique(Cq_mean_sd$time[Cq_mean_sd$sample == strain])){
      
      # calc resected fraction (with error propagation)
      tmp <- CalcResectedFraction(Cq_digest = Cq_mean_sd$mean[Cq_mean_sd$sample == strain & Cq_mean_sd$primers == resect_primers & Cq_mean_sd$time == t & Cq_mean_sd$digest != "none"], 
                                  Cq_mock = Cq_mean_sd$mean[Cq_mean_sd$sample == strain & Cq_mean_sd$primers == resect_primers & Cq_mean_sd$time == t & Cq_mean_sd$digest == "none"], 
                                  cut_fraction = Cut_mean_sd$mean[Cut_mean_sd$strain == strain & Cut_mean_sd$time == t], 
                                  resect_primer_eff = efficiencies$efficiency[efficiencies$primers == resect_primers])
      
      # save results
      Resected_fraction <- rbind(Resected_fraction, 
                                 data.frame(strain = strain, primers = resect_primers, time = t, mean = as.numeric(tmp), sd = errors(tmp)))
      
    }
    
  }
  
}

# replace Inf with 0 for t = 0
Resected_fraction$mean[Resected_fraction$time == 0 & Resected_fraction$mean == Inf] <- 0  
Resected_fraction$sd[Resected_fraction$time == 0 & Resected_fraction$sd == Inf] <- 0  # replace Inf with 0 for t = 0

# write data to file
write.table(x = Resected_fraction, file = "Moments/Resection.txt", row.names = FALSE, col.names = TRUE)
