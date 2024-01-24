# info --------------------------------------------------------------------
# purpose: collect all resection assay data for averaging
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 01/24/24
# last modified: 01/24/24

# Collect cutting and resection data --------------------------------------
file_paths <- grep(pattern = "/Moments/", x = list.files(path = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin", full.names = TRUE, recursive = TRUE), value = TRUE)

# remove unwanted dirs and files
file_paths <- file_paths[!(grepl(pattern = "Nicking_at_HMR", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "Nicking_at_PHO_promoters", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "snf5", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "rsc2", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "sth1-aid_snf2-aid", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "22-12-23-qPCR_resection_assay_nicking", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "22-10-12-pGPD-OsTir1_artifact_resection_assay", x = file_paths))]
file_paths <- file_paths[(grepl(pattern = "Cutting", x = file_paths) | grepl(pattern = "Resection", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "Resection_reduced_RsaI_activity.txt", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "Resection_reduced_RsaI_activity.txt", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "Resection_wt_sth1-aa_repetition", x = file_paths))]
file_paths <- file_paths[!(grepl(pattern = "Resection_wt_sth1_repetition.txt", x = file_paths))]


# how are the file names called?
unique(basename(file_paths))
# [1] "Cutting.txt"        "Resection.txt"


# iterate through dirs and save data in dfs "cut" and "resect"
cut <- data.frame()
resect <- data.frame()
experiment_ID <- 0
for(dir in unique(dirname(file_paths))){
  
    # step experiment_ID counter
    experiment_ID <- experiment_ID + 1
    
    # cutting
    tmp <- read.table(file = paste(dir, "Cutting.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
    tmp <- tmp[!grepl(pattern = "conf", x = colnames(tmp))]  # remove "conf.2.5" and "conf.97.5" columns
    cut <- rbind(cut, cbind(tmp, experiment_ID = experiment_ID))
    
    # resection
    tmp <- read.table(file = paste(dir, "Resection.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE)
    tmp <- tmp[!grepl(pattern = "conf", x = colnames(tmp))]  # remove "conf.2.5" and "conf.97.5" columns
    resect <- rbind(resect, cbind(tmp, experiment_ID = experiment_ID))

}

# sanity check that there are no duplicated rows
all(!duplicated(cut))
all(!duplicated(resect))


# save collected data -----------------------------------------------------
write.table(x = cut, file = "01_Raw_data/All_cutting_data.txt", row.names = FALSE, col.names = TRUE)
write.table(x = resect, file = "01_Raw_data/All_resection_data.txt", row.names = FALSE, col.names = TRUE)
