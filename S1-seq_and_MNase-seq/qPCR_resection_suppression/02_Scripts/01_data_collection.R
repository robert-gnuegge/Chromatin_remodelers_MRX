# info --------------------------------------------------------------------
# purpose: collect qPCR data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/31/24
# last modified: 07/13/25

# collect resection data ------------------------------------------------------------
tmp <- read.table(file = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin/FUN30_RAD9/24-02-08-S1-seq_MNase-seq_sample_collection/24-02-12-qPCR_resection_assay/Moments/Resection.txt", header = TRUE)
tmp$strain <- paste0(tmp$strain, "_1")

tmp_2  <- read.table(file = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin/FUN30_RAD9/24-03-26-S1-seq_MNase-seq_sample_collection/24-03-28-qPCR_resection_assay/Moments/Resection.txt", header = TRUE)
tmp_2$strain <- paste0(tmp_2$strain, "_2")

out <- rbind(tmp, tmp_2)

# write data to file
write.table(x = out, file = "03_Processed_data/Resection.txt", row.names = FALSE, col.names = TRUE)

# collect cutting data ------------------------------------------------------------
tmp <- read.table(file = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin/FUN30_RAD9/24-02-08-S1-seq_MNase-seq_sample_collection/24-02-12-qPCR_resection_assay/Moments/Cutting.txt", header = TRUE)
tmp$strain <- paste0(tmp$strain, "_1")

tmp_2  <- read.table(file = "/home/robert/Research/LabNoteBook/Projects/MRX_nicking_chromatin/FUN30_RAD9/24-03-26-S1-seq_MNase-seq_sample_collection/24-03-28-qPCR_resection_assay/Moments/Cutting.txt", header = TRUE)
tmp_2$strain <- paste0(tmp_2$strain, "_2")

out <- rbind(tmp, tmp_2)

# write data to file
write.table(x = out, file = "03_Processed_data/Cutting.txt", row.names = FALSE, col.names = TRUE)
