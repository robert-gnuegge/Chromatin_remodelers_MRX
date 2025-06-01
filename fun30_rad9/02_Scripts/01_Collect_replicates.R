# info --------------------------------------------------------------------
# purpose: collect resection data for the relevant strains
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 01/24/24
# last modified: 01/24/24


# read data ---------------------------------------------------------------
cut <- read.table(file = "../Misc/03_Processed_data/All_cutting_data.txt", header = TRUE)
resect <- read.table(file = "../Misc/03_Processed_data/All_resection_data.txt", header = TRUE)

# only keep relevant strains
unique(cut$strain)
to_keep <- c("4518-13B", "5415", "5452", "5451")  # WT, fun30, rad9, fun30 rad9
cut <- cut[(cut$strain %in% to_keep), ]
resect <- resect[(resect$strain %in% to_keep), ]

# save collected data -----------------------------------------------------
write.table(x = cut, file = "01_Raw_data/Cutting_replicates.txt", row.names = FALSE, col.names = TRUE)
write.table(x = resect, file = "01_Raw_data/Resection_replicates.txt", row.names = FALSE, col.names = TRUE)
