# info --------------------------------------------------------------------
# purpose: collect resection data for the relevant strains
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 01/24/24
# last modified: 01/24/24


# read data ---------------------------------------------------------------
cut <- read.table(file = "../Misc/01_Raw_data/All_cutting_data.txt", header = TRUE)
resect <- read.table(file = "../Misc/01_Raw_data/All_resection_data.txt", header = TRUE)

# only keep relevant strains
unique(cut$strain)
to_keep <- c("4518-13B", "5415", "4822-3C", "4994-89A", "5758", "5983", "5985", "5986", "5987")
# WT, fun30, sth1-frb, snf2-frb, chd1, ino80-frb, swr1, isw1, isw2
cut <- cut[(cut$strain %in% to_keep), ]
resect <- resect[(resect$strain %in% to_keep), ]

# save collected data -----------------------------------------------------
write.table(x = cut, file = "01_Raw_data/Cutting_replicates.txt", row.names = FALSE, col.names = TRUE)
write.table(x = resect, file = "01_Raw_data/Resection_replicates.txt", row.names = FALSE, col.names = TRUE)
