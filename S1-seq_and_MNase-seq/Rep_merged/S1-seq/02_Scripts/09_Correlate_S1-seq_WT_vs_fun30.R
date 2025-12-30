# info --------------------------------------------------------------------
# purpose: check correlation between WT and fun30
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 12/30/25
# last modified: 12/30/25

# load libraries ----------------------------------------------------------
library(GenomicRanges)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

plot_dir <- "04_Plots/S1-seq_WT_vs_fun30_correlation"
dir.create(path = plot_dir)


# read and process data ---------------------------------------------------

# function to process data
process_data <- function(coverage, roi = DSB_regions(DSBs = SrfIcs, region_width = 4000)){
  sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = coverage, query = roi)), ignore.strand = FALSE)
}

# load and process data
load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h <- process_data(coverage = LSY4518_13B_0h_S1_seq)
LSY4518_13B_1h <- process_data(coverage = LSY4518_13B_1h_S1_seq)
LSY4518_13B_2h <- process_data(coverage = LSY4518_13B_2h_S1_seq)
LSY4518_13B_4h <- process_data(coverage = LSY4518_13B_4h_S1_seq)

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h <- process_data(coverage = LSY5415_0h_S1_seq)
LSY5415_1h <- process_data(coverage = LSY5415_1h_S1_seq)
LSY5415_2h <- process_data(coverage = LSY5415_2h_S1_seq)
LSY5415_4h <- process_data(coverage = LSY5415_4h_S1_seq)

# make sure all GRanges have same order
all(granges(LSY4518_13B_0h) == granges(LSY5415_0h))
all(granges(LSY4518_13B_1h) == granges(LSY5415_1h))
all(granges(LSY4518_13B_2h) == granges(LSY5415_2h))
all(granges(LSY4518_13B_4h) == granges(LSY5415_4h))



# plot correlation --------------------------------------------------------

png(filename = paste0(plot_dir, "/WT_1h_vs_fun30_4h.png"), width = 2.5, height = 2.5, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.2, 1.4, 4.0, 2.0), tcl = -0.3, mgp = c(2, 0.5, 0), las = 1)

# start empty plot
axis_max <-  max(c(LSY4518_13B_1h$score, LSY5415_4h$score))
plot(x = NA, y = NA, xlim = c(0, log10(x = axis_max)), ylim = c(0, log10(x = axis_max)), xaxt = "n", yaxt = "n",
     xlab = expression("log"["10"]*"( WT 1 h )"), ylab = NA)
title(ylab = expression("log"["10"]*"( "*italic("fun30")*Delta~"4 h )"), line = 1.5)
axis(side = 1, at = 0:3)
axis(side = 2, at = 0:3)

# add points
symbols(x = log10(LSY4518_13B_1h$score + 1), y = log10(LSY5415_4h$score + 1), circles = rep(1, length(LSY5415_4h$score)), 
        inches = 0.01, add = TRUE, fg = NA, bg = adjustcolor(col = "black", alpha.f = 0.1))

# calculate cor
cor(x = log10(LSY4518_13B_1h$score + 1), y = log10(LSY5415_4h$score + 1), method = "pearson")
# 0.9182618

# add cor
text(x = par("usr")[1] + 0.9 * (par("usr")[2] - par("usr")[1]), 
     y = par("usr")[3] + 0.1  * (par("usr")[4] - par("usr")[3]), 
     adj = c(1, 0), 
     labels = substitute(italic("r")~"="~0.92))

dev.off()
