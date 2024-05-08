# info --------------------------------------------------------------------
# purpose: check correlation between replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/01/24
# last modified: 05/01/24

# load libraries ----------------------------------------------------------
library(GenomicRanges)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

plot_dir <- "04_Plots/Rep_correlation"
dir.create(path = plot_dir)

# functions ---------------------------------------------------------------

# read and process data
process_data <- function(coverage, roi = DSB_regions(DSBs = SrfIcs, region_width = 4000)){
  sort(x = as_nt_resolved_GRanges(GRanges = subsetByIntersect(subject = coverage, query = roi)), ignore.strand = FALSE)
}

# plotting
add_points_and_corr <- function(rep1, rep2, col, alpha.f = 0.1, x_cor, y_cor){
  # add points
  symbols(x = rep1, y = rep2, circles = rep(1, length(rep1)), inches = 0.01, add = TRUE, fg = NA, bg = adjustcolor(col = col, alpha.f = alpha.f))
  # calculate cor
  r <- round(cor(x = rep1, y = rep2, method = "pearson"), 2)
  # add cor
  text(x = x_cor, y = y_cor, adj = c(1, 0), labels = substitute(italic("r")~"="~r), col = col)
}


# LSY4518-13B correlation -------------------------------------------------

# load and process Rep1 data
load(file = "../../Rep1/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY4518-13B_MNase-seq.RData")
LSY4518_13B_0h_rep1 <- process_data(coverage = LSY4518_13B_0h_MNase_seq)
LSY4518_13B_1h_rep1 <- process_data(coverage = LSY4518_13B_1h_MNase_seq)
LSY4518_13B_2h_rep1 <- process_data(coverage = LSY4518_13B_2h_MNase_seq)
LSY4518_13B_4h_rep1 <- process_data(coverage = LSY4518_13B_4h_MNase_seq)

# load and process Rep2 data
load(file = "../../Rep2/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY4518-13B_MNase-seq.RData")
LSY4518_13B_0h_rep2 <- process_data(coverage = LSY4518_13B_0h_MNase_seq)
LSY4518_13B_1h_rep2 <- process_data(coverage = LSY4518_13B_1h_MNase_seq)
LSY4518_13B_2h_rep2 <- process_data(coverage = LSY4518_13B_2h_MNase_seq)
LSY4518_13B_4h_rep2 <- process_data(coverage = LSY4518_13B_4h_MNase_seq)

# make sure all GRanges have same order
all(granges(LSY4518_13B_0h_rep1) == granges(LSY4518_13B_0h_rep2))
all(granges(LSY4518_13B_1h_rep1) == granges(LSY4518_13B_1h_rep2))
all(granges(LSY4518_13B_2h_rep1) == granges(LSY4518_13B_2h_rep2))
all(granges(LSY4518_13B_4h_rep1) == granges(LSY4518_13B_4h_rep2))

# plot correlation
png(filename = paste0(plot_dir, "/LSY4518-13B.png"), width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.8, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = NA, y = NA, xlim = c(0, 27), ylim = c(0, 27),
     xlab = "Replicate 1", ylab = "Replicate 2")

x_cor <- par("usr")[1] + 0.95 * (par("usr")[2] - par("usr")[1])
y_shift <- 0.085

add_points_and_corr(rep1 = LSY4518_13B_0h_rep1$score, rep2 = LSY4518_13B_0h_rep2$score, 
                    col = JFly_colors[1], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + 0.075 * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY4518_13B_1h_rep1$score, rep2 = LSY4518_13B_1h_rep2$score, 
                    col = JFly_colors[2], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 1 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY4518_13B_2h_rep1$score, rep2 = LSY4518_13B_2h_rep2$score, 
                    col = JFly_colors[3], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 2 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY4518_13B_4h_rep1$score, rep2 = LSY4518_13B_4h_rep2$score, 
                    col = JFly_colors[4], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 3 * y_shift) * (par("usr")[4]-par("usr")[3]))

leg <- legend(x = "topleft", inset = c(-0.025, 0.03), legend = paste0(c(0, 1, 2, 4), " h"), text.col = JFly_colors[1:4], bty = "n")
leg_box_w <- 6
leg_box_h <- leg$rect$h - 1
rect(xleft = leg$text$x[1] - 0.2 * leg_box_w, xright = leg$text$x[1] + 0.8 * leg_box_w, 
     ybottom = mean(range(leg$text$y)) - 0.5 * leg_box_h, ytop = mean(range(leg$text$y)) + 0.5 * leg_box_h)

dev.off()


# LSY5415 correlation -------------------------------------------------

# load and process Rep1 data
load(file = "../../Rep1/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq.RData")
LSY5415_0h_rep1 <- process_data(coverage = LSY5415_0h_MNase_seq)
LSY5415_1h_rep1 <- process_data(coverage = LSY5415_1h_MNase_seq)
LSY5415_2h_rep1 <- process_data(coverage = LSY5415_2h_MNase_seq)
LSY5415_4h_rep1 <- process_data(coverage = LSY5415_4h_MNase_seq)

# load and process Rep2 data
load(file = "../../Rep2/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq.RData")
LSY5415_0h_rep2 <- process_data(coverage = LSY5415_0h_MNase_seq)
LSY5415_1h_rep2 <- process_data(coverage = LSY5415_1h_MNase_seq)
LSY5415_2h_rep2 <- process_data(coverage = LSY5415_2h_MNase_seq)
LSY5415_4h_rep2 <- process_data(coverage = LSY5415_4h_MNase_seq)

# make sure all GRanges have same order
all(granges(LSY5415_0h_rep1) == granges(LSY5415_0h_rep2))
all(granges(LSY5415_1h_rep1) == granges(LSY5415_1h_rep2))
all(granges(LSY5415_2h_rep1) == granges(LSY5415_2h_rep2))
all(granges(LSY5415_4h_rep1) == granges(LSY5415_4h_rep2))

# plot correlation
png(filename = paste0(plot_dir, "/LSY5415.png"), width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.8, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = NA, y = NA, xlim = c(0, 27), ylim = c(0, 27),
     xlab = "Replicate 1", ylab = "Replicate 2")

x_cor <- par("usr")[1] + 0.95 * (par("usr")[2] - par("usr")[1])
y_shift <- 0.085

add_points_and_corr(rep1 = LSY5415_0h_rep1$score, rep2 = LSY5415_0h_rep2$score, 
                    col = JFly_colors[1], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + 0.075 * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5415_1h_rep1$score, rep2 = LSY5415_1h_rep2$score, 
                    col = JFly_colors[2], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 1 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5415_2h_rep1$score, rep2 = LSY5415_2h_rep2$score, 
                    col = JFly_colors[3], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 2 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5415_4h_rep1$score, rep2 = LSY5415_4h_rep2$score, 
                    col = JFly_colors[4], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 3 * y_shift) * (par("usr")[4]-par("usr")[3]))

leg <- legend(x = "topleft", inset = c(-0.025, 0.03), legend = paste0(c(0, 1, 2, 4), " h"), text.col = JFly_colors[1:4], bty = "n")
leg_box_w <- 6
leg_box_h <- leg$rect$h - 1
rect(xleft = leg$text$x[1] - 0.2 * leg_box_w, xright = leg$text$x[1] + 0.8 * leg_box_w, 
     ybottom = mean(range(leg$text$y)) - 0.5 * leg_box_h, ytop = mean(range(leg$text$y)) + 0.5 * leg_box_h)

dev.off()


# LSY5934 correlation -------------------------------------------------

# load and process Rep1 data
load(file = "../../Rep1/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5934_MNase-seq.RData")
LSY5934_0h_rep1 <- process_data(coverage = LSY5934_0h_MNase_seq)
LSY5934_1h_rep1 <- process_data(coverage = LSY5934_1h_MNase_seq)
LSY5934_2h_rep1 <- process_data(coverage = LSY5934_2h_MNase_seq)
LSY5934_4h_rep1 <- process_data(coverage = LSY5934_4h_MNase_seq)

# load and process Rep2 data
load(file = "../../Rep2/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5934_MNase-seq.RData")
LSY5934_0h_rep2 <- process_data(coverage = LSY5934_0h_MNase_seq)
LSY5934_1h_rep2 <- process_data(coverage = LSY5934_1h_MNase_seq)
LSY5934_2h_rep2 <- process_data(coverage = LSY5934_2h_MNase_seq)
LSY5934_4h_rep2 <- process_data(coverage = LSY5934_4h_MNase_seq)

# make sure all GRanges have same order
all(granges(LSY5934_0h_rep1) == granges(LSY5934_0h_rep2))
all(granges(LSY5934_1h_rep1) == granges(LSY5934_1h_rep2))
all(granges(LSY5934_2h_rep1) == granges(LSY5934_2h_rep2))
all(granges(LSY5934_4h_rep1) == granges(LSY5934_4h_rep2))

# plot correlation
png(filename = paste0(plot_dir, "/LSY5934.png"), width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.8, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = NA, y = NA, xlim = c(0, 27), ylim = c(0, 27),
     xlab = "Replicate 1", ylab = "Replicate 2")

x_cor <- par("usr")[1] + 0.95 * (par("usr")[2] - par("usr")[1])
y_shift <- 0.085

add_points_and_corr(rep1 = LSY5934_0h_rep1$score, rep2 = LSY5934_0h_rep2$score, 
                    col = JFly_colors[1], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + 0.075 * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5934_1h_rep1$score, rep2 = LSY5934_1h_rep2$score, 
                    col = JFly_colors[2], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 1 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5934_2h_rep1$score, rep2 = LSY5934_2h_rep2$score, 
                    col = JFly_colors[3], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 2 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5934_4h_rep1$score, rep2 = LSY5934_4h_rep2$score, 
                    col = JFly_colors[4], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 3 * y_shift) * (par("usr")[4]-par("usr")[3]))

leg <- legend(x = "topleft", inset = c(-0.025, 0.03), legend = paste0(c(0, 1, 2, 4), " h"), text.col = JFly_colors[1:4], bty = "n")
leg_box_w <- 6
leg_box_h <- leg$rect$h - 1
rect(xleft = leg$text$x[1] - 0.2 * leg_box_w, xright = leg$text$x[1] + 0.8 * leg_box_w, 
     ybottom = mean(range(leg$text$y)) - 0.5 * leg_box_h, ytop = mean(range(leg$text$y)) + 0.5 * leg_box_h)

dev.off()

# LSY5935 correlation -------------------------------------------------

# load and process Rep1 data
load(file = "../../Rep1/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5935_MNase-seq.RData")
LSY5935_0h_rep1 <- process_data(coverage = LSY5935_0h_MNase_seq)
LSY5935_1h_rep1 <- process_data(coverage = LSY5935_1h_MNase_seq)
LSY5935_2h_rep1 <- process_data(coverage = LSY5935_2h_MNase_seq)
LSY5935_4h_rep1 <- process_data(coverage = LSY5935_4h_MNase_seq)

# load and process Rep2 data
load(file = "../../Rep2/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5935_MNase-seq.RData")
LSY5935_0h_rep2 <- process_data(coverage = LSY5935_0h_MNase_seq)
LSY5935_1h_rep2 <- process_data(coverage = LSY5935_1h_MNase_seq)
LSY5935_2h_rep2 <- process_data(coverage = LSY5935_2h_MNase_seq)
LSY5935_4h_rep2 <- process_data(coverage = LSY5935_4h_MNase_seq)

# make sure all GRanges have same order
all(granges(LSY5935_0h_rep1) == granges(LSY5935_0h_rep2))
all(granges(LSY5935_1h_rep1) == granges(LSY5935_1h_rep2))
all(granges(LSY5935_2h_rep1) == granges(LSY5935_2h_rep2))
all(granges(LSY5935_4h_rep1) == granges(LSY5935_4h_rep2))

# plot correlation
png(filename = paste0(plot_dir, "/LSY5935.png"), width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.8, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = NA, y = NA, xlim = c(0, 27), ylim = c(0, 27),
     xlab = "Replicate 1", ylab = "Replicate 2")

x_cor <- par("usr")[1] + 0.95 * (par("usr")[2] - par("usr")[1])
y_shift <- 0.085

add_points_and_corr(rep1 = LSY5935_0h_rep1$score, rep2 = LSY5935_0h_rep2$score, 
                    col = JFly_colors[1], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + 0.075 * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5935_1h_rep1$score, rep2 = LSY5935_1h_rep2$score, 
                    col = JFly_colors[2], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 1 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5935_2h_rep1$score, rep2 = LSY5935_2h_rep2$score, 
                    col = JFly_colors[3], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 2 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = LSY5935_4h_rep1$score, rep2 = LSY5935_4h_rep2$score, 
                    col = JFly_colors[4], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 3 * y_shift) * (par("usr")[4]-par("usr")[3]))

leg <- legend(x = "topleft", inset = c(-0.025, 0.03), legend = paste0(c(0, 1, 2, 4), " h"), text.col = JFly_colors[1:4], bty = "n")
leg_box_w <- 6
leg_box_h <- leg$rect$h - 1
rect(xleft = leg$text$x[1] - 0.2 * leg_box_w, xright = leg$text$x[1] + 0.8 * leg_box_w, 
     ybottom = mean(range(leg$text$y)) - 0.5 * leg_box_h, ytop = mean(range(leg$text$y)) + 0.5 * leg_box_h)

dev.off()
