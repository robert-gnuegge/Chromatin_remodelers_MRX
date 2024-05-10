# info --------------------------------------------------------------------
# purpose: check correlation between replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/10/24
# last modified: 05/10/24

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
load(file = "../../Rep1/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h_rep1 <- process_data(coverage = LSY4518_13B_0h_S1_seq)
LSY4518_13B_1h_rep1 <- process_data(coverage = LSY4518_13B_1h_S1_seq)
LSY4518_13B_2h_rep1 <- process_data(coverage = LSY4518_13B_2h_S1_seq)
LSY4518_13B_4h_rep1 <- process_data(coverage = LSY4518_13B_4h_S1_seq)

# load and process Rep2 data
load(file = "../../Rep2/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h_rep2 <- process_data(coverage = LSY4518_13B_0h_S1_seq)
LSY4518_13B_1h_rep2 <- process_data(coverage = LSY4518_13B_1h_S1_seq)
LSY4518_13B_2h_rep2 <- process_data(coverage = LSY4518_13B_2h_S1_seq)
LSY4518_13B_4h_rep2 <- process_data(coverage = LSY4518_13B_4h_S1_seq)

# make sure all GRanges have same order
all(granges(LSY4518_13B_0h_rep1) == granges(LSY4518_13B_0h_rep2))
all(granges(LSY4518_13B_1h_rep1) == granges(LSY4518_13B_1h_rep2))
all(granges(LSY4518_13B_2h_rep1) == granges(LSY4518_13B_2h_rep2))
all(granges(LSY4518_13B_4h_rep1) == granges(LSY4518_13B_4h_rep2))

axis_max <-  max(c(LSY4518_13B_0h_rep1$score, LSY4518_13B_0h_rep2$score, 
                   LSY4518_13B_1h_rep1$score, LSY4518_13B_1h_rep2$score, 
                   LSY4518_13B_2h_rep1$score, LSY4518_13B_2h_rep2$score,
                   LSY4518_13B_4h_rep1$score, LSY4518_13B_4h_rep2$score))

# plot correlation
png(filename = paste0(plot_dir, "/LSY4518-13B.png"), width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.6, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = NA, y = NA, xlim = c(0, log10(x = axis_max)), ylim = c(0, log10(x = axis_max)),
     xlab = expression("log"["10"]*"(1"^"st"~"Replicate + 1)"), ylab = expression("log"["10"]*"(2"^"nd"~"Replicate + 1)"))

x_cor <- par("usr")[1] + 0.95 * (par("usr")[2] - par("usr")[1])
y_shift <- 0.085

add_points_and_corr(rep1 = log10(LSY4518_13B_1h_rep1$score + 1), rep2 = log10(LSY4518_13B_1h_rep2$score + 1), 
                    col = JFly_colors[2], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 1 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = log10(LSY4518_13B_2h_rep1$score + 1), rep2 = log10(LSY4518_13B_2h_rep2$score + 1), 
                    col = JFly_colors[3], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 2 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = log10(LSY4518_13B_4h_rep1$score + 1), rep2 = log10(LSY4518_13B_4h_rep2$score + 1), 
                    col = JFly_colors[4], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 3 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = log10(LSY4518_13B_0h_rep1$score + 1), rep2 = log10(LSY4518_13B_0h_rep2$score + 1),
                    col = JFly_colors[1], alpha.f = 0.1,
                    x_cor = x_cor, y_cor = par("usr")[3] + 0.075 * (par("usr")[4]-par("usr")[3]))

leg <- legend(x = "topleft", inset = c(-0.025, 0.03), legend = paste0(c(0, 1, 2, 4), " h"), text.col = JFly_colors[1:4], bty = "n")
leg_box_w <- 0.7
leg_box_h <- leg$rect$h - 0.1
rect(xleft = leg$text$x[1] - 0.25 * leg_box_w, xright = leg$text$x[1] + 0.75 * leg_box_w, 
     ybottom = mean(range(leg$text$y)) - 0.5 * leg_box_h, ytop = mean(range(leg$text$y)) + 0.5 * leg_box_h)

dev.off()


# LSY5415 correlation -----------------------------------------------------

# load and process Rep1 data
load(file = "../../Rep1/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h_rep1 <- process_data(coverage = LSY5415_0h_S1_seq)
LSY5415_1h_rep1 <- process_data(coverage = LSY5415_1h_S1_seq)
LSY5415_2h_rep1 <- process_data(coverage = LSY5415_2h_S1_seq)
LSY5415_4h_rep1 <- process_data(coverage = LSY5415_4h_S1_seq)

# load and process Rep2 data
load(file = "../../Rep2/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h_rep2 <- process_data(coverage = LSY5415_0h_S1_seq)
LSY5415_1h_rep2 <- process_data(coverage = LSY5415_1h_S1_seq)
LSY5415_2h_rep2 <- process_data(coverage = LSY5415_2h_S1_seq)
LSY5415_4h_rep2 <- process_data(coverage = LSY5415_4h_S1_seq)

# make sure all GRanges have same order
all(granges(LSY5415_0h_rep1) == granges(LSY5415_0h_rep2))
all(granges(LSY5415_1h_rep1) == granges(LSY5415_1h_rep2))
all(granges(LSY5415_2h_rep1) == granges(LSY5415_2h_rep2))
all(granges(LSY5415_4h_rep1) == granges(LSY5415_4h_rep2))

axis_max <-  max(c(LSY5415_0h_rep1$score, LSY5415_0h_rep2$score, 
                   LSY5415_1h_rep1$score, LSY5415_1h_rep2$score, 
                   LSY5415_2h_rep1$score, LSY5415_2h_rep2$score,
                   LSY5415_4h_rep1$score, LSY5415_4h_rep2$score))

# plot correlation
png(filename = paste0(plot_dir, "/LSY5415.png"), width = 3, height = 3, units = "in", res = 900)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.8, 0.6, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.5, 0), las = 1)

plot(x = NA, y = NA, xlim = c(0, log10(x = axis_max)), ylim = c(0, log10(x = axis_max)),
     xlab = expression("log"["10"]*"(1"^"st"~"Replicate + 1)"), ylab = expression("log"["10"]*"(2"^"nd"~"Replicate + 1)"))

x_cor <- par("usr")[1] + 0.95 * (par("usr")[2] - par("usr")[1])
y_shift <- 0.085

add_points_and_corr(rep1 = log10(LSY5415_1h_rep1$score + 1), rep2 = log10(LSY5415_1h_rep2$score + 1), 
                    col = JFly_colors[2], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 1 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = log10(LSY5415_2h_rep1$score + 1), rep2 = log10(LSY5415_2h_rep2$score + 1), 
                    col = JFly_colors[3], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 2 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = log10(LSY5415_4h_rep1$score + 1), rep2 = log10(LSY5415_4h_rep2$score + 1), 
                    col = JFly_colors[4], alpha.f = 0.1, 
                    x_cor = x_cor, y_cor = par("usr")[3] + (0.075 + 3 * y_shift) * (par("usr")[4]-par("usr")[3]))

add_points_and_corr(rep1 = log10(LSY5415_0h_rep1$score + 1), rep2 = log10(LSY5415_0h_rep2$score + 1),
                    col = JFly_colors[1], alpha.f = 0.1,
                    x_cor = x_cor, y_cor = par("usr")[3] + 0.075 * (par("usr")[4]-par("usr")[3]))

leg <- legend(x = "topleft", inset = c(-0.025, 0.03), legend = paste0(c(0, 1, 2, 4), " h"), text.col = JFly_colors[1:4], bty = "n")
leg_box_w <- 0.7
leg_box_h <- leg$rect$h - 0.1
rect(xleft = leg$text$x[1] - 0.25 * leg_box_w, xright = leg$text$x[1] + 0.75 * leg_box_w, 
     ybottom = mean(range(leg$text$y)) - 0.5 * leg_box_h, ytop = mean(range(leg$text$y)) + 0.5 * leg_box_h)

dev.off()
