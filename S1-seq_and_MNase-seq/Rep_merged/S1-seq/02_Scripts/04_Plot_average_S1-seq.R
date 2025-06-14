# info --------------------------------------------------------------------
# purpose: plot average S1-seq spreading over time
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 06/01/25
# last modified: 06/14/25


# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(plotrix)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")


# function definitions ----------------------------------------------------

# find x for which AUC is half of total AUC
# argument: numeric vectors
# result: double
find_x_where_half_AUC <- function(x, y){
  AUC <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2  # area under curve (source: https://stackoverflow.com/a/30280873/11705274)
  # find x where half AUC is reached  
  AUC.x <- 0
  z <- 1
  while (AUC.x <= 0.5 * AUC) {
    z <- z + 1
    AUC.x <- sum(diff(x[1:z]) * (head(y[1:z],-1)+tail(y[1:z],-1)))/2
  }
  return(x[z])
}

# load and process S1-seq data --------------------------------------------
SrfIcs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
SrfIcs <- SrfIcs[SrfIcs$DSB_kinetics_rank < 11]  # include only efficiently cut sites
roi <- DSB_regions(DSBs = SrfIcs, region_width = 4000, up_rev_down_fw = TRUE)
# restrict data processing to regions around DSBs
# chose a bigger region size than used for plotting to allow smoothing

process_S1_seq <- function(GRanges, roi){
  tmp <- subsetByIntersect(subject = GRanges, query = roi)  # only keep S1-seq coverage in DSB regions
  tmp <- sort(as_nt_resolved_GRanges(tmp), ignore.strand = TRUE)  # nt resolution, and sort
  # add distance to DSB
  # immediately adjacent positions are also resulting in a 0 distance, therefore 1 is added to all distances and then the SrfIcs positions are set to 0
  tmp$distance_to_DSB <- mcols(distanceToNearest(x = tmp, subject = SrfIcs))$distance + 1
  mcols(tmp[nearest(x = SrfIcs, subject = tmp) - 1])$distance_to_DSB <- 0
  return(tmp)
}

load(file = "03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
WT_1h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_1h_S1_seq, roi = roi)
WT_2h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_2h_S1_seq, roi = roi)
WT_4h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_4h_S1_seq, roi = roi)

load(file = "03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
fun30_1h_S1_seq <- process_S1_seq(GRanges = LSY5415_1h_S1_seq, roi = roi)
fun30_2h_S1_seq <- process_S1_seq(GRanges = LSY5415_2h_S1_seq, roi = roi)
fun30_4h_S1_seq <- process_S1_seq(GRanges = LSY5415_4h_S1_seq, roi = roi)

# calculate average (mean) and moving average as function of distance from DSB and grouped by DSB kinetics
k <- 51
keep <- 2

WT_1h_agg <- aggregate(score ~ distance_to_DSB, data = WT_1h_S1_seq, FUN = mean)
WT_2h_agg <- aggregate(score ~ distance_to_DSB, data = WT_2h_S1_seq, FUN = mean)
WT_4h_agg <- aggregate(score ~ distance_to_DSB, data = WT_4h_S1_seq, FUN = mean)

WT_1h_agg$mov_med <- moving_average(x = WT_1h_agg$score, k = k, keep = keep)
WT_2h_agg$mov_med <- moving_average(x = WT_2h_agg$score, k = k, keep = keep)
WT_4h_agg$mov_med <- moving_average(x = WT_4h_agg$score, k = k, keep = keep)

fun30_1h_agg <- aggregate(score ~ distance_to_DSB, data = fun30_1h_S1_seq, FUN = mean)
fun30_2h_agg <- aggregate(score ~ distance_to_DSB, data = fun30_2h_S1_seq, FUN = mean)
fun30_4h_agg <- aggregate(score ~ distance_to_DSB, data = fun30_4h_S1_seq, FUN = mean)

fun30_1h_agg$mov_med <- moving_average(x = fun30_1h_agg$score, k = k, keep = keep)
fun30_2h_agg$mov_med <- moving_average(x = fun30_2h_agg$score, k = k, keep = keep)
fun30_4h_agg$mov_med <- moving_average(x = fun30_4h_agg$score, k = k, keep = keep)


# plotting ================================================================
plot_dir <- "04_Plots/Avg_S1-seq_spreading"
dir.create(path = plot_dir, recursive = TRUE, showWarnings = FALSE)

my_colors <- JFly_colors[c(1, 2)]

# 1h ----------------------------------------------------------------------
fun30_1h_agg$mov_med[1:3]
fun30_1h_agg$mov_med[1:2] <- fun30_1h_agg$mov_med[1:2] - 170  # to fit into part above y-axis break (see below)

pdf(file = "tmp.pdf", width=2.5, height=2.25)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

# plot fun30 moving average
plot(x = fun30_1h_agg$distance_to_DSB, y = fun30_1h_agg$mov_med, 
     xlim = c(0, 1300), ylim = c(0, 60), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[2], xaxt = "n", yaxt = "n")
title(xlab = "Distance from DSB (nt)", line = 2)

# axes ticks
axis(side = 1, at = c(0, 4, 8, 12) * 100)
axis(side = 2, at = 0:6 * 10, labels = c(0:4 * 10, 220, 230))

# y-axis break
axis.break(axis = 2, breakpos = 45, style = "slash", bgcol = "white")
rect(xleft = -10, ybottom = 45 - 0.5, xright = 10, ytop =  45 + 0.5, col = "white", border = "white")

# # add average spread (half AUC)
# i <- find_x_where_half_AUC(x = fun30_1h_agg$distance_to_DSB, y = fun30_1h_agg$mov_med)
# yy <- par("usr")[4] - 0.15 * (par("usr")[4] - par("usr")[3]) - y_shift
# lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
# text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[2])

# add WT moving average
points(x = WT_1h_agg$distance_to_DSB, y = WT_1h_agg$mov_med, 
       type = "l", lwd = 1.5, col = my_colors[1])

# # add "average" spread (half AUC)
# i <- find_x_where_half_AUC(x = WT_1h_agg$distance_to_DSB, y = WT_1h_agg$mov_med)
# yy <- par("usr")[4] - 0.325 * (par("usr")[4] - par("usr")[3]) - y_shift
# lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
# text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

# legend
legend(x = "right", legend = c("WT", expression(italic("fun30"*Delta))), 
       col = my_colors, lwd = 1.5, seg.len = 1, inset = 0.05)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_1h.pdf"))


# 2h ----------------------------------------------------------------------
pdf(file = "tmp.pdf", width=2.5, height=2.25)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

# plot fun30 moving average
plot(x = fun30_2h_agg$distance_to_DSB, y = fun30_2h_agg$mov_med, 
     xlim = c(0, 1300), ylim = c(0, 60), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[2], xaxt = "n")
title(xlab = "Distance from DSB (nt)", line = 2)

# x axis ticks
axis(side = 1, at = c(0, 4, 8, 12) * 100)

# # add average spread (half AUC)
# i <- find_x_where_half_AUC(x = fun30_2h_agg$distance_to_DSB, y = fun30_2h_agg$mov_med)
# yy <- par("usr")[4] - 0.15 * (par("usr")[4] - par("usr")[3]) - y_shift
# lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
# text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[2])

# add WT moving average
points(x = WT_2h_agg$distance_to_DSB, y = WT_2h_agg$mov_med, 
       type = "l", lwd = 1.5, col = my_colors[1])

# # add "average" spread (half AUC)
# i <- find_x_where_half_AUC(x = WT_2h_agg$distance_to_DSB, y = WT_2h_agg$mov_med)
# yy <- par("usr")[4] - 0.325 * (par("usr")[4] - par("usr")[3]) - y_shift
# lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
# text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_2h.pdf"))


# 4h ----------------------------------------------------------------------
pdf(file = "tmp.pdf", width=2.5, height=2.25)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 2), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)

# plot fun30 moving average
plot(x = fun30_4h_agg$distance_to_DSB, y = fun30_4h_agg$mov_med, 
     xlim = c(0, 1300), ylim = c(0, 60), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[2], xaxt = "n")
title(xlab = "Distance from DSB (nt)", line = 2)

# x axis ticks
axis(side = 1, at = c(0, 4, 8, 12) * 100)

# # add average spread (half AUC)
# i <- find_x_where_half_AUC(x = fun30_4h_agg$distance_to_DSB, y = fun30_4h_agg$mov_med)
# yy <- par("usr")[4] - 0.15 * (par("usr")[4] - par("usr")[3]) - y_shift
# lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
# text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[2])

# add WT moving average
points(x = WT_4h_agg$distance_to_DSB, y = WT_4h_agg$mov_med, 
       type = "l", lwd = 1.5, col = my_colors[1])

# # add "average" spread (half AUC)
# i <- find_x_where_half_AUC(x = WT_4h_agg$distance_to_DSB, y = WT_4h_agg$mov_med)
# yy <- par("usr")[4] - 0.3 * (par("usr")[4] - par("usr")[3]) - y_shift
# lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
# text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_4h.pdf"))
