# info --------------------------------------------------------------------
# purpose: plot average S1-seq spreading over time
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 02/06/26
# last modified: 02/06/26

plot_dir <- "04_Plots/Avg_S1-seq_spreading"

# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(plotrix)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

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

load(file = "../../Rep1/S1-seq/03_Processed_data/S1-seq_coverage/LSY5935_S1-seq.RData")
WT_1h_S1_seq <- process_S1_seq(GRanges = LSY5935_1h_S1_seq, roi = roi)
WT_2h_S1_seq <- process_S1_seq(GRanges = LSY5935_2h_S1_seq, roi = roi)
WT_4h_S1_seq <- process_S1_seq(GRanges = LSY5935_4h_S1_seq, roi = roi)

load(file = "../../Rep1/S1-seq/03_Processed_data/S1-seq_coverage/LSY5934_S1-seq.RData")
fun30_1h_S1_seq <- process_S1_seq(GRanges = LSY5934_1h_S1_seq, roi = roi)
fun30_2h_S1_seq <- process_S1_seq(GRanges = LSY5934_2h_S1_seq, roi = roi)
fun30_4h_S1_seq <- process_S1_seq(GRanges = LSY5934_4h_S1_seq, roi = roi)

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


# plot all time points WT ======================================================
x <- WT_1h_agg$distance_to_DSB
my_colors <- gray(level = c(0, 0.5, 0.75))

range(c(WT_1h_agg$score, WT_2h_agg$score, WT_4h_agg$score))
range(c(fun30_1h_agg$score, fun30_2h_agg$score, fun30_4h_agg$score))

pdf(file = "tmp.pdf", width=2.25, height=2)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.7, 1.7, 4, 2), tcl = -0.3, mgp = c(2.5, 0.5, 0), las = 1)

# start empty plot
plot(x = NA, y = NA, xlim = c(0, 1300), ylim = c(0, 4500), ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
axis(side = 1, at = c(0, 4, 8, 12) * 100) 
axis(side = 2, at = 0:4 * 1000)

# add data points
points(x = WT_4h_agg$distance_to_DSB, y = WT_4h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[3])
points(x = WT_2h_agg$distance_to_DSB, y = WT_2h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[2])
points(x = WT_1h_agg$distance_to_DSB, y = WT_1h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[1])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_mre11-nd.pdf"))



# plot zoom-in ------------------------------------------------------------
pdf(file = "tmp.pdf", width=1.4, height=1.3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.9, 2.6, 3.9, 2), tcl = -0.2, mgp = c(2.5, 0.3, 0), las = 1)

# start empty plot
plot(x = NA, y = NA, xlim = c(0, 1300), ylim = c(0, 0.8), ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
axis(side = 1, at = 0:1 * 500)
axis(side = 1, at = 1000)
axis(side = 2, at = 0:4 * 0.2)

# add data points
points(x = WT_4h_agg$distance_to_DSB, y = WT_4h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[3])
points(x = WT_2h_agg$distance_to_DSB, y = WT_2h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[2])
points(x = WT_1h_agg$distance_to_DSB, y = WT_1h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[1])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_mre11-nd_zoom.pdf"))


# plot all time points fun30 ======================================================
x <- fun30_1h_agg$distance_to_DSB
my_colors <- gray(level = c(0, 0.5, 0.75))

range(c(fun30_1h_agg$score, fun30_2h_agg$score, fun30_4h_agg$score))
range(c(fun30_1h_agg$score, fun30_2h_agg$score, fun30_4h_agg$score))

pdf(file = "tmp.pdf", width=2.25, height=2)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.7, 1.7, 4, 2), tcl = -0.3, mgp = c(2.5, 0.5, 0), las = 1)

# start empty plot
plot(x = NA, y = NA, xlim = c(0, 1300), ylim = c(0, 4500), ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
axis(side = 1, at = c(0, 4, 8, 12) * 100) # , labels = NA
axis(side = 2, at = 0:4 * 1000)

# add data points
points(x = fun30_4h_agg$distance_to_DSB, y = fun30_4h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[3])
points(x = fun30_2h_agg$distance_to_DSB, y = fun30_2h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[2])
points(x = fun30_1h_agg$distance_to_DSB, y = fun30_1h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[1])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_mre11-nd_fun30.pdf"))



# plot zoom-in ------------------------------------------------------------
pdf(file = "tmp.pdf", width=1.4, height=1.3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.9, 2.6, 3.9, 2), tcl = -0.2, mgp = c(2.5, 0.3, 0), las = 1)

# start empty plot
plot(x = NA, y = NA, xlim = c(0, 1300), ylim = c(0, 0.8), ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
axis(side = 1, at = 0:1 * 500)
axis(side = 1, at = 1000)
axis(side = 2, at = 0:4 * 0.2)

# add data points
points(x = fun30_4h_agg$distance_to_DSB, y = fun30_4h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[3])
points(x = fun30_2h_agg$distance_to_DSB, y = fun30_2h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[2])
points(x = fun30_1h_agg$distance_to_DSB, y = fun30_1h_agg$mov_med, type = "l", lwd = 1.5, col = my_colors[1])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Avg_S1-seq_spreading_mre11-nd_fun30_zoom.pdf"))



# legend ====================================================================

pdf(file = "tmp.pdf", width=0.8, height=0.85)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = paste0(c(1, 2, 4), " h"), col = my_colors, lwd = 1.5, seg.len = 1)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Legend.pdf"))

pdf(file = "tmp.pdf", width=2.2, height=0.45)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = paste0(c(1, 2, 4), " h"), col = my_colors, lwd = 1.5, ncol = 3, 
       seg.len = 1, # adjust line lengths 
       x.intersp = 0.75, # adjust space between symbol and text
       text.width = 1.25 * strwidth("4 h") # adjust space for text
       )
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Legend_horiz.pdf"))
