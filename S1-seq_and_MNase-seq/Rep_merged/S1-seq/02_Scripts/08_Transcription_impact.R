# info --------------------------------------------------------------------
# purpose: compare S1-seq of transcribed regions that are co-directional or converging with resection direction
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 06/29/25
# last modified: 06/29/25


# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(rtracklayer)
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

# transform data for axis break with scale change -------------------------
# argument: numeric vector
# result: numeric vector
trans <- function(x, threshold, threshold_trans, factor){
  pmin(x, threshold) + factor * pmax(x - threshold_trans, 0)
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


# load RNA-seq data -------------------------------------------------------
# data are from Maya-Miles et al, 2019 (pmid: 31331360)

# there are two biological replicates
roi <- DSB_regions(DSBs = SrfIcs, region_width = 20000)
RNA_seq_1 <- import(con = "../../../Misc/01_Raw_data/Maya-Miles2019_RNA-seq/GSM3567364_w303_rep1.bigwig",
                    which = roi)
RNA_seq_2 <- import(con = "../../../Misc/01_Raw_data/Maya-Miles2019_RNA-seq/GSM3567364_w303_rep1.bigwig",
                    which = roi)
# let's take the average RNA-seq score
RNA_seq_1 <- subsetByIntersect(subject = RNA_seq_1, query = RNA_seq_2)  # to let both data sets have the same granges
RNA_seq_2 <- subsetByIntersect(subject = RNA_seq_2, query = RNA_seq_1)
all(granges(RNA_seq_1) == granges(RNA_seq_2))
RNA_seq <- RNA_seq_1
RNA_seq$score <- apply(X = cbind(RNA_seq_1$score, RNA_seq_2$score), MARGIN = 1, FUN = mean)


# load genome features ----------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and curated with Adjust_chromosomal_features.R

load(file = "../../../Misc/01_Raw_data/S_cerevisiae_chromosomal_features/S_cerevisiae_genome_features.RData")
all_features <- all_features[all_features$type == "ORF"]


# plot transcribed vs untranscribed =======================================
plot_dir <- "04_Plots/Transcription_impact"
dir.create(path = plot_dir, showWarnings = FALSE)

k <- 51  # for moving average
keep <- 2

# plot WT t = 1 h ------------------------------------------------------------

my_colors <- gray(level = c(0, 0.67))
tmp <- WT_1h_S1_seq

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp, FUN = mean)
# use tmp[tmp$DSB_kinetics_rank < 17] to exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.2, 0.9, 4, 1.7), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
  x <- agg$distance_to_DSB[agg$transcribed]
  y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
  plot(x = x, y = y, xlim = c(0, 800), ylim = c(0, 35), ylab = "Average S1-seq (RPM)", xlab = NA, 
       type = "l", lwd = 1.5, col = my_colors[1])
  title(xlab = "Distance from DSB (nt)", line = 1.75)
  
  y_shift <- 5
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])
  
  x <- agg$distance_to_DSB[!agg$transcribed]
  y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
  points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])
  
  y <- y[x > keep]
  x <- x[x > keep]
  i <- find_x_where_half_AUC(x = x, y = y)
  yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
  lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
  text(x = x[i], y = yy, labels = x[i], adj = c(0.05, -0.3), col = my_colors[2])
  
    # legend(x = "topright", legend = c("Transcribed", "Untranscribed"), col = my_colors, lwd = 1.5, bty = "n")
  
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/WT_1h.pdf"))


# plot WT t = 2 h ------------------------------------------------------------

my_colors <- gray(level = c(0, 0.67))
tmp <- WT_2h_S1_seq

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp, FUN = mean)
# use tmp[tmp$DSB_kinetics_rank < 17] to exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 1.5), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, xlim = c(0, 1000), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[1])
title(xlab = "Distance from DSB (nt)", line = 2)

y_shift <- 6
y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.1, -0.3), col = my_colors[2])

# legend(x = "topright", legend = c("Transcribed", "Untranscribed"), col = my_colors, lwd = 1.5, bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/WT_2h.pdf"))


# plot WT t = 4 h ------------------------------------------------------------

my_colors <- gray(level = c(0, 0.67))
tmp <- WT_4h_S1_seq

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp, FUN = mean)
# use tmp[tmp$DSB_kinetics_rank < 17] to exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 1.5), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, xlim = c(0, 1000), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[1])
title(xlab = "Distance from DSB (nt)", line = 2)

y_shift <- 6
y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.1, -0.3), col = my_colors[2])

# legend(x = "topright", legend = c("Transcribed", "Untranscribed"), col = my_colors, lwd = 1.5, bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/WT_4h.pdf"))


# plot fun30 t = 1 h ------------------------------------------------------------

my_colors <- c(JFly_colors[2], adjustcolor(col = JFly_colors[2], alpha.f = 0.67))
tmp <- fun30_1h_S1_seq

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp, FUN = mean)
# use tmp[tmp$DSB_kinetics_rank < 17] to exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 1.5), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, xlim = c(0, 1000), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[1])
title(xlab = "Distance from DSB (nt)", line = 2)

y_shift <- 6
y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.1, -0.3), col = my_colors[2])

# legend(x = "topright", legend = c("Transcribed", "Untranscribed"), col = my_colors, lwd = 1.5, bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/fun30_1h.pdf"))


# plot fun30 t = 2 h ------------------------------------------------------------

my_colors <- c(JFly_colors[2], adjustcolor(col = JFly_colors[2], alpha.f = 0.67))
tmp <- fun30_2h_S1_seq

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp, FUN = mean)
# use tmp[tmp$DSB_kinetics_rank < 17] to exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.7, 4, 1.5), tcl = -0.3, mgp = c(2.5, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, xlim = c(0, 1000), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[1])
title(xlab = "Distance from DSB (nt)", line = 2)

y_shift <- 6
y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.1, -0.3), col = my_colors[2])

# legend(x = "topright", legend = c("Transcribed", "Untranscribed"), col = my_colors, lwd = 1.5, bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/fun30_2h.pdf"))


# plot fun30 t = 4 h ------------------------------------------------------------

my_colors <- c(JFly_colors[2], adjustcolor(col = JFly_colors[2], alpha.f = 0.5))
tmp <- fun30_4h_S1_seq

hits <- findOverlaps(subject = RNA_seq, query = tmp, ignore.strand = TRUE)
tmp$RNA_seq <- mcols(RNA_seq[subjectHits(hits)])$score
tmp$transcribed <- tmp$RNA_seq > 0

agg <- aggregate(score ~ distance_to_DSB + transcribed, data = tmp, FUN = mean)
# use tmp[tmp$DSB_kinetics_rank < 17] to exclude slowly formed DSBs

pdf(file = "tmp.pdf", width=2.75, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.2, 0.9, 4, 1.7), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
x <- agg$distance_to_DSB[agg$transcribed]
y <- moving_average(x = agg$score[agg$transcribed], k = k, keep = keep)
plot(x = x, y = y, xlim = c(0, 800), ylim = c(0, 35), ylab = "Average S1-seq (RPM)", xlab = NA, 
     type = "l", lwd = 1.5, col = my_colors[1])
title(xlab = "Distance from DSB (nt)", line = 1.75)

y_shift <- 5
y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) + 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[1], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.5, -0.3), col = my_colors[1])

x <- agg$distance_to_DSB[!agg$transcribed]
y <- moving_average(x = agg$score[!agg$transcribed], k = k, keep = keep)
points(x = x, y = y, type = "l", lwd = 1.5, col = my_colors[2])

y <- y[x > keep]
x <- x[x > keep]
i <- find_x_where_half_AUC(x = x, y = y)
yy <- par("usr")[4] - 0.25 * (par("usr")[4] - par("usr")[3]) - 0.5 * y_shift
lines(x = c(x[i], x[i]), y = c(0, yy), col = my_colors[2], lty = "dashed")
text(x = x[i], y = yy, labels = x[i], adj = c(0.05, -0.3), col = my_colors[2])

# legend(x = "topright", legend = c("Transcribed", "Untranscribed"), col = my_colors, lwd = 1.5, bty = "n")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/fun30_4h.pdf"))
