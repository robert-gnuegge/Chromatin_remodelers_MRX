# info --------------------------------------------------------------------
# purpose: align MNase-seq and S1-seq data at DSB-proximal nucleosome, average and plot
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/26/24
# last modified: 05/27/24

# load libraries ----------------------------------------------------------
library(GenomicRanges)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

# function definitions ----------------------------------------------------

# find DSB proximal nucleosome
# arguments: GRanges objects and logical
# result: GRanges object
find_DSB_prox_nuc <- function(DSB, nuc_pos_GRanges, upstream = TRUE){
  roi <- flank(x = DSB, width = 3000, start = upstream) 
  tmp <- subsetByIntersect(subject = nuc_pos_GRanges, query = roi)
  idx <- nearest(x = DSB, subject = tmp)
  prox_nuc <- range(tmp[idx])
  mcols(prox_nuc) <- data.frame(distance = distance(x = DSB, y = prox_nuc), upstream = upstream)
  return(prox_nuc)
}


# load and process S1-seq data --------------------------------------------
SrfIcs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
roi <- DSB_regions(DSBs = SrfIcs, region_width = 5000)
# restrict data processing to regions around DSBs
# chose a bigger region size than used for plotting to allow smoothing

process_S1_seq <- function(GRanges, roi){
  tmp <- subsetByIntersect(subject = GRanges, query = roi)  # only keep S1-seq coverage in DSB regions
  tmp <- sort(as_nt_resolved_GRanges(tmp), ignore.strand = TRUE)  # nt resolution, and sort
  return(tmp)
}

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_0h_S1_seq, roi = roi)
LSY4518_13B_1h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_1h_S1_seq, roi = roi)
LSY4518_13B_2h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_2h_S1_seq, roi = roi)
LSY4518_13B_4h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_4h_S1_seq, roi = roi)

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h_S1_seq <- process_S1_seq(GRanges = LSY5415_0h_S1_seq, roi = roi)
LSY5415_1h_S1_seq <- process_S1_seq(GRanges = LSY5415_1h_S1_seq, roi = roi)
LSY5415_2h_S1_seq <- process_S1_seq(GRanges = LSY5415_2h_S1_seq, roi = roi)
LSY5415_4h_S1_seq <- process_S1_seq(GRanges = LSY5415_4h_S1_seq, roi = roi)


# load and process MNase-seq data -----------------------------------------
process_MNase_seq <- function(GRanges, roi){
  tmp <- subsetByIntersect(subject = GRanges, query = roi)  # only keep MNase-seq coverage in DSB regions
  tmp <- sort(as_nt_resolved_GRanges(tmp), ignore.strand = TRUE)  # nt resolution, and sort
  return(tmp)
}

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY4518-13B_MNase-seq.RData")
LSY4518_13B_0h_MNase_seq <- process_MNase_seq(GRanges = LSY4518_13B_0h_MNase_seq, roi = roi)
LSY4518_13B_1h_MNase_seq <- process_MNase_seq(GRanges = LSY4518_13B_1h_MNase_seq, roi = roi)
LSY4518_13B_2h_MNase_seq <- process_MNase_seq(GRanges = LSY4518_13B_2h_MNase_seq, roi = roi)
LSY4518_13B_4h_MNase_seq <- process_MNase_seq(GRanges = LSY4518_13B_4h_MNase_seq, roi = roi)

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq.RData")
LSY5415_0h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_0h_MNase_seq, roi = roi)
LSY5415_1h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_1h_MNase_seq, roi = roi)
LSY5415_2h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_2h_MNase_seq, roi = roi)
LSY5415_4h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_4h_MNase_seq, roi = roi)



# find DSB-proximal nucleosomes -------------------------------------------

# load nucleosome position data 
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY4518-13B_nucleosome_positions.RData")
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5415_nucleosome_positions.RData")
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5934_nucleosome_positions.RData")
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5935_nucleosome_positions.RData")

# function to find DSB proximal DSBs
find_DSB_prox_nuc <- function(DSB, nuc_pos_GRanges){
  
  roi <- flank(x = DSB, width = 3000, start = TRUE) 
  tmp <- subsetByIntersect(subject = nuc_pos_GRanges, query = roi)
  idx <- nearest(x = DSB, subject = tmp)
  up_prox_nuc <- granges(tmp[idx])
  mcols(up_prox_nuc) <- data.frame(distance = distance(x = DSB, y = up_prox_nuc), upstream = TRUE)
  
  roi <- flank(x = DSB, width = 3000, start = FALSE) 
  tmp <- subsetByIntersect(subject = nuc_pos_GRanges, query = roi)
  idx <- nearest(x = DSB, subject = tmp)
  dn_prox_nuc <- granges(tmp[idx])
  mcols(dn_prox_nuc) <- data.frame(distance = distance(x = DSB, y = dn_prox_nuc), upstream = FALSE)
  
  prox_nuc <- sort(x = c(up_prox_nuc, dn_prox_nuc))
  
  return(prox_nuc)

}

# find DSB-proximal nucleosomes
LSY4518_13B_0h_DSB_prox_nuc <- find_DSB_prox_nuc(DSB = SrfIcs, nuc_pos_GRanges = LSY4518_13B_0h_nucleosome_positions)
LSY5415_0h_DSB_prox_nuc <- find_DSB_prox_nuc(DSB = SrfIcs, nuc_pos_GRanges = LSY5415_0h_nucleosome_positions)
LSY5934_0h_DSB_prox_nuc <- find_DSB_prox_nuc(DSB = SrfIcs, nuc_pos_GRanges = LSY5934_0h_nucleosome_positions)
LSY5935_0h_DSB_prox_nuc <- find_DSB_prox_nuc(DSB = SrfIcs, nuc_pos_GRanges = LSY5935_0h_nucleosome_positions)

# add distance to DSB-proximal nucleosome to coverage data
add_distance_to_nucleosome <- function(GRanges, DSBs, DSB_prox_nucs, nuc_dist = 165){
  
  roi <- flank(x = DSBs, width = 3000, start = TRUE) 
  tmp <- subsetByIntersect(subject = GRanges, query = roi)
  distanceToNearest(x = )
  
  idx <- nearest(x = DSB, subject = tmp)
  up_prox_nuc <- granges(tmp[idx])
  mcols(up_prox_nuc) <- data.frame(distance = distance(x = DSB, y = up_prox_nuc), upstream = TRUE)
  
  
  out <- GRanges()
  for(n in 1:length(DSB_prox_nucs)){
    nuc_region <- subsetByOverlaps(x = nuc_centers, ranges = roi[n])
    nuc_pos <- nuc_region[which.min(abs(nuc_region$nuc_number) - nuc)]
    tmp <- subsetByOverlaps(x = GRanges, ranges = roi[n])
    if(unique(strand(tmp)) == "-"){
      tmp$dist_from_nuc <- start(nuc_pos) - start(tmp) + (abs(nuc_pos$nuc_number) - nuc) * nuc_dist
    }else{
      tmp$dist_from_nuc <- start(tmp) - start(nuc_pos) + (nuc_pos$nuc_number - nuc) * nuc_dist
    }
    tmp$nuc_number <- nuc_pos$nuc_number 
    out <- c(out, tmp)
  }
  return(out)
}


roi <- flank(x = SrfIcs, width = 3000, start = TRUE) 
tmp <- subsetByIntersect(subject = LSY4518_13B_0h_S1_seq, query = roi)
DSB_prox_nuc <- LSY4518_13B_0h_DSB_prox_nuc[LSY4518_13B_0h_DSB_prox_nuc$upstream == TRUE]

tmp <- tmp[strand(tmp) == "-"]
dist_to_nearest <- distanceToNearest(x = tmp, subject = )

# adjust for nucleosomes != -1
adjustment_value <- (ceiling(mcols(LSY4518_13B_0h_DSB_prox_nuc[LSY4518_13B_0h_DSB_prox_nuc$upstream == TRUE])$distance / nuc_dist) - 1) * nuc_dist
rep(adjustment_value, rle(subjectHits(dist_to_nearest))$lengths)


idx <- nearest(x = DSB, subject = tmp)
up_prox_nuc <- granges(tmp[idx])
mcols(up_prox_nuc) <- data.frame(distance = distance(x = DSB, y = up_prox_nuc), upstream = TRUE)



# plotting ================================================================


# all t -------------------------------------------------------------------
k <- 51

pdf(file = "tmp.pdf", width=3.75, height=3)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 1, 2.9, -1), las = 1, tcl = -0.25, mgp = c(1.75, 0.4, 0))

# start empty plot
plot(x = NA, y = NA, xlim = c(-80, 1100), ylim = c(0,30.65), xaxt = "n",
     xlab = "Distance from +1 dyade [nt]", ylab = NA, yaxt = "n")
axis(side = 1, at = 0:4 * 250)
axis_col <- gray(level = 0.6)
axis(side = 4, at = pretty(c(0,30.65)), col = axis_col, col.axis = axis_col)
text(x = par("usr")[2]*1.2, y = 0.9* mean(par("usr")[3:4]), labels = "Average MNase-seq [RPM]", srt = -90, pos = 3, col = axis_col, xpd = TRUE)

# plot MNase-seq
MNase_seq <- MNase_agg_0
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]), 
        y = c(0, y, 0), col = gray(level = 0.9), border = NA)

MNase_seq <- MNase_agg_1
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + nuc_dist, 
        y = c(0, y, 0), col = gray(level = 0.8), border = NA)

MNase_seq <- MNase_agg_2
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 2 * nuc_dist, 
        y = c(0, y, 0), col = gray(level = 0.7), border = NA)

MNase_seq <- MNase_agg_4
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 3 * nuc_dist, 
        y = c(0, y, 0), col = gray(level = 0.6), border = NA)

# add nucleosome labels 
x <- 1:6
nuc <- 1
segments(x0 = (x - nuc) * nuc_dist, y0 = 1.01 * par("usr")[4], x1 = (x - nuc) * nuc_dist, y1 = 0, col = gray(level = 0.4), lty = "dashed", xpd = TRUE)
text(x = (x - nuc) * nuc_dist, y = par("usr")[4], labels = paste0("+", x), pos = 3, xpd = TRUE, col = gray(level = 0.4))

# add S1-seq plots
par(new = TRUE)

S1_seq <- S1_agg_1
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
plot(x = S1_seq$dist_from_nuc + nuc_dist, y = y, ylim = c(0, 26), xlim = c(-80, 1100), type = "l", col = JFly_colors[1], axes = FALSE, ann = FALSE)

S1_seq <- S1_agg_2
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 2 * nuc_dist, y = y, type = "l", col = JFly_colors[4])

S1_seq <- S1_agg_4
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 3 * nuc_dist, y = y, type = "l", col = JFly_colors[5])

axis(side = 2, at = pretty(c(0,26)))
mtext("Average S1-seq [RPM]", side = 2, line = 2, las = 0)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Aligned_S1_seq_and_MNase-seq.pdf"))
