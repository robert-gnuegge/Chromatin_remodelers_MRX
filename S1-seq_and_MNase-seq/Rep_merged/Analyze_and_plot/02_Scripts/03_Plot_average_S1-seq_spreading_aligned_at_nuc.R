# info --------------------------------------------------------------------
# purpose: align MNase-seq and S1-seq data at DSB-proximal nucleosome, average and plot
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/27/24
# last modified: 11/16/24


# load libraries ----------------------------------------------------------
library(GenomicRanges)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")


# function definitions ----------------------------------------------------

# add distance to specified nucleosome
# arguments: GRanges and numeric
# result: GRanges with added mcols "dist_to_nuc" and "chosen_nuc" = "nuc"
add_distance_to_nuc <- function(seq_data, nuc_pos, nuc = 1, nuc_dist = 165){
  
  # determine if data are on fw or rev strand (needed for distance calculation etc.)
  plus_minus <- ifelse(test = as.character(unique(strand(nuc_pos))) == "+", yes = 1, no = -1)
  
  # sometimes the same idx is assigned to several nucs
  # in this case merge nucs and assign mean coordinate
  if(any(duplicated(nuc_pos$idx))){
    duplicates <- unique(nuc_pos$idx[duplicated(nuc_pos$idx)])
    for(d in duplicates){
      tmp <- nuc_pos[nuc_pos$idx == d]
      coord <- round(mean(start(tmp)))
      tmp <- tmp[1]
      ranges(tmp) <- IRanges(start = coord, width = 1)
      nuc_pos <- c(tmp, nuc_pos[nuc_pos$idx != d])
    }
    nuc_pos <- sort(nuc_pos)
  }
  
  # if "(-)nuc" cannot be found, add it placing it at the expected distance
  # from the next available nuc
  if(!(nuc %in% abs(nuc_pos$idx))){
    picked_nuc <- nuc + 1
    while(!(picked_nuc %in% abs(nuc_pos$idx))){
      picked_nuc <- picked_nuc + 1
      if(picked_nuc - nuc > 3){
        stop("No nucleosome in reasonable vicinity could be identified.")
      }
    }
    tmp <- nuc_pos[1]
    tmp$idx <- plus_minus * nuc
    coord <- start(nuc_pos[abs(nuc_pos$idx) == picked_nuc]) - plus_minus * (picked_nuc - 1) * nuc_dist
    ranges(tmp) <- IRanges(start = coord, width = 1)
    nuc_pos <- sort(c(tmp, nuc_pos))
  }
  
  seq_data$dist_to_nuc <- plus_minus * (start(seq_data) - start(nuc_pos[abs(nuc_pos$idx) == nuc]))
  seq_data$chosen_nuc <- plus_minus * nuc
  return(seq_data)
  
}

# dirs --------------------------------------------------------------------
plot_dir <- "04_Plots/MRE11/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)


# LSY4518-13B =============================================================

# load and process data ---------------------------------------------------
DSBs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
roi <- DSB_regions(DSBs = DSBs, region_width = 4000, up_rev_down_fw = TRUE)

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h_S1_seq <- subsetByIntersect(subject = LSY4518_13B_0h_S1_seq, query = roi)
LSY4518_13B_1h_S1_seq <- subsetByIntersect(subject = LSY4518_13B_1h_S1_seq, query = roi)
LSY4518_13B_2h_S1_seq <- subsetByIntersect(subject = LSY4518_13B_2h_S1_seq, query = roi)
LSY4518_13B_4h_S1_seq <- subsetByIntersect(subject = LSY4518_13B_4h_S1_seq, query = roi)

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY4518-13B_MNase-seq.RData")
LSY4518_13B_0h_MNase_seq <- subsetByIntersect(subject = LSY4518_13B_0h_MNase_seq, query = roi)
LSY4518_13B_1h_MNase_seq <- subsetByIntersect(subject = LSY4518_13B_1h_MNase_seq, query = roi)
LSY4518_13B_2h_MNase_seq <- subsetByIntersect(subject = LSY4518_13B_2h_MNase_seq, query = roi)
LSY4518_13B_4h_MNase_seq <- subsetByIntersect(subject = LSY4518_13B_4h_MNase_seq, query = roi)

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY4518-13B_nucleosome_positions.RData")
LSY4518_13B_0h_nucleosome_positions <- subsetByIntersect(subject = LSY4518_13B_0h_nucleosome_positions, query = roi)
LSY4518_13B_1h_nucleosome_positions <- subsetByIntersect(subject = LSY4518_13B_1h_nucleosome_positions, query = roi)
LSY4518_13B_2h_nucleosome_positions <- subsetByIntersect(subject = LSY4518_13B_2h_nucleosome_positions, query = roi)
LSY4518_13B_4h_nucleosome_positions <- subsetByIntersect(subject = LSY4518_13B_4h_nucleosome_positions, query = roi)


# add distance to DSB-proximal nucleosome ---------------------------------

# helper function (using the same nuc_pos object [t=0] for all seq_data objects)
add_distance_to_nuc_helper <- function(object_name){
  tmp <- GRanges()
  for(r in 1:length(roi)){
    seq_data <- subsetByIntersect(subject = get(object_name), query = roi[r])
    nuc_pos <- subsetByIntersect(subject = LSY4518_13B_0h_nucleosome_positions, query = roi[r])
    tmp <- c(tmp, add_distance_to_nuc(seq_data = seq_data, nuc_pos = nuc_pos, nuc = 1))
  }
  assign(x = object_name, value = tmp, envir = .GlobalEnv)
}

add_distance_to_nuc_helper(object_name = "LSY4518_13B_0h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_1h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_2h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_4h_S1_seq")

add_distance_to_nuc_helper(object_name = "LSY4518_13B_0h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_1h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_2h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_4h_MNase_seq")


# plotting ----------------------------------------------------------------
S1_seq_1 <- LSY4518_13B_0h_S1_seq
MNase_seq_0 <- LSY4518_13B_0h_MNase_seq

plotting_function <- function(S1_seq_1, S1_seq_2, S1_seq_4, 
                              MNase_seq_0, MNase_seq_1, MNase_seq_2, MNase_seq_4,
                              k = 51, xlim = NULL, ylim = NULL, 
                              MNase_seq_cols = gray(level = 9:6 * 0.1), 
                              S1_seq_cols = JFly_colors[c(1, 4, 5)],
                              nuc_marks = 1:6, nuc_mark_col = gray(level = 0.4)){
  
  MNase_agg_0 <- aggregate(score ~ dist_to_nuc, data = MNase_seq_0, FUN = mean)
  MNase_agg_0$score <- runmed(x = MNase_agg_0$score, k = k)
  
  plot_polygon <- function(MNase_seq, col){
    polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]), 
            y = c(0, y, 0), col = col, border = NA)
  }
  
  
  S1_agg_1 <- aggregate(score ~ dist_to_nuc, data = S1_seq_1, FUN = mean)
  S1_agg_2 <- aggregate(score ~ dist_to_nuc, data = S1_seq_2, FUN = mean)
  S1_agg_4 <- aggregate(score ~ dist_to_nuc, data = S1_seq_4, FUN = mean)
  
  
  
}


pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

# start empty plot
plot(x = NA, y = NA, xlim = c(-80, 1100), ylim = c(0, 11), 
     axes = FALSE, ann = FALSE)

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

# add nucleosome marks
x <- 1:6
nuc <- 1
segments(x0 = (x - nuc) * nuc_dist, y0 = 11, x1 = (x - nuc) * nuc_dist, y1 = 0, col = gray(level = 0.4), lty = "dashed")

# add S1-seq plots
par(new = TRUE)

S1_seq <- S1_agg_1
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
plot(x = S1_seq$dist_from_nuc + nuc_dist, y = y, ylim = c(0, 15), xlim = c(-80, 1100), type = "l", col = JFly_colors[1], 
     axes = FALSE, ann = FALSE)

S1_seq <- S1_agg_2
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 2 * nuc_dist, y = y, type = "l", col = JFly_colors[4])

S1_seq <- S1_agg_4
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 3 * nuc_dist, y = y, type = "l", col = JFly_colors[5])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4518-13B.pdf"))


# LSY5415 =============================================================

# load and process S1-seq data --------------------------------------------
DSBs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)  # keep only regions with correct orientation w.r.t DSBs

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h_S1_seq <- process_S1_seq(GRanges = LSY5415_0h_S1_seq, roi = roi)
LSY5415_1h_S1_seq <- process_S1_seq(GRanges = LSY5415_1h_S1_seq, roi = roi)
LSY5415_2h_S1_seq <- process_S1_seq(GRanges = LSY5415_2h_S1_seq, roi = roi)
LSY5415_4h_S1_seq <- process_S1_seq(GRanges = LSY5415_4h_S1_seq, roi = roi)


# load and process MNase-seq data -----------------------------------------
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq.RData")
LSY5415_0h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_0h_MNase_seq, roi = roi)
LSY5415_1h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_1h_MNase_seq, roi = roi)
LSY5415_2h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_2h_MNase_seq, roi = roi)
LSY5415_4h_MNase_seq <- process_MNase_seq(GRanges = LSY5415_4h_MNase_seq, roi = roi)


# process nucleosome center data ------------------------------------------

# load nucleosome position data 
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5415_nucleosome_positions.RData")

# nucleosome properties (according to Jansen et al., 2011; pmid: 21646431)
nuc_width <- 147
nuc_dist <- 165  # nuc_width + 18 (average linker length)

# define ideal nucleosome positions relative to DSBs
ideal_nuc_centers <- GRanges()
for(n in 0:19){
  tmp <- shift(x = DSBs, shift = -round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- -(n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
  tmp <- shift(x = DSBs, shift = round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- (n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
}
ideal_nuc_centers <- sort(ideal_nuc_centers)
strand(ideal_nuc_centers) <- ifelse(test = ideal_nuc_centers$nuc_number < 0, yes = "-", no = "+")

# add nucleosome numbers and strand information to Nuc_centers GRanges  
LSY5415_0h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5415_0h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5415_1h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5415_1h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5415_2h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5415_2h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5415_4h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5415_4h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)

# calc distance from DSB-proximal nucleosome for sequencing data sets ------
DSBs <- SrfIcs[-c(9, 17)]
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)

# t = 0
nuc <- 1
MNase <- add_distance_to_nucleosome(GRanges = LSY5415_0h_MNase_seq, nuc_centers = LSY5415_0h_nucleosome_positions, roi = roi, nuc = nuc)
S1 <- add_distance_to_nucleosome(GRanges = LSY5415_0h_S1_seq, nuc_centers = LSY5415_0h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_0 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_0 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 1
nuc <- 2
MNase <- add_distance_to_nucleosome(GRanges = LSY5415_1h_MNase_seq, nuc_centers = LSY5415_1h_nucleosome_positions, roi = roi, nuc = nuc)
S1 <- add_distance_to_nucleosome(GRanges = LSY5415_1h_S1_seq, nuc_centers = LSY5415_1h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_1 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_1 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 2
nuc <- 3
MNase <- add_distance_to_nucleosome(GRanges = LSY5415_2h_MNase_seq, nuc_centers = LSY5415_2h_nucleosome_positions, roi = roi, nuc = nuc)
S1 <- add_distance_to_nucleosome(GRanges = LSY5415_2h_S1_seq, nuc_centers = LSY5415_2h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_2 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_2 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 4
nuc <- 4
MNase <- add_distance_to_nucleosome(GRanges = LSY5415_4h_MNase_seq, nuc_centers = LSY5415_4h_nucleosome_positions, roi = roi, nuc = nuc)
S1 <- add_distance_to_nucleosome(GRanges = LSY5415_4h_S1_seq, nuc_centers = LSY5415_4h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- S1$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
S1_agg_4 <- aggregate(score ~ dist_from_nuc, data = S1[idx], FUN = mean)
MNase_agg_4 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)


# plotting ----------------------------------------------------------------
k <- 51

pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

# start empty plot
plot(x = NA, y = NA, xlim = c(-80, 1100), ylim = c(-11.5, 0), 
     axes = FALSE, ann = FALSE)

# plot MNase-seq
MNase_seq <- MNase_agg_0
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]), 
        y = c(0, -y, 0), col = gray(level = 0.9), border = NA)

MNase_seq <- MNase_agg_1
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + nuc_dist,
        y = c(0, -y, 0), col = gray(level = 0.8), border = NA)

MNase_seq <- MNase_agg_2
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 2 * nuc_dist,
        y = c(0, -y, 0), col = gray(level = 0.7), border = NA)

MNase_seq <- MNase_agg_4
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 3 * nuc_dist,
        y = c(0, -y, 0), col = gray(level = 0.6), border = NA)

# add nucleosome marks
x <- 1:6
nuc <- 1
segments(x0 = (x - nuc) * nuc_dist, y0 = 0, x1 = (x - nuc) * nuc_dist, y1 = -11.5, col = gray(level = 0.4), lty = "dashed")

# add S1-seq plots
par(new = TRUE)

S1_seq <- S1_agg_1
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
plot(x = S1_seq$dist_from_nuc + nuc_dist, y = -y, ylim = c(-28,0), xlim = c(-80, 1100), type = "l", col = JFly_colors[1],
     axes = FALSE, ann = FALSE)

S1_seq <- S1_agg_2
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 2 * nuc_dist, y = -y, type = "l", col = JFly_colors[4])

S1_seq <- S1_agg_4
S1_seq <- S1_seq[S1_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = S1_seq$score, k = k, keep = 0)
points(x = S1_seq$dist_from_nuc + 3 * nuc_dist, y = -y, type = "l", col = JFly_colors[5])

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5415.pdf"))


# LSY5935 =============================================================

# load and process MNase-seq data -----------------------------------------
DSBs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)  # keep only regions with correct orientation w.r.t DSBs

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5935_MNase-seq.RData")
LSY5935_0h_MNase_seq <- process_MNase_seq(GRanges = LSY5935_0h_MNase_seq, roi = roi)
LSY5935_1h_MNase_seq <- process_MNase_seq(GRanges = LSY5935_1h_MNase_seq, roi = roi)
LSY5935_2h_MNase_seq <- process_MNase_seq(GRanges = LSY5935_2h_MNase_seq, roi = roi)
LSY5935_4h_MNase_seq <- process_MNase_seq(GRanges = LSY5935_4h_MNase_seq, roi = roi)


# process nucleosome center data ------------------------------------------

# load nucleosome position data 
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5935_nucleosome_positions.RData")

# nucleosome properties (according to Jansen et al., 2011; pmid: 21646431)
nuc_width <- 147
nuc_dist <- 165  # nuc_width + 18 (average linker length)

# define ideal nucleosome positions relative to DSBs
ideal_nuc_centers <- GRanges()
for(n in 0:19){
  tmp <- shift(x = DSBs, shift = -round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- -(n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
  tmp <- shift(x = DSBs, shift = round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- (n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
}
ideal_nuc_centers <- sort(ideal_nuc_centers)
strand(ideal_nuc_centers) <- ifelse(test = ideal_nuc_centers$nuc_number < 0, yes = "-", no = "+")

# add nucleosome numbers and strand information to Nuc_centers GRanges  
LSY5935_0h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5935_0h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5935_1h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5935_1h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5935_2h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5935_2h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5935_4h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5935_4h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)

# calc distance from DSB-proximal nucleosome for sequencing data sets ------
DSBs <- SrfIcs[-c(9, 17)]
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)

# t = 0
nuc <- 1
MNase <- add_distance_to_nucleosome(GRanges = LSY5935_0h_MNase_seq, nuc_centers = LSY5935_0h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_0 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 1
nuc <- 2
MNase <- add_distance_to_nucleosome(GRanges = LSY5935_1h_MNase_seq, nuc_centers = LSY5935_1h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_1 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 2
nuc <- 3
MNase <- add_distance_to_nucleosome(GRanges = LSY5935_2h_MNase_seq, nuc_centers = LSY5935_2h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_2 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 4
nuc <- 4
MNase <- add_distance_to_nucleosome(GRanges = LSY5935_4h_MNase_seq, nuc_centers = LSY5935_4h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_4 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)


# plotting ----------------------------------------------------------------
k <- 51

# dirs --------------------------------------------------------------------
plot_dir <- "04_Plots/mre11-nd/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)

pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

# start empty plot
plot(x = NA, y = NA, xlim = c(-80, 1100), ylim = c(0, 11), 
     axes = FALSE, ann = FALSE)

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

# add nucleosome marks
x <- 1:6
nuc <- 1
segments(x0 = (x - nuc) * nuc_dist, y0 = 11, x1 = (x - nuc) * nuc_dist, y1 = 0, col = gray(level = 0.4), lty = "dashed")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5935.pdf"))


# LSY5934 =============================================================

# load and process MNase-seq data -----------------------------------------
DSBs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)  # keep only regions with correct orientation w.r.t DSBs

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5934_MNase-seq.RData")
LSY5934_0h_MNase_seq <- process_MNase_seq(GRanges = LSY5934_0h_MNase_seq, roi = roi)
LSY5934_1h_MNase_seq <- process_MNase_seq(GRanges = LSY5934_1h_MNase_seq, roi = roi)
LSY5934_2h_MNase_seq <- process_MNase_seq(GRanges = LSY5934_2h_MNase_seq, roi = roi)
LSY5934_4h_MNase_seq <- process_MNase_seq(GRanges = LSY5934_4h_MNase_seq, roi = roi)


# process nucleosome center data ------------------------------------------

# load nucleosome position data 
load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5934_nucleosome_positions.RData")

# nucleosome properties (according to Jansen et al., 2011; pmid: 21646431)
nuc_width <- 147
nuc_dist <- 165  # nuc_width + 18 (average linker length)

# define ideal nucleosome positions relative to DSBs
ideal_nuc_centers <- GRanges()
for(n in 0:19){
  tmp <- shift(x = DSBs, shift = -round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- -(n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
  tmp <- shift(x = DSBs, shift = round((n + 0.5) * nuc_dist))
  tmp$nuc_number <- (n + 1)
  ideal_nuc_centers <- c(ideal_nuc_centers, tmp)
}
ideal_nuc_centers <- sort(ideal_nuc_centers)
strand(ideal_nuc_centers) <- ifelse(test = ideal_nuc_centers$nuc_number < 0, yes = "-", no = "+")

# add nucleosome numbers and strand information to Nuc_centers GRanges  
LSY5934_0h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5934_0h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5934_1h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5934_1h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5934_2h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5934_2h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)
LSY5934_4h_nucleosome_positions <- add_nuc_numbers_and_strand(DSBs = DSBs, nuc_centers = LSY5934_4h_nucleosome_positions, ideal_nuc_centers = ideal_nuc_centers)

# calc distance from DSB-proximal nucleosome for sequencing data sets ------
DSBs <- SrfIcs[-c(9, 17)]
roi <- DSB_regions(DSBs = DSBs, region_width = 6000, up_rev_down_fw = TRUE)

# t = 0
nuc <- 1
MNase <- add_distance_to_nucleosome(GRanges = LSY5934_0h_MNase_seq, nuc_centers = LSY5934_0h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_0 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 1
nuc <- 2
MNase <- add_distance_to_nucleosome(GRanges = LSY5934_1h_MNase_seq, nuc_centers = LSY5934_1h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_1 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 2
nuc <- 3
MNase <- add_distance_to_nucleosome(GRanges = LSY5934_2h_MNase_seq, nuc_centers = LSY5934_2h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_2 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)

# t = 4
nuc <- 4
MNase <- add_distance_to_nucleosome(GRanges = LSY5934_4h_MNase_seq, nuc_centers = LSY5934_4h_nucleosome_positions, roi = roi, nuc = nuc)
idx <- MNase$dist_from_nuc > -(nuc - 0.5) * nuc_dist  # & S1$dist_from_nuc <= 1200
MNase_agg_4 <- aggregate(score ~ dist_from_nuc, data = MNase[idx], FUN = mean)


# plotting ----------------------------------------------------------------
k <- 51

# dirs --------------------------------------------------------------------
plot_dir <- "04_Plots/mre11-nd/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)

pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

# start empty plot
plot(x = NA, y = NA, xlim = c(-80, 1100), ylim = c(-11,0), 
     axes = FALSE, ann = FALSE)

# plot MNase-seq
MNase_seq <- MNase_agg_0
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]), 
        y = c(0, -y, 0), col = gray(level = 0.9), border = NA)

MNase_seq <- MNase_agg_1
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + nuc_dist, 
        y = c(0, -y, 0), col = gray(level = 0.8), border = NA)

MNase_seq <- MNase_agg_2
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 2 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 2 * nuc_dist, 
        y = c(0, -y, 0), col = gray(level = 0.7), border = NA)

MNase_seq <- MNase_agg_4
MNase_seq <- MNase_seq[MNase_seq$dist_from_nuc < 1100 - 3 * nuc_dist, ]
y <- moving_average(x = MNase_seq$score, k = k, keep = 0)
polygon(x = c(MNase_seq$dist_from_nuc[1], MNase_seq$dist_from_nuc, MNase_seq$dist_from_nuc[length(MNase_seq$dist_from_nuc)]) + 3 * nuc_dist, 
        y = c(0, -y, 0), col = gray(level = 0.6), border = NA)

# add nucleosome marks
x <- 1:6
nuc <- 1
segments(x0 = (x - nuc) * nuc_dist, y0 = -11, x1 = (x - nuc) * nuc_dist, y1 = 0, col = gray(level = 0.4), lty = "dashed")

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5934.pdf"))
