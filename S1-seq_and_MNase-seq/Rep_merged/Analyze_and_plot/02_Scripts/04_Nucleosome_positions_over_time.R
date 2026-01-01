# info --------------------------------------------------------------------
# purpose: align MNase-seq and S1-seq data at DSB-proximal nucleosome, average and plot
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 12/31/25
# last modified: 01/01/26


# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(nucleR)

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
  
  # if "nuc" cannot be found, add it placing it at the expected distance
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

# define DSBs to be analyzed ==============================================
DSBs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
DSBs <- DSBs[DSBs$DSB_kinetics_rank < 17] # exclude very slowly formed DSBs
roi <- DSB_regions(DSBs = DSBs, region_width = 4000, up_rev_down_fw = TRUE)

# LSY4518-13B =============================================================

# load and process data ---------------------------------------------------

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
  # add also distance to nearest DSB
  # immediately adjacent positions are also resulting in a 0 distance, therefore 1 is added to all distances and then the SrfIcs positions are set to 0
  tmp$distance_to_DSB <- mcols(distanceToNearest(x = tmp, subject = DSBs))$distance + 1
  mcols(tmp[nearest(x = DSBs, subject = tmp) - 1])$distance_to_DSB <- 0
  assign(x = object_name, value = tmp, envir = .GlobalEnv)
}

add_distance_to_nuc_helper(object_name = "LSY4518_13B_0h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_1h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_2h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_4h_MNase_seq")


# average and smoothen ----------------------------------------------------
k <- 51

# MNase-seq data
for(t in c(0, 1, 2, 4)){
  tmp <- aggregate(cbind(score, distance_to_DSB) ~ dist_to_nuc, 
                   data = get(paste0("LSY4518_13B_", t, "h_MNase_seq")), FUN = mean)
  tmp$score <- runmed(x = tmp$score, k = k)
  assign(x = paste0("WT_MNase_seq_", t), value = tmp)
}


# find t = 0 nucleosome positions -----------------------------------------
tmp <- WT_MNase_seq_0[WT_MNase_seq_0$dist_to_nuc >= -100 & WT_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.009)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

WT <- tmp[idx, ]
WT$nuc_ID <- 1:8
WT$time <- 0

ideal_nucs <- WT[WT$time == 0, c("dist_to_nuc", "nuc_ID")]

# find t = 1 nucleosome positions -----------------------------------------
tmp <- WT_MNase_seq_1[WT_MNase_seq_1$dist_to_nuc >= -100 & WT_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.01)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

tmp <- tmp[idx, ]

for(n in (1:nrow(tmp))){
  idx <- which.min(abs(tmp$dist_to_nuc[n] - ideal_nucs$dist_to_nuc))
  tmp$nuc_ID[n] <- ideal_nucs$nuc_ID[idx]
}

tmp$time <- 1

WT <- rbind(WT, tmp)


# find t = 2 nucleosome positions -----------------------------------------
tmp <- WT_MNase_seq_2[WT_MNase_seq_2$dist_to_nuc >= -100 & WT_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.01)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

tmp <- tmp[idx, ]

for(n in (1:nrow(tmp))){
  idx <- which.min(abs(tmp$dist_to_nuc[n] - ideal_nucs$dist_to_nuc))
  tmp$nuc_ID[n] <- ideal_nucs$nuc_ID[idx]
}

tmp$time <- 2

WT <- rbind(WT, tmp)


# find t = 4 nucleosome positions -----------------------------------------
tmp <- WT_MNase_seq_4[WT_MNase_seq_4$dist_to_nuc >= -100 & WT_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.01)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

tmp <- tmp[idx, ]

for(n in (1:nrow(tmp))){
  idx <- which.min(abs(tmp$dist_to_nuc[n] - ideal_nucs$dist_to_nuc))
  tmp$nuc_ID[n] <- ideal_nucs$nuc_ID[idx]
}

tmp$time <- 4

WT <- rbind(WT, tmp)
WT <- WT[WT$nuc_ID <= 7, ]



# plotting ----------------------------------------------------------------
library(shape)

plot_dir <- "04_Plots/MRE11/DSB_prxml_nucs_over_time"
dir.create(path = plot_dir, showWarnings = FALSE)

WT <- WT[!(WT$nuc_ID == 1 & WT$time > 0), ]

pdf(file = "tmp.pdf", width=3, height=1.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

y <- c(4:1) * 1
plot(x = NA, y = NA, xlim = c(-200, 1080), ylim = c(0.75 * min(y), 1.25 * max(y)), axes = FALSE, ann = FALSE)

x <- WT$dist_to_nuc[WT$time == 0 & WT$nuc_ID <= 7]
text(x = x, y = (max(y) + 0.33 * min(y)), labels = paste0("+", 1:7), pos = 3)
segments(x0 = x, x1 = x, y0 = min(y), y1 = max(y), lty = "dashed")

text(x = -80, y = y, labels = paste0(c(0, 1, 2, 4), " h"), pos = 2)

mapping <- data.frame(t = c(0, 1, 2, 4), idx = 1:4)
WT$y <- y[match(x = WT$time, table = mapping$t)]

points(x = WT$dist_to_nuc, y = WT$y, pch = 21, bg = gray(level = 1 - (WT$score - min(WT$score)) / (max(WT$score) - min(WT$score))), cex = 2)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4518-13B.pdf"))



# LSY5415 =============================================================

# load and process data ---------------------------------------------------

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq.RData")
LSY5415_0h_MNase_seq <- subsetByIntersect(subject = LSY5415_0h_MNase_seq, query = roi)
LSY5415_1h_MNase_seq <- subsetByIntersect(subject = LSY5415_1h_MNase_seq, query = roi)
LSY5415_2h_MNase_seq <- subsetByIntersect(subject = LSY5415_2h_MNase_seq, query = roi)
LSY5415_4h_MNase_seq <- subsetByIntersect(subject = LSY5415_4h_MNase_seq, query = roi)

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5415_nucleosome_positions.RData")
LSY5415_0h_nucleosome_positions <- subsetByIntersect(subject = LSY5415_0h_nucleosome_positions, query = roi)
LSY5415_1h_nucleosome_positions <- subsetByIntersect(subject = LSY5415_1h_nucleosome_positions, query = roi)
LSY5415_2h_nucleosome_positions <- subsetByIntersect(subject = LSY5415_2h_nucleosome_positions, query = roi)
LSY5415_4h_nucleosome_positions <- subsetByIntersect(subject = LSY5415_4h_nucleosome_positions, query = roi)


# add distance to DSB-proximal nucleosome ---------------------------------

# helper function (using the same nuc_pos object [t=0] for all seq_data objects)
add_distance_to_nuc_helper <- function(object_name){
  tmp <- GRanges()
  for(r in 1:length(roi)){
    seq_data <- subsetByIntersect(subject = get(object_name), query = roi[r])
    nuc_pos <- subsetByIntersect(subject = LSY5415_0h_nucleosome_positions, query = roi[r])
    tmp <- c(tmp, add_distance_to_nuc(seq_data = seq_data, nuc_pos = nuc_pos, nuc = 1))
  }
  assign(x = object_name, value = tmp, envir = .GlobalEnv)
}

add_distance_to_nuc_helper(object_name = "LSY5415_0h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5415_1h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5415_2h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5415_4h_MNase_seq")


# average and smoothen ----------------------------------------------------
k <- 51

# MNase-seq data
for(t in c(0, 1, 2, 4)){
  tmp <- aggregate(score ~ dist_to_nuc, 
                   data = get(paste0("LSY5415_", t, "h_MNase_seq")), FUN = mean)
  tmp$score <- runmed(x = tmp$score, k = k)
  assign(x = paste0("fun30_MNase_seq_", t), value = tmp)
}

# find t = 0 nucleosome positions -----------------------------------------
tmp <- fun30_MNase_seq_0[fun30_MNase_seq_0$dist_to_nuc >= -100 & fun30_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.009)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

fun30 <- tmp[idx, ]
fun30$nuc_ID <- 1:8
fun30$time <- 0

ideal_nucs <- fun30[fun30$time == 0, c("dist_to_nuc", "nuc_ID")]

# find t = 1 nucleosome positions -----------------------------------------
tmp <- fun30_MNase_seq_1[fun30_MNase_seq_1$dist_to_nuc >= -100 & fun30_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.01)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

tmp <- tmp[idx, ]

for(n in (1:nrow(tmp))){
  idx <- which.min(abs(tmp$dist_to_nuc[n] - ideal_nucs$dist_to_nuc))
  tmp$nuc_ID[n] <- ideal_nucs$nuc_ID[idx]
}

tmp$time <- 1

fun30 <- rbind(fun30, tmp)


# find t = 2 nucleosome positions -----------------------------------------
tmp <- fun30_MNase_seq_2[fun30_MNase_seq_2$dist_to_nuc >= -100 & fun30_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.01)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

tmp <- tmp[idx, ]

for(n in (1:nrow(tmp))){
  idx <- which.min(abs(tmp$dist_to_nuc[n] - ideal_nucs$dist_to_nuc))
  tmp$nuc_ID[n] <- ideal_nucs$nuc_ID[idx]
}

tmp$time <- 2

fun30 <- rbind(fun30, tmp)


# find t = 4 nucleosome positions -----------------------------------------
tmp <- fun30_MNase_seq_4[fun30_MNase_seq_4$dist_to_nuc >= -100 & fun30_MNase_seq_0$dist_to_nuc <= 1200, 1:2]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.0105)  # filter noise
points(x = tmp$dist_to_nuc, y = fft, type = "l", col = "red")

idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov=0.1)
abline(v = tmp$dist_to_nuc[idx], col = "blue")

tmp <- tmp[idx, ]

for(n in (1:nrow(tmp))){
  idx <- which.min(abs(tmp$dist_to_nuc[n] - ideal_nucs$dist_to_nuc))
  tmp$nuc_ID[n] <- ideal_nucs$nuc_ID[idx]
}

tmp$time <- 4

fun30 <- rbind(fun30, tmp)
old_fun30 <- fun30
fun30 <- fun30[fun30$nuc_ID <= 7, ]



# plotting ----------------------------------------------------------------
library(shape)

plot_dir <- "04_Plots/MRE11/DSB_prxml_nucs_over_time"
dir.create(path = plot_dir, showWarnings = FALSE)

fun30 <- fun30[!(fun30$nuc_ID == 1 & fun30$time > 1), ]

pdf(file = "tmp.pdf", width=3, height=1.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

y <- c(4:1) * 1
plot(x = NA, y = NA, xlim = c(-200, 1080), ylim = c(0.75 * min(y), 1.25 * max(y)), axes = FALSE, ann = FALSE)

x <- fun30$dist_to_nuc[fun30$time == 0 & fun30$nuc_ID <= 7]
text(x = x, y = (max(y) + 0.33 * min(y)), labels = paste0("+", 1:7), pos = 3)
segments(x0 = x, x1 = x, y0 = min(y), y1 = max(y), lty = "dashed")

text(x = -80, y = y, labels = paste0(c(0, 1, 2, 4), " h"), pos = 2)

mapping <- data.frame(t = c(0, 1, 2, 4), idx = 1:4)
fun30$y <- y[match(x = fun30$time, table = mapping$t)]

points(x = fun30$dist_to_nuc, y = fun30$y, pch = 21, bg = gray(level = 1 - (fun30$score - min(fun30$score)) / (max(fun30$score) - min(fun30$score))), cex = 2)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5415.pdf"))
