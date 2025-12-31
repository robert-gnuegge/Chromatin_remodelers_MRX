# info --------------------------------------------------------------------
# purpose: align MNase-seq and S1-seq data at DSB-proximal nucleosome, average and plot
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 12/31/25
# last modified: 12/31/25


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
  assign(x = paste0("MNase_seq_", t), value = tmp)
}


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
  assign(x = paste0("MNase_seq_", t), value = tmp)
}

tmp <- MNase_seq_4[MNase_seq_4$dist_to_nuc <= 1100, ]
plot(x = tmp$dist_to_nuc, y = tmp$score, type = "l")

fft <- filterFFT(data = tmp$score, pcKeepComp = 0.009)  # filter noise
idx <- peakDetection(fft, threshold="10%", score=FALSE, min.cov = 0.1)

tmp$dist_to_nuc[idx]
abline(v = tmp$dist_to_nuc[idx])
