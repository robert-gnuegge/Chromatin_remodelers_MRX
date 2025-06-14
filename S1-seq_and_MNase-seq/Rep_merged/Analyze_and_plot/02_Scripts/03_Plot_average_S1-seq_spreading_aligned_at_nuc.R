# info --------------------------------------------------------------------
# purpose: align MNase-seq and S1-seq data at DSB-proximal nucleosome, average and plot
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/27/24
# last modified: 06/14/25


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

# plotting function
# arguments: dataframes with columns "dist_to_nuc" and "score", numeric, colors
# result: plot
plotting_function <- function(MNase_seq_0, MNase_seq_1, MNase_seq_2, MNase_seq_4,
                              S1_seq_1 = NULL, S1_seq_2 = NULL, S1_seq_4 = NULL,
                              xlim = NULL, ylim_MNase = NULL, ylim_S1 = NULL,
                              MNase_seq_cols = gray(level = 9:6 * 0.1),
                              S1_seq_cols = JFly_colors[c(1, 4, 5)],
                              nuc_marks = NULL, nuc_dist = 165, nuc_mark_col = gray(level = 0.4)){
  
  # function to plot MNase-seq data (as polygon)
  plot_MNase_seq <- function(data, col){
    polygon(x = c(data$dist_to_nuc[1], data$dist_to_nuc, data$dist_to_nuc[length(data$dist_to_nuc)]), 
            y = c(0, data$score, 0), col = col, border = NA)
  }
  
  # function to plot S1-seq data (as line plot)
  plot_S1_seq <- function(data, col){
    points(x = data$dist_to_nuc, y = data$score, type = "l", col = col)
  }
  
  # calculate axis ranges if necessary
  if(is.null(xlim)){
    xlim <- range(c(0, S1_seq_1$dist_to_nuc, S1_seq_2$dist_to_nuc, S1_seq_4$dist_to_nuc,
                    MNase_seq_0$dist_to_nuc, MNase_seq_1$dist_to_nuc, MNase_seq_2$dist_to_nuc, MNase_seq_4$dist_to_nuc))
    cat("\nCalculated xlim = (", xlim[1], ", ", xlim[2], ")", sep = "")
  }
  
  if(is.null(ylim_MNase)){
    ylim_MNase <- range(c(0, MNase_seq_0$score[MNase_seq_0$dist_to_nuc >= xlim[1] & MNase_seq_0$dist_to_nuc <= xlim[2]], 
                          MNase_seq_1$score[MNase_seq_1$dist_to_nuc >= xlim[1] & MNase_seq_1$dist_to_nuc <= xlim[2]], 
                          MNase_seq_2$score[MNase_seq_2$dist_to_nuc >= xlim[1] & MNase_seq_2$dist_to_nuc <= xlim[2]], 
                          MNase_seq_4$score[MNase_seq_4$dist_to_nuc >= xlim[1] & MNase_seq_4$dist_to_nuc <= xlim[2]]))
    cat("\nCalculated ylim_MNase = (", ylim_MNase[1], ", ", ylim_MNase[2], ")", sep = "")
  }
  
  # start empty plot
  plot(x = NA, y = NA, xlim = xlim, ylim = ylim_MNase, axes = FALSE, ann = FALSE)
  
  # add MNase-seq data
  plot_MNase_seq(data = MNase_seq_0, col = MNase_seq_cols[1])
  plot_MNase_seq(data = MNase_seq_1, col = MNase_seq_cols[2])
  plot_MNase_seq(data = MNase_seq_2, col = MNase_seq_cols[3])
  plot_MNase_seq(data = MNase_seq_4, col = MNase_seq_cols[4])
  
  if(!is.null(nuc_marks)){
    # add nucleosome marks
    segments(x0 = (nuc_marks - 1) * nuc_dist, y0 = ylim_MNase[1], 
             x1 = (nuc_marks - 1) * nuc_dist, y1 = ylim_MNase[2], 
             col = nuc_mark_col, lty = "dashed")
  }
  
  if(!is.null(S1_seq_1) | !is.null(S1_seq_2) | !is.null(S1_seq_4)){
    # add S1-seq data as new plot
    par(new = TRUE)
    
    # calculate axis ranges if necessary
    if(is.null(ylim_S1)){
      ylim_S1 <- range(c(0, S1_seq_1$score[S1_seq_1$dist_to_nuc >= xlim[1] & S1_seq_1$dist_to_nuc <= xlim[2]], 
                         S1_seq_2$score[S1_seq_2$dist_to_nuc >= xlim[1] & S1_seq_2$dist_to_nuc <= xlim[2]], 
                         S1_seq_4$score[S1_seq_4$dist_to_nuc >= xlim[1] & S1_seq_4$dist_to_nuc <= xlim[2]]))
      cat("\nCalculated ylim_S1 = (", ylim_S1[1], ", ", ylim_S1[2], ")", sep = "")
    }
    
    # start empty plot
    plot(x = NA, y = NA, xlim = xlim, ylim = ylim_S1, axes = FALSE, ann = FALSE)
    
    # add MNase-seq data
    if(!is.null(S1_seq_1)){plot_S1_seq(data = S1_seq_1, col = S1_seq_cols[1])}
    if(!is.null(S1_seq_2)){plot_S1_seq(data = S1_seq_2, col = S1_seq_cols[2])}
    if(!is.null(S1_seq_4)){plot_S1_seq(data = S1_seq_4, col = S1_seq_cols[3])}
  }
  
}


# define DSBs to be analyzed ==============================================
DSBs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
DSBs <- DSBs[DSBs$DSB_kinetics_rank < 17] # exclude very slowly formed DSBs
roi <- DSB_regions(DSBs = DSBs, region_width = 4000, up_rev_down_fw = TRUE)

# LSY4518-13B =============================================================

# load and process data ---------------------------------------------------

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
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

add_distance_to_nuc_helper(object_name = "LSY4518_13B_1h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_2h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_4h_S1_seq")

add_distance_to_nuc_helper(object_name = "LSY4518_13B_0h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_1h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_2h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY4518_13B_4h_MNase_seq")


# average and smoothen ----------------------------------------------------
k <- 51

# MNase-seq data
for(t in c(0, 1, 2, 4)){
  tmp <- aggregate(score ~ dist_to_nuc, 
                   data = get(paste0("LSY4518_13B_", t, "h_MNase_seq")), FUN = mean)
  tmp$score <- runmed(x = tmp$score, k = k)
  assign(x = paste0("MNase_seq_", t), value = tmp)
}

# S1-seq data
for(t in c(1, 2, 4)){
  tmp <- aggregate(score ~ dist_to_nuc, 
                   data = get(paste0("LSY4518_13B_", t, "h_S1_seq")), FUN = mean)
  tmp$score <- runmed(x = tmp$score, k = k)
  assign(x = paste0("S1_seq_", t), value = tmp)
}


# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/MRE11/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)


pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))
  
  plotting_function(MNase_seq_0 = MNase_seq_0, MNase_seq_1 = MNase_seq_1, MNase_seq_2 = MNase_seq_2, MNase_seq_4 = MNase_seq_4,
                    S1_seq_1 = S1_seq_1, S1_seq_2 = S1_seq_2, S1_seq_4 = S1_seq_4,
                    xlim = c(-160, 990))  

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY4518-13B.pdf"))


# LSY5415 =============================================================

# load and process data ---------------------------------------------------

load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_1h_S1_seq <- subsetByIntersect(subject = LSY5415_1h_S1_seq, query = roi)
LSY5415_2h_S1_seq <- subsetByIntersect(subject = LSY5415_2h_S1_seq, query = roi)
LSY5415_4h_S1_seq <- subsetByIntersect(subject = LSY5415_4h_S1_seq, query = roi)

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

add_distance_to_nuc_helper(object_name = "LSY5415_1h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY5415_2h_S1_seq")
add_distance_to_nuc_helper(object_name = "LSY5415_4h_S1_seq")

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
  tmp$score <- runmed(x = -tmp$score, k = k)
  assign(x = paste0("MNase_seq_", t), value = tmp)
}

# S1-seq data
for(t in c(1, 2, 4)){
  tmp <- aggregate(score ~ dist_to_nuc, 
                   data = get(paste0("LSY5415_", t, "h_S1_seq")), FUN = mean)
  tmp$score <- runmed(x = -tmp$score, k = k)
  assign(x = paste0("S1_seq_", t), value = tmp)
}


# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/MRE11/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)


pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

plotting_function(MNase_seq_0 = MNase_seq_0, MNase_seq_1 = MNase_seq_1, MNase_seq_2 = MNase_seq_2, MNase_seq_4 = MNase_seq_4,
                  S1_seq_1 = S1_seq_1, S1_seq_2 = S1_seq_2, S1_seq_4 = S1_seq_4,
                  xlim = c(-160, 990))  

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5415.pdf"))


# LSY5935 =============================================================

# load and process data ---------------------------------------------------

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5935_MNase-seq.RData")
LSY5935_0h_MNase_seq <- subsetByIntersect(subject = LSY5935_0h_MNase_seq, query = roi)
LSY5935_1h_MNase_seq <- subsetByIntersect(subject = LSY5935_1h_MNase_seq, query = roi)
LSY5935_2h_MNase_seq <- subsetByIntersect(subject = LSY5935_2h_MNase_seq, query = roi)
LSY5935_4h_MNase_seq <- subsetByIntersect(subject = LSY5935_4h_MNase_seq, query = roi)

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5935_nucleosome_positions.RData")
LSY5935_0h_nucleosome_positions <- subsetByIntersect(subject = LSY5935_0h_nucleosome_positions, query = roi)
LSY5935_1h_nucleosome_positions <- subsetByIntersect(subject = LSY5935_1h_nucleosome_positions, query = roi)
LSY5935_2h_nucleosome_positions <- subsetByIntersect(subject = LSY5935_2h_nucleosome_positions, query = roi)
LSY5935_4h_nucleosome_positions <- subsetByIntersect(subject = LSY5935_4h_nucleosome_positions, query = roi)


# add distance to DSB-proximal nucleosome ---------------------------------

# helper function (using the same nuc_pos object [t=0] for all seq_data objects)
add_distance_to_nuc_helper <- function(object_name){
  tmp <- GRanges()
  for(r in 1:length(roi)){
    seq_data <- subsetByIntersect(subject = get(object_name), query = roi[r])
    nuc_pos <- subsetByIntersect(subject = LSY5935_0h_nucleosome_positions, query = roi[r])
    tmp <- c(tmp, add_distance_to_nuc(seq_data = seq_data, nuc_pos = nuc_pos, nuc = 1))
  }
  assign(x = object_name, value = tmp, envir = .GlobalEnv)
}

add_distance_to_nuc_helper(object_name = "LSY5935_0h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5935_1h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5935_2h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5935_4h_MNase_seq")


# average and smoothen ----------------------------------------------------
k <- 51

# MNase-seq data
for(t in c(0, 1, 2, 4)){
  tmp <- aggregate(score ~ dist_to_nuc, 
                   data = get(paste0("LSY5935_", t, "h_MNase_seq")), FUN = mean)
  tmp$score <- runmed(x = tmp$score, k = k)
  assign(x = paste0("MNase_seq_", t), value = tmp)
}


# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/mre11-nd/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)


pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

plotting_function(MNase_seq_0 = MNase_seq_0, MNase_seq_1 = MNase_seq_1, MNase_seq_2 = MNase_seq_2, MNase_seq_4 = MNase_seq_4,
                  xlim = c(-160, 990))  # , nuc_marks = 1:6  

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5935.pdf"))

# save coverage at +1 Nuc
tmp <- data.frame(time = c(0, 1, 2, 4),
                  score = c(MNase_seq_0$score[MNase_seq_0$dist_to_nuc == 0],
                            MNase_seq_1$score[MNase_seq_0$dist_to_nuc == 0],
                            MNase_seq_2$score[MNase_seq_0$dist_to_nuc == 0],
                            MNase_seq_4$score[MNase_seq_0$dist_to_nuc == 0]))

write.table(x = tmp, file = paste0(plot_dir, "/Nuc+1_cov_LSY5935.txt"), row.names = FALSE)


# LSY5934 =============================================================

# load and process data ---------------------------------------------------

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5934_MNase-seq.RData")
LSY5934_0h_MNase_seq <- subsetByIntersect(subject = LSY5934_0h_MNase_seq, query = roi)
LSY5934_1h_MNase_seq <- subsetByIntersect(subject = LSY5934_1h_MNase_seq, query = roi)
LSY5934_2h_MNase_seq <- subsetByIntersect(subject = LSY5934_2h_MNase_seq, query = roi)
LSY5934_4h_MNase_seq <- subsetByIntersect(subject = LSY5934_4h_MNase_seq, query = roi)

load(file = "../../Rep_merged/MNase-seq/03_Processed_data/Nucleosome_positions/LSY5934_nucleosome_positions.RData")
LSY5934_0h_nucleosome_positions <- subsetByIntersect(subject = LSY5934_0h_nucleosome_positions, query = roi)
LSY5934_1h_nucleosome_positions <- subsetByIntersect(subject = LSY5934_1h_nucleosome_positions, query = roi)
LSY5934_2h_nucleosome_positions <- subsetByIntersect(subject = LSY5934_2h_nucleosome_positions, query = roi)
LSY5934_4h_nucleosome_positions <- subsetByIntersect(subject = LSY5934_4h_nucleosome_positions, query = roi)


# add distance to DSB-proximal nucleosome ---------------------------------

# helper function (using the same nuc_pos object [t=0] for all seq_data objects)
add_distance_to_nuc_helper <- function(object_name){
  tmp <- GRanges()
  for(r in 1:length(roi)){
    seq_data <- subsetByIntersect(subject = get(object_name), query = roi[r])
    nuc_pos <- subsetByIntersect(subject = LSY5934_0h_nucleosome_positions, query = roi[r])
    tmp <- c(tmp, add_distance_to_nuc(seq_data = seq_data, nuc_pos = nuc_pos, nuc = 1))
  }
  assign(x = object_name, value = tmp, envir = .GlobalEnv)
}

add_distance_to_nuc_helper(object_name = "LSY5934_0h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5934_1h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5934_2h_MNase_seq")
add_distance_to_nuc_helper(object_name = "LSY5934_4h_MNase_seq")


# average and smoothen ----------------------------------------------------
k <- 51

# MNase-seq data
for(t in c(0, 1, 2, 4)){
  tmp <- aggregate(score ~ dist_to_nuc, 
                   data = get(paste0("LSY5934_", t, "h_MNase_seq")), FUN = mean)
  tmp$score <- runmed(x = -tmp$score, k = k)
  assign(x = paste0("MNase_seq_", t), value = tmp)
}


# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/mre11-nd/Avgs_algned_at_DSB_prxml_nuc"
dir.create(path = plot_dir, showWarnings = FALSE)


pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = rep(0,4), oma = rep(0, 4))

plotting_function(MNase_seq_0 = MNase_seq_0, MNase_seq_1 = MNase_seq_1, MNase_seq_2 = MNase_seq_2, MNase_seq_4 = MNase_seq_4,
                  xlim = c(-160, 990))

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/LSY5934.pdf"))

# save coverage at +1 Nuc
tmp <- data.frame(time = c(0, 1, 2, 4),
                  score = c(MNase_seq_0$score[MNase_seq_0$dist_to_nuc == 0],
                            MNase_seq_1$score[MNase_seq_0$dist_to_nuc == 0],
                            MNase_seq_2$score[MNase_seq_0$dist_to_nuc == 0],
                            MNase_seq_4$score[MNase_seq_0$dist_to_nuc == 0]))

tmp$score <- -tmp$score

write.table(x = tmp, file = paste0(plot_dir, "/Nuc+1_cov_LSY5934.txt"), row.names = FALSE)
