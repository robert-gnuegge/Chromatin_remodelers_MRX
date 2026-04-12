# info --------------------------------------------------------------------
# purpose: analyze meltability impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 06/22/25
# last modified: 06/22/25

# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(BSgenome)
library(DNAshapeR)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

# read modified S. cerevisiae genome
genome <- import(con = "../../../Reference_genome/03_Processed_data/S288C_R64-4-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")


# function definitions ----------------------------------------------------

# find resection extend
# argument: GRanges
# result: GRanges
crop_to_resection_extend <- function(GRanges, DSBs, max_DSB_dist = 2500, threshold = 0.95){
  out <- DSBs
  out$DSB_id <- as.character(DSBs)
  for(n in 1:length(out)){
    # find resection tract end on - strand
    roi <- DSBs[n]
    strand(roi) <- "-"
    roi <- resize(x = roi, width = max_DSB_dist, fix = "start")
    tmp <- subsetByIntersect(subject = GRanges, query = roi)
    tmp <- rev(tmp)  # for increasing distance_to_DSB values
    idx <- which.max(cumsum(tmp$score) >= threshold * max(cumsum(tmp$score)))
    start(out[n]) <- start(tmp[idx])
    # find tract end on + strand
    roi <- DSBs[n]
    strand(roi) <- "+"
    roi <- resize(x = roi, width = max_DSB_dist, fix = "start")
    tmp <- subsetByIntersect(subject = GRanges, query = roi)
    idx <- which.max(cumsum(tmp$score) >= threshold * max(cumsum(tmp$score)))
    end(out[n]) <- end(tmp[idx])
  }
  subsetByIntersect(subject = GRanges, query = out)
}

# calculate average meltability profiles
# argument: GRanges, logical, character
# result: numeric vector
average_shape <- function(GRanges_DNA_shape, GRanges_S1, bg = FALSE, shape = "MGW"){
  hits <- findOverlaps(query = GRanges_S1, subject = GRanges_DNA_shape)
  S1 <- GRanges_S1[queryHits(hits)]
  DNA_shape <- GRanges_DNA_shape[subjectHits(hits)]
  tmp <- mcols(DNA_shape)[shape]
  if(bg){
    out <- apply(X = as.matrix(tmp), MARGIN = 2, FUN = mean)
  }else{
    out <- apply(X = as.matrix(tmp) * S1$score, MARGIN = 2, FUN = sum) / sum(S1$score)
  }
  return(out)
}


# # Calculate shape for each nt in DSB regions ========================
# 
# # get GRanges around all potential nick sites
# DSBs <- SrfIcs[- c(9, 17)]  # exclude SrfIcs in duplicated region
# pos <- as_nt_resolved_GRanges(DSB_regions(DSBs = DSBs, region_width = 3000, up_rev_down_fw = TRUE))
# shape_width <- 70
# tiles <- resize(x = pos, width = shape_width, fix = "center")
# # remove tiles that overlap with DSBs 
# hits <- findOverlaps(query = DSBs, subject = tiles)
# tiles <- tiles[-subjectHits(hits)]
# pos <- pos[-subjectHits(hits)]
# 
# # get sequences
# sequences <- as.character(getSeq(x = genome, names = tiles))
# 
# # save as fasta files (required for using getShape command) 
# names(sequences) <- paste0(seqnames(pos), ":", start(pos))  # add names (will be fasta entry IDs)
# sequences <- DNAStringSet(sequences)
# writeXStringSet(x = sequences, filepath = "03_Processed_data/Nt_sequences_70nt_windows_around_SrfIcs.fasta")
# 
# # calculate DNA shape features for each sequence
# res <- getShape(filename = "03_Processed_data/Nt_sequences_70nt_windows_around_SrfIcs.fasta")
# lapply(X = res, FUN = dim)
# 
# # save results
# DNA_shape <- pos
# DNA_shape$MGW <- res$MGW
# DNA_shape$Roll <- cbind(res$Roll, NA)  # add NA column for equal ncol
# DNA_shape$HelT <- cbind(res$HelT, NA)  # add NA column for equal ncol
# DNA_shape$ProT <- res$ProT
# DNA_shape$EP <- res$EP
# 
# save(DNA_shape, file = "03_Processed_data/Shape/DNA_shape_around_SrfIcs.RData")

# load shape data
load(file = "03_Processed_data/Shape/DNA_shape_around_SrfIcs.RData")


# process data for plotting ===============================================

# read unnormalized S1-seq coverage data
load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq_unnormalized.RData")
load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq_unnormalized.RData")

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
WT_1h <- crop_to_resection_extend(GRanges = LSY4518_13B_1h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
WT_2h <- crop_to_resection_extend(GRanges = LSY4518_13B_2h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
WT_4h <- crop_to_resection_extend(GRanges = LSY4518_13B_4h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])

fun30_1h <- crop_to_resection_extend(GRanges = LSY5415_1h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
fun30_2h <- crop_to_resection_extend(GRanges = LSY5415_2h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
fun30_4h <- crop_to_resection_extend(GRanges = LSY5415_4h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])

# remove SrfIcs sequence regions
# they give a lot of S1-seq signal, but derive from unprocessed DSBs and not from MRX nicking
DSBs <- DSB_regions(DSBs = SrfIcs, region_width = 8)
unique(getSeq(x = genome, names = DSBs))  # they are indeed all the SrfIcs sequence

WT_1h <- WT_1h[-subjectHits(findOverlaps(query = DSBs, subject = WT_1h))]
WT_2h <- WT_2h[-subjectHits(findOverlaps(query = DSBs, subject = WT_2h))]
WT_4h <- WT_4h[-subjectHits(findOverlaps(query = DSBs, subject = WT_4h))]

fun30_1h <- fun30_1h[-subjectHits(findOverlaps(query = DSBs, subject = fun30_1h))]
fun30_2h <- fun30_2h[-subjectHits(findOverlaps(query = DSBs, subject = fun30_2h))]
fun30_4h <- fun30_4h[-subjectHits(findOverlaps(query = DSBs, subject = fun30_4h))]

# calculate average MGW profiles
WT <- data.frame(MGW = average_shape(GRanges_DNA_shape = DNA_shape, GRanges_S1 = c(WT_1h, WT_2h, WT_4h), shape = "MGW"))
WT_bg <- data.frame(MGW = average_shape(GRanges_DNA_shape = DNA_shape, GRanges_S1 = c(WT_1h, WT_2h, WT_4h), shape = "MGW", bg = TRUE))

fun30 <- data.frame(MGW = average_shape(GRanges_DNA_shape = DNA_shape, GRanges_S1 = c(fun30_1h, fun30_2h, fun30_4h), shape = "MGW"))
fun30_bg <- data.frame(MGW = average_shape(GRanges_DNA_shape = DNA_shape, GRanges_S1 = c(fun30_1h, fun30_2h, fun30_4h), shape = "MGW", bg = TRUE))

# plotting ================================================================

dir.create(path = "04_Plots/Shape_impact/")

# plotting helper function ------------------------------------------------

# add custom x axis
add_custom_x_axis <- function(){
  at <- c(0:6 * 5 + 1, 8:14 * 5)
  labels <- c(-7:-1, 1:7) * 5
  axis(side = 1, at = 1:70, labels = FALSE, tcl = -0.15)
  axis(side = 1, at = c(at, 35, 36), labels = NA, tcl = -0.3)
  axis(side = 1, at = at, labels = labels, tick = FALSE)
}

# WT -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.2, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(WT$MGW, WT_bg$MGW), na.rm = TRUE),
     xlab = NA, ylab = expression("MGW ["*ring(A)*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site [nt]", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = WT$MGW, type = "l")
points(x = 1:70, y = WT_bg$MGW, type = "l", lty = "dashed")

legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = "04_Plots/Shape_impact/WT_MWG_preference.pdf")

# fun30 -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.2, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(fun30$MGW, fun30_bg$MGW), na.rm = TRUE),
     xlab = NA, ylab = expression("MGW ["*ring(A)*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site [nt]", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = fun30$MGW, type = "l")
points(x = 1:70, y = fun30_bg$MGW, type = "l", lty = "dashed")

legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = "04_Plots/Shape_impact/fun30_MWG_preference.pdf")

# both -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.2, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(WT$MGW, WT_bg$MGW, fun30$MGW, fun30_bg$MGW), na.rm = TRUE),
     xlab = NA, ylab = expression("MGW ["*ring(A)*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site (nt)", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = WT$MGW, type = "l")
points(x = 1:70, y = fun30$MGW, type = "l", col = JFly_colors[2])
points(x = 1:70, y = WT_bg$MGW, type = "l", lty = "dashed")
points(x = 1:70, y = fun30_bg$MGW, type = "l", lty = "dashed", col = JFly_colors[2])

Leg_txt <- c(expression(italic("FUN30")~"background"), expression(italic("FUN30")~"observed"))
legend(x = "bottomleft", legend = Leg_txt, lty = c("dashed", "solid"), bty = "n", inset = 0.03)
Leg_txt <- c(expression(italic("fun30"*Delta)~"background"), expression(italic("fun30"*Delta)~"observed"))
legend(x = "bottomright", legend = Leg_txt, lty = c("dashed", "solid"), col = JFly_colors[2],  bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = "04_Plots/Shape_impact/WT_and_fun30_MWG_preference.pdf")
