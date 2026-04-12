# info --------------------------------------------------------------------
# purpose: analyze meltability impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 06/20/25
# last modified: 06/22/25

# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(BSgenome)
library(rmelting)
library(parallel)
# library(tictoc)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

# read modified S. cerevisiae genome
genome <- import(con = "../../../Reference_genome/03_Processed_data/S288C_R64-4-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")


# function definitions ----------------------------------------------------

# calculate meltability for DNA sequence
# argument: character string
# result: data.frame with rows enthalpy.J, entropy.J, and temperature.C
calc_melting <- function(sequence, nucleic.acid.conc = 5.54e-10, Na.conc = 0.02, K.conc = 0.3, Mg.conc = 0.002, hybridisation.type = "dnadna"){
  tmp <- melting(sequence = sequence, nucleic.acid.conc = nucleic.acid.conc, 
                 Na.conc = Na.conc, K.conc = K.conc, Mg.conc = Mg.conc, hybridisation.type = hybridisation.type)
  out <- unlist(tmp$Results[5])
  return(out)
}

# extract meltability profiles
# argument: GRanges
# result: list with GRanges and matrix
extract_T_m_profiles <- function(GRanges, GRanges_w_T_m, width, DSBs){
  tiles <- resize(x = GRanges, width = width, fix = "center")
  hits <- findOverlaps(query = DSBs, subject = tiles)
  if(length(hits) > 0){
    tiles <- tiles[-subjectHits(hits)]  # remove tiles that overlap DSBs
    GRanges <- GRanges[-subjectHits(hits)]  # also remove corresponding GRanges
  }
  profiles <- matrix(data = 0, nrow = length(tiles), ncol = width)
  pb <- txtProgressBar(min = 0, max = length(tiles), initial = 0, style = 3, width = 76)  # show progress bar
  for(n in 1:length(tiles)){
    hits <- findOverlaps(query = tiles[n], subject = GRanges_w_T_m)
    profiles[n, ] <- GRanges_w_T_m[subjectHits(hits)]$T_m
    setTxtProgressBar(pb, n)
  }
  GRanges$T_m_profile <- profiles
  return(GRanges)
}

# calculate average meltability profiles
# argument: GRanges, logical, character
# result: numeric vector
average_Tm_profile <- function(GRanges_Tm, GRanges_S1, bg = FALSE){
  hits <- findOverlaps(query = GRanges_S1, subject = GRanges_Tm)
  S1 <- GRanges_S1[queryHits(hits)]
  Tm <- GRanges_Tm[subjectHits(hits)]
  tmp <- Tm$T_m_profile
  if(bg){
    out <- apply(X = as.matrix(tmp), MARGIN = 2, FUN = mean)
  }else{
    out <- apply(X = as.matrix(tmp) * S1$score, MARGIN = 2, FUN = sum) / sum(S1$score)
  }
  return(out)
}

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

# # Calculate meltability for each nt in DSB regions ========================
# 
# # get GRanges around all potential nick sites
# DSBs <- SrfIcs[- c(9, 17)]  # exclude SrfIcs in duplicated region
# pos <- as_nt_resolved_GRanges(DSB_regions(DSBs = DSBs, region_width = 3000, up_rev_down_fw = TRUE))
# melting_window <- 15  # define how long the molten DNA stretch is
# tiles <- resize(x = pos, width = melting_window, fix = "center")
# 
# # get sequences
# sequences <- as.character(getSeq(x = genome, names = tiles))
# 
# # calculate melting temperature for each sequence
# # let's use parallel computation for speed-up
# tictoc::tic(msg = "Melting calculation")
# n_cores <- detectCores()
# cl <- makeCluster(spec = n_cores - 2)  # start cluster
# clusterEvalQ(cl = cl, expr = {library(rmelting)})  # load required libraries
# rmelting_out <- parSapply(cl = cl, X = sequences, FUN = calc_melting)
# stopCluster(cl)  # stop cluster
# tictoc::toc()  # Melting calculation: 4705.847 sec elapsed
# 
# # construct GRanges object with calculated melting parameters for each nt
# Melting_15bp <- pos
# Melting_15bp$T_m <- rmelting_out
# 
# # derive T_m profiles around all potential nick sites
# pos <- as_nt_resolved_GRanges(DSB_regions(DSBs = DSBs, region_width = 2900, up_rev_down_fw = TRUE))
# Melting_profiles_15bp <- extract_T_m_profiles(GRanges = pos, GRanges_w_T_m = Melting_15bp, width = 70, DSBs = DSBs)
# 
# # save results
# save_dir <- "03_Processed_data/Melting/"
# dir.create(path = save_dir)
# save(Melting_15bp, Melting_profiles_15bp, file = paste0(save_dir, "Melting_15bp.RData"))

# load melting data
load(file = "03_Processed_data/Melting/Melting_15bp.RData")


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

# calculate average Tm profiles
WT <- data.frame(Tm_15bp = average_Tm_profile(GRanges_Tm = Melting_profiles_15bp, GRanges_S1 = c(WT_1h, WT_2h, WT_4h)))
WT_bg <- data.frame(Tm_15bp = average_Tm_profile(GRanges_Tm = Melting_profiles_15bp, GRanges_S1 = c(WT_1h, WT_2h, WT_4h), bg = TRUE))

fun30 <- data.frame(Tm_15bp = average_Tm_profile(GRanges_Tm = Melting_profiles_15bp, GRanges_S1 = c(fun30_1h, fun30_2h, fun30_4h)))
fun30_bg <- data.frame(Tm_15bp = average_Tm_profile(GRanges_Tm = Melting_profiles_15bp, GRanges_S1 = c(fun30_1h, fun30_2h, fun30_4h), bg = TRUE))

# plotting ================================================================

dir.create(path = "04_Plots/Melting_impact/")

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
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(WT$Tm_15bp, WT_bg$Tm_15bp)),
     xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site [nt]", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = WT$Tm_15bp, type = "l")
points(x = 1:70, y = WT_bg$Tm_15bp, type = "l", lty = "dashed")

legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = "04_Plots/Melting_impact/WT_melting_preference_15bp.pdf")

# fun30 -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(fun30$Tm_15bp, fun30_bg$Tm_15bp)),
     xlab = NA, ylab = expression("Melting Temperature ["*degree*"C"*"]"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site [nt]", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = fun30$Tm_15bp, type = "l")
points(x = 1:70, y = fun30_bg$Tm_15bp, type = "l", lty = "dashed")

legend(x = "bottomleft", legend = c("Background", "Observed"), lty = c("dashed", "solid"), bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = "04_Plots/Melting_impact/fun30_melting_preference_15bp.pdf")

# both -------------------------------------------------------------------
pdf(file = "tmp.pdf", width = 8, height = 2.75)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 0.5, 4, 2), tcl = -0.25, mgp = c(2.5, 0.5, 0), las = 1)

plot(x = 1:70, y = rep(NA, 70), ylim = range(c(WT$Tm_15bp, WT_bg$Tm_15bp, fun30$Tm_15bp, fun30_bg$Tm_15bp)),
     xlab = NA, ylab = expression("Melting Temperature ("*degree*"C"*")"), xaxt = "n")
add_custom_x_axis()
title(xlab = "Distance from Nick Site (nt)", line = 1.75)
abline(v = 35.5, col = "gray")

points(x = 1:70, y = WT$Tm_15bp, type = "l")
points(x = 1:70, y = fun30$Tm_15bp, type = "l", col = JFly_colors[2])
points(x = 1:70, y = WT_bg$Tm_15bp, type = "l", lty = "dashed")
points(x = 1:70, y = fun30_bg$Tm_15bp, type = "l", lty = "dashed", col = JFly_colors[2])

Leg_txt <- c(expression(italic("FUN30")~"background"), expression(italic("FUN30")~"observed"))
legend(x = "bottomleft", legend = Leg_txt, lty = c("dashed", "solid"), bty = "n", inset = 0.03)
Leg_txt <- c(expression(italic("fun30"*Delta)~"background"), expression(italic("fun30"*Delta)~"observed"))
legend(x = "bottomright", legend = Leg_txt, lty = c("dashed", "solid"), col = JFly_colors[2],  bty = "n", inset = 0.03)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = "04_Plots/Melting_impact/WT_and_fun30_melting_preference_15bp.pdf")
