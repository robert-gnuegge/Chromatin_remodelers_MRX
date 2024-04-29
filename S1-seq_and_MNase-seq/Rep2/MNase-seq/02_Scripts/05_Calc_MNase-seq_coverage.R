# info --------------------------------------------------------------------
# purpose: filter alignments, calculate MNase-seq coverage, and extract mapping stats
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 04/25/24
# last modified: 04/29/24

# load libraries ----------------------------------------------------------
library(GenomicAlignments)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")


# process all samples =====================================================

# file base paths
BAM_dir <- "03_Processed_data/BAM/"
# BAM_dir <- "/media/robert/Elements/Deep_sequencing_data/24-04-24-MNase-seq/BAM"
save_dir <- "03_Processed_data/MNase-seq_coverage"
dir.create(path = save_dir, showWarnings = FALSE)
plot_dir_edit <- "04_Plots/Mapping_edit_distance_distributions"
dir.create(path = plot_dir_edit, showWarnings = FALSE)
plot_dir_insert <- "04_Plots/Insert_size_distributions"
dir.create(path = plot_dir_insert, showWarnings = FALSE)

# initialize data.frames for result recording
mapping_stats <- data.frame()

# iterate through samples
samples <- paste0(rep(c("LSY4518-13B", "LSY5415", "LSY5934", "LSY5935"), rep(4, 4)), "_", c(0, 1, 2, 4), "h")

for(sample in samples){
  
  cat("\n\nProcessing ", sample, "...", sep = "")
  
  # read BAM file -----------------------------------------------------------
  cat("\nReading BAM file...")
  tmp <- readGAlignmentPairs(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = ScanBamParam(tag = c("NM")))
  # NM: edit distance
  total_algns <- length(tmp) # for mapping statistics
  cat(" read", total_algns, "alignment pairs.")
  
  # plot edit distance distribution -----------------------------------------
  cat("\nPlotting edit distance distribution...")
  edit_dist <- c(mcols(first(tmp))$NM, mcols(last(tmp))$NM)
  h <- hist(x = edit_dist, breaks = 0:max(edit_dist), plot = FALSE)  # use hist function to calculate densities (fractions)
  # plot fractions of alignments with edit distance = 0, 1, 2, ..., 5, >5
  pdf(file = "tmp.pdf", width=2.5, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.2, 0.9, 3.6, 2.1), las = 1, tcl = -0.3, mgp = c(2.3, 0.5, 0))
  bp <- barplot(height = c(h$density[1:5], sum(h$density[6:length(h$density)])), ylim = c(0, 1),
                names.arg = c(as.character(0:4), expression("">="5")), ylab = "Fraction of Alignments", xlab = NA)  # plot density values as bar plot
  title(xlab = "Edit Distance", line = 2)
  # add edit distance = 0 percentage
  text(x = bp[1, 1], y = 0.5, labels = paste0(round(x = 100 * h$density[1], digits = 2), "%"), srt = 90, col = "red")
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir_edit, "/", sample, ".pdf"))
  
  # process alignments ------------------------------------------------------
  cat("\nConverting to GRanges...")
  # conversion to GRanges is necessary for trimming
  # and for coverage calculation considering the complete insert sequence
  tmp <- GRanges(tmp)
  
  # Remove alignments with too large insert size (probably disomes) ---------
  cat("\nRemoving alignments with insert size >250 bp...")
  old_length <- length(tmp)
  tmp <- tmp[width(tmp) <= 250]
  new_length <- length(tmp)
  greater_250 <- old_length - new_length
  cat(" kept ",  round(100 * new_length / old_length, digits = 2), "% of initial alignments.", sep = "")
  
  # plot insert size distribution -----------------------------------------
  cat("\nPlotting insert size distribution...")
  insert_size <- width(tmp)
  h <- hist(x = insert_size, breaks = 100, plot = FALSE)  # to calculate ylim of histogram
  # plot distribution and median
  pdf(file = "tmp.pdf", width=2.75, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.2, 0.2, 4.1, 1.6), las = 1, tcl = -0.3, mgp = c(1.8, 0.5, 0))
  hist(x = insert_size, breaks = 100, xlim = c(50, 250), probability = TRUE,
       ylim = c(0, 1.2 * max(h$density)), xlab = "Insert size (bp)", ylab = NA, main = NA)
  title(ylab = "Probability", line = 3)
  med <- median(insert_size)
  segments(x0 = med, y0 = 0, y1 = 1.1 * max(h$density), col = "red")
  text(x = med, y = 1.1 * max(h$density), labels = med, pos = 3, offset = 0.4, col = "red")
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir_insert, "/", sample, ".pdf"))

  # Trim to "max_insert" bp (for sharper nucleosome peaks) ------------------
  max_insert <- 117
  cat("\nTrimming to max.", max_insert, "bp insert size...")
  idx <- width(tmp) > max_insert
  tmp[idx] <- resize(x = tmp[idx], width = max_insert, fix = "center")
  
  # calculate MNase-seq coverage --------------------------------------------
  cat("\nCalculating coverage...")
  tmp_coverage <- GRanges(coverage(tmp))
  tmp_coverage$score <- tmp_coverage$score / length(tmp) * 1e6  # convert to RPM
  assign(x = paste0(dash_to_underscore(sample), "_MNase_seq"), value = tmp_coverage)
  
  # record mapping statistics data ------------------------------------------
  mapping_stats <- rbind(mapping_stats,
                         data.frame(sample = sample,
                                    mapped_pairs = total_algns,
                                    greater250 = greater_250,
                                    insert_size_median = med,
                                    alignments_trimmed = sum(idx)))
}

# save data ---------------------------------------------------------------
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY4518-13B", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY4518-13B_MNase-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5415", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY5415_MNase-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5934", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY5934_MNase-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5935", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY5935_MNase-seq.RData"))

write.table(x = mapping_stats, file = paste0(plot_dir_edit, "/Mapping_stats.txt"), row.names = FALSE)
