# info --------------------------------------------------------------------
# purpose: filter alignments, calculate MNase-seq coverage, and extract mapping stats
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/20/24
# last modified: 03/21/24

# load libraries ----------------------------------------------------------
library(GenomicAlignments)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")


# process all samples =====================================================

# file base paths
BAM_dir <- "/media/robert/Elements/Deep_sequencing_data/24-03-20-MNase-seq/BAM"
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
  tmp <- readGAlignmentPairs(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = ScanBamParam(tag = c("AS", "YS", "NM"), what = "isize"))
  # AS: alignment score for first mate (max. read length * match bonus [--ma, default: 2])
  # YS: alignment score for opposite mate
  # NM: edit distance
  # isize: insert size
  all_mapped <- length(tmp) / 2  # paired-end reads, for mapping statistics
  cat(all_mapped, "alignments read.")
  
  # Remove alignments with too large insert size ----------------------------
  cat("\nRemoving alignments with insert size >250 bp...")
  tmp <- tmp[mcols(first(tmp))$isize <= 250]
  mapped <- length(tmp) / 2  # paired-end reads
  cat(" kept ", mapped, " alignments (",  round(100 * mapped / all_mapped, digits = 2), "%).", sep = "")

  # plot edit distance distribution -----------------------------------------
  cat("\nPlotting edit distance distribution...")
  pdf(file = "tmp.pdf", width=2.5, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.6, 0.6, 3.6, 2.1), las = 1, tcl = -0.3, mgp = c(2.5, 0.6, 0))
  # plot fractions of alignments with edit distance = 0, 1, 2, ..., 5, >5
  edit_dist <- c(mcols(first(tmp))$NM, mcols(last(tmp))$NM)
  h <- hist(x = edit_dist, breaks = 0:max(edit_dist), plot = FALSE)  # use hist function to calculate densities (fractions)
  bp <- barplot(height = c(h$density[1:6], sum(h$density[7:length(h$density)])), ylim = c(0, 1),
                names.arg = c(as.character(0:5), ">5"), ylab = "Fraction of Alignments", xlab = "Edit Distance")  # plot density values as bar plot
  # add edit distance <= 3 threshold
  th <- mean(bp[4:5, ])
  abline(v = th, col = "red")
  text(x = 0.95 * th, y = 1, adj = c(1,1), xpd = TRUE, col = "red", labels = paste0(round(x = 100 * ecdf(edit_dist)(3), digits = 2), "%"))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir_edit, "/", sample, ".pdf"))

  # plot insert size distribution -----------------------------------------
  cat("\nPlotting insert size distribution...")
  pdf(file = "tmp.pdf", width=2.75, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2, 0.2, 4.1, 1.6), las = 1, tcl = -0.3, mgp = c(2, 0.5, 0))
  insert_size <- abs(mcols(first(tmp))$isize)
  h <- hist(x = insert_size, breaks = 100, plot = FALSE)  # to calculate ylim of histogram
  hist(x = insert_size, breaks = 100, xlim = c(50, 250), probability = TRUE,
       ylim = c(0, 1.2 * max(h$density)), xlab = "Insert size [bp]", ylab = NA, main = NA)
  title(ylab = "Probability", line = 3)
  med <- median(insert_size)
  segments(x0 = med, y0 = 0, y1 = 1.1 * max(h$density), col = "red")
  text(x = med, y = 1.1 * max(h$density), labels = med, pos = 3, offset = 0.4, col = "red")
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir_insert, "/", sample, ".pdf"))

  # Convert to GRanges ------------------------------------------------------
  # conversion to GRanges is necessary for trimming (setting start and end is not possible for GAlignmentPairs class)
  # and for coverage calculation to consider the complete insert sequence
  cat("\nConverting to GRanges...")
  tmp <- GRanges(tmp)

  # Adjust to median insert size = 147 bp -----------------------------------
  cat("\nMedian insert size is ", med , " bp.", sep = "")
  # trim
  if(med != 147){
      cat(" Adjusting to max insert size = 147 bp...")
      tmp[width(tmp) > 147] <- resize(x = tmp[width(tmp) > 147], width = 147, fix = "center")
  }

  # calculate MNase-seq coverage --------------------------------------------
  cat("\nCalculating coverage...")
  tmp_coverage <- GRanges(coverage(tmp))
  tmp_coverage$score <- tmp_coverage$score / mapped * 1e6  # convert to RPM
  assign(x = paste0(dash_to_underscore(sample), "_MNase_seq"), value = tmp_coverage)
  tmp_trimmed_coverage <- GRanges(coverage(tmp_trimmed))
  tmp_trimmed_coverage$score <- tmp_trimmed_coverage$score / mapped * 1e6  # convert to RPM
  assign(x = paste0(dash_to_underscore(sample), "_MNase_seq_trimmed"), value = tmp_trimmed_coverage)
  
  # record mapping statistics data ------------------------------------------
  mapping_stats <- rbind(mapping_stats,
                         data.frame(sample = sample,
                                    reads_mapped = mapped,
                                    insert_size_median = med))
  
}

# save data ---------------------------------------------------------------
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY4518-13B", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY4518-13B_MNase-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5415", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY5415_MNase-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5934", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY5934_MNase-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5935", x = samples)]), "_MNase_seq"), file = paste0(save_dir, "/LSY5935_MNase-seq.RData"))

save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY4518-13B", x = samples)]), "_MNase_seq_trimmed"), file = paste0(save_dir, "/LSY4518-13B_MNase-seq_trimmed.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5415", x = samples)]), "_MNase_seq_trimmed"), file = paste0(save_dir, "/LSY5415_MNase-seq_trimmed.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5934", x = samples)]), "_MNase_seq_trimmed"), file = paste0(save_dir, "/LSY5934_MNase-seq_trimmed.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5935", x = samples)]), "_MNase_seq_trimmed"), file = paste0(save_dir, "/LSY5935_MNase-seq_trimmed.RData"))

write.table(x = mapping_stats, file = paste0(plot_dir_edit, "/Mapping_stats.txt"), row.names = FALSE)
