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

# iterate through samples
for(strain in c("LSY4518-13B", "LSY5415", "LSY5934", "LSY5935")[1]){
  for(t in c(0, 1, 2, 4)[1]){
    
    samples <- paste0(strain, "_", t, "h")
    
    for(sample in samples){
      
      cat("\n\nProcessing ", sample, "...", sep = "")
      
      # read BAM file -----------------------------------------------------------
      cat("\nReading BAM file...")
      tmp <- readGAlignmentPairs(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = ScanBamParam(tag = c("AS", "YS", "NM"), what = "isize"))
      # AS: alignment score for first mate (max. read length * match bonus [--ma, default: 2])
      # YS: alignment score for opposite mate
      # NM: edit distance
      # isize: insert size
      all_mapped <- length(tmp)  # for mapping statistics
      cat(all_mapped, "alignments read.")
      
      # Remove alignments with too large insert size ----------------------------
      cat("\nRemoving alignments with insert size >250 bp...")
      tmp <- tmp[mcols(first(tmp))$isize <= 250]
      mapped <- length(tmp)
      cat(" kept ", mapped, " alignments (",  round(100 * mapped / all_mapped, digits = 2), "%).", sep = "")
      
      # Convert to GRanges ------------------------------------------------------
      # conversion to GRanges is necessary for trimming (setting start and end is not possible for GAlignmentPairs class)
      # and for coverage calculation to consider the complete insert sequence
      cat("\nConverting to GRanges...")
      tmp <- GRanges(tmp)
      
      # Adjust to insert size = 100 bp ------------------------------------------
      cat("\nTrimming to max. 100 bp insert sizes...")
      tmp[width(tmp) > 100] <- resize(x = tmp[width(tmp) > 100], width = 100, fix = "center")

      # calculate MNase-seq coverage --------------------------------------------
      cat("\nCalculating coverage...")
      tmp_coverage <- GRanges(coverage(tmp))
      tmp_coverage$score <- tmp_coverage$score / length(tmp_coverage) * 1e6  # convert to RPM
      assign(x = paste0(dash_to_underscore(sample), "_MNase_seq"), value = tmp_coverage)
    }
    
  }
  
  # save and free up storage
  cat("\n\nSaving coverage data...")
  save(list = paste0(dash_to_underscore(samples), "_MNase_seq"), file = paste0(save_dir, "/", strain, "_MNase-seq_trimmed_100bp.RData"))
  rm(samples)
  gc()
  
}
