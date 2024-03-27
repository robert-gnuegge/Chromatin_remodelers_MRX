# info --------------------------------------------------------------------
# purpose: read alignments in DSB regions, trim to x bp, and calculate coverage
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/27/24
# last modified: 03/27/24

# load libraries ----------------------------------------------------------
library(GenomicAlignments)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")


# define DSB regions ------------------------------------------------------
# DSBs <- SrfIcs[-c(9, 17)]  # exclude duplicated genome regions
roi <- DSB_regions(DSBs = SrfIcs, region_width = 20000)

# process all samples =====================================================

# file base paths
BAM_dir <- "/media/robert/One Touch/Deep_sequencing_data/Robert_Gnuegge/24-03-20-MNase-seq/BAM"
save_dir <- "03_Processed_data/MNase-seq_coverage"
dir.create(path = save_dir, showWarnings = FALSE)


# iterate through strains
for(strain in c("LSY4518-13B", "LSY5415", "LSY5934", "LSY5935")){
  
  # iterate through time points
  for(t in c(0, 1, 2, 4)){
    
    sample <- paste0(strain, "_", t, "h")
    
    cat("\n\nProcessing ", sample, "...", sep = "")
    
    # read BAM file -----------------------------------------------------------
    cat("\nReading BAM file...")
    
    # read alignments in DSB regions
    params <- ScanBamParam(which = roi)  # only read alignments in region of interest
    tmp <- readGAlignmentPairs(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = params)
    
    # count total mapped alignments for normalization
    total_algns <- sum(idxstatsBam(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"))$mapped)
    total_algns <- total_algns / 2  # paired reads
    
    # process alignments ------------------------------------------------------
    cat("\nConverting to GRanges...")
    # conversion to GRanges is necessary for trimming
    # and for coverage calculation considering the complete insert sequence
    tmp <- GRanges(tmp)
    
    cat("\nRemoving alignments with insert size >250 bp...")
    tmp_length <- length(tmp)
    tmp <- tmp[width(tmp) <= 250]
    algns_le_250 <- length(tmp) / tmp_length
    cat(" kept ",  round(100 * algns_le_250, digits = 2), "% of initial alignments.", sep = "")
  
    # iterate through max insert sizes
    for(max_insert in c(147, 127, 117, 107, 97, 87, 74)){
    
      # Trim to "max_insert" bp -------------------------------------------------
      cat("\nTrimming to max.", max_insert, "bp insert size...")
      tmp[width(tmp) > max_insert] <- resize(x = tmp[width(tmp) > max_insert], width = max_insert, fix = "center")
      
      # calculate MNase-seq coverage --------------------------------------------
      cat("\nCalculating coverage...")
      tmp_coverage <- GRanges(coverage(tmp))
      tmp_coverage <- subsetByIntersect(subject = tmp_coverage, query = roi)
      tmp_coverage$score <- tmp_coverage$score / (total_algns * algns_le_250) * 1e6  # convert to RPM
    
      # save
      file_name <- paste0(dash_to_underscore(sample), "_MNase_seq_", max_insert, "bp")
      cat("\nSaving as", file_name, "...")
      assign(x = file_name, value = tmp_coverage)
    }
  }
}

cat("\n")

# save to file
for(max_insert in c(147, 127, 117, 107, 97, 87, 74)){
  file_list <- paste0(rep(c("LSY4518_13B", "LSY5415", "LSY5934", "LSY5935"), rep(4, 4)), "_", c(0, 1, 2, 4), "h_MNase_seq_", max_insert, "bp")
  file_name <- paste0(save_dir, "/MNase-seq_", max_insert, "bp.RData")
  cat("\nSaving coverage data to", file_name, "...")
  save(list = file_list, file = file_name)
}
