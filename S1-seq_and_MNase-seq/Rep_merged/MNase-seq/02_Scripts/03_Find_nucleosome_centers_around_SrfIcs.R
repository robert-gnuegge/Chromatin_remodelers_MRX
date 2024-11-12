# info --------------------------------------------------------------------
# purpose: identify nucleosome positions next to DSBs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/26/24
# last modified: 11/11/24

# Comment: For nucleosome peak identification, we need to artificially shorten the 
# alignments more aggressively (for sharper peaks). Therefore, we cannot use the merged 
# MNase-seq data, where the alingments where shorten to a longer length.
# Thus, we read the original BAM data in the DSB regions, shorten the alignments,
# and average.

# load libraries ----------------------------------------------------------
library(GenomicAlignments)
library(nucleR)
library(Gviz)
options(ucscChromosomeNames=FALSE)  # for using custom chromosome names (e.g. "micron")

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")


# define functions --------------------------------------------------------

# process paired-end alignments - remove disomes, trim to defined size, calculate coverage
# arguments: GAlignmentPairs, numeric, numeric
# result: GRanges
process_alnmts <- function(GAlignmentPairs, max_width, total_algns){
  tmp <- GRanges(GAlignmentPairs)  # convert to GRanges object (necessary for size manipulations)
  tmp <- tmp[width(tmp) <= 250]  # remove alignments larger 250 bp (probably disomes)
  tmp[width(tmp) > max_width] <- resize(x = tmp[width(tmp) > max_width], width = max_width, fix = "center")
  tmp <- GRanges(coverage(tmp))
  tmp$score <- tmp$score / total_algns * 1e6  # convert to RPM
  tmp <- as_nt_resolved_GRanges(GRanges = tmp)
  subsetByIntersect(subject = tmp, query = roi)
}

# average two GRanges
# argument: GRanges
# result: GRanges
# note: only GRanges with numeric mcols (and identical mcol names) are supported
average_GRanges <- function(rep1, rep2, verbose = FALSE){
  stopifnot(names(mcols(rep1)) %in% names(mcols(rep2)))
  if(verbose){
    cat("\nMaking sorted nt-resolved GRanges...")
  }
  rep1 <- sort(as_nt_resolved_GRanges(GRanges = rep1))
  rep2 <- sort(as_nt_resolved_GRanges(GRanges = rep2))
  if(verbose){
    cat("\nChecking identity of all granges of replicates...")
  }
  stopifnot(all(granges(rep1) == granges(rep2)))
  if(verbose){
    cat("\nSetting up output GRanges object...")
  }
  out <- granges(rep1)
  mcol_names <- names(mcols(rep1))
  mcols(out) <- matrix(data = rep(NA, length(out) * length(mcol_names)), ncol = length(mcol_names))
  names(mcols(out)) <- mcol_names
  for(n in 1:length(mcol_names)){
    if(verbose){
      cat("\nAveraging mcol '", mcol_names[n], "'...", sep = "")
    }
    mcol1 <- mcols(rep1)[names(mcols(rep1)) == mcol_names[n]]
    mcol2 <- mcols(rep2)[names(mcols(rep2)) == mcol_names[n]]
    mcols(out)[names(mcols(out)) == mcol_names[n]] <- apply(X = cbind(mcol1$score, mcol2$score), MARGIN = 1, FUN = mean)
  }
  if(verbose){
    cat("\nDone.\n")
  }
  return(out)
}

# plotting function
# arguments: GPos objects and color vector
# result: GRanges
my_plot <- function(raw_GPos, fft_GPos, nuc_GPos, DSB_GPos, colors = c(gray(level = c(0.75, 0.5)), JFly_colors[c(1, 8)])){
  
  x_raw <- start(raw_GPos)
  y_raw <- raw_GPos$score
  x_fft <- start(fft_GPos)
  y_fft <- fft_GPos$score
  
  # start empty plot
  plot(x = NA, y = NA, xlim = range(c(x_raw, x_fft)), ylim = range(c(y_raw, y_fft)), 
       axes = FALSE, ann = FALSE)
  
  # plot raw MNase-seq coverage
  polygon(x = c(x_raw[1], x_raw, x_raw[length(x_raw)]), 
          y = c(0, y_raw, 0), col = colors[1], border = NA)
  
  # plot fft MNase-seq coverage
  points(x = x_fft, y = y_fft, type = "l", col = colors[2])
  
  # plot DSB location
  abline(v = start(DSB_GPos), col = colors[3], lty = "dotted", xpd = TRUE)
  
  # plot nucleosome locations
  segments(x0 = start(nuc_GPos), y0 = 0, x1 = start(nuc_GPos), y1 = fft_GPos$score[start(fft_GPos) %in% start(nuc_GPos)], col = colors[4], lty = "dashed")
  text(x = start(nuc_GPos), y = par("usr")[3], labels = nuc_GPos$idx, pos = 1, offset = 0.06, cex = 0.7, col = colors[4], xpd = TRUE)
  
}

# file base paths ---------------------------------------------------------
BAM_dir_Rep1 <- "/media/robert/Elements/Deep_sequencing_data/24-03-20-MNase-seq/BAM"
BAM_dir_Rep2 <- "/media/robert/Elements/Deep_sequencing_data/24-04-24-MNase-seq/BAM"
save_dir <- "03_Processed_data/Nucleosome_positions"
dir.create(path = save_dir, showWarnings = FALSE)

# define DSB regions ------------------------------------------------------
DSBs <- SrfIcs[-c(9, 17)]  # exclude duplicated genome regions
roi <- DSB_regions(DSBs = DSBs, region_width = 3000)


# iterate through all strains
for(strain in c("LSY4518-13B", "LSY5415", "LSY5934", "LSY5935")){
  
  cat("\n\nProcessing ", strain, "...", sep = "")  

  # read and process BAM files ----------------------------------------------
  cat("\n\nReading and processing all BAM files...")
  
  for(t in c(0, 1, 2, 4)){
    
    cat("\n\nProcessing t = ", t, " h...", sep = "")
    sample <- paste0(strain, "_", t, "h")

    cat("\nReading and processing Rep1 BAM file...")
    params <- ScanBamParam(which = roi, what = scanBamWhat())  # only read alignments in roi
    tmp <- readGAlignmentPairs(file = paste0(BAM_dir_Rep1, "/", sample, "/", sample, ".bam"), param = params)
    total_algns <- sum(idxstatsBam(file = paste0(BAM_dir_Rep1, "/", sample, "/", sample, ".bam"))$mapped)  # count total mapped alignments for normalization
    # trim and calculate coverage
    coverage_117_Rep1 <- process_alnmts(GAlignmentPairs = tmp, max_width = 117, total_algns = total_algns)
    coverage_40_Rep1 <- process_alnmts(GAlignmentPairs = tmp, max_width = 40, total_algns = total_algns)
    
    cat("\nReading and processing Rep2 BAM file...")
    params <- ScanBamParam(which = roi, what = scanBamWhat())  # only read alignments in roi
    tmp <- readGAlignmentPairs(file = paste0(BAM_dir_Rep2, "/", sample, "/", sample, ".bam"), param = params)
    total_algns <- sum(idxstatsBam(file = paste0(BAM_dir_Rep2, "/", sample, "/", sample, ".bam"))$mapped)  # count total mapped alignments for normalization
    # trim and calculate coverage
    coverage_117_Rep2 <- process_alnmts(GAlignmentPairs = tmp, max_width = 117, total_algns = total_algns)
    coverage_40_Rep2 <- process_alnmts(GAlignmentPairs = tmp, max_width = 40, total_algns = total_algns)
    
    cat("\nAveraging replicates...")
    assign(x = paste0("coverage_117_", t, "h"), value = average_GRanges(rep1 = coverage_117_Rep1, rep2 = coverage_117_Rep2))
    assign(x = paste0("coverage_40_", t, "h"), value = average_GRanges(rep1 = coverage_40_Rep1, rep2 = coverage_40_Rep2))
    
  }
  
  # identify nucleosome positions -------------------------------------------
  
  # set up GRanges objects to collect identified nucleosome centers
  nuc_pos_0h <- GRanges()
  nuc_pos_1h <- GRanges()
  nuc_pos_2h <- GRanges()
  nuc_pos_4h <- GRanges()
  
  # set up plotting
  plot_dir <- paste0("04_Plots/Nucleosome_centers/", strain)
  dir.create(path = plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for(r in 1:length(roi)){
    
    cat("\n\nIdentifying nucleosome positions at", as.character(SrfIcs[r]), "...")
    
    # start PDF device for plotting all time points
    pdf(file = "tmp.pdf", width=3, height=5)
    par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(4.5, 4.1, 4.1, 2.1), oma = rep(0, 4), mfrow = c(4, 1))
    
    for(t in c(0, 1, 2, 4)){
      
      cat("\nIdentifying nucleosome positions for t =", t, "h...")
      
      # find MNase-seq coverage peaks
      coverage_40_roi <- subsetByIntersect(subject = get(paste0("coverage_40_", t, "h")), query = roi[r])
      fft <- filterFFT(data = coverage_40_roi$score, pcKeepComp = 0.01)  # filter noise
      idx <- peakDetection(fft, threshold="25%", score=FALSE, min.cov = 0.1)
      peaks <- granges(coverage_40_roi)[idx]
      
      if(t == 0){
        # identify DSB-proximal nucleosome
        tmp_hits <- distanceToNearest(x = DSBs[r], subject = peaks)  
        DSB_prox_nuc <- peaks[subjectHits(tmp_hits)]
        DSB_prox_nuc_dist <- mcols(tmp_hits)$distance
        
        # specify where neighboring nucleosomes would be expected ("ideal nucleosomes")
        nuc_dist <- 165
        ideal_nucs <- rep(granges(DSB_prox_nuc), 21)
        ranges(ideal_nucs[1:11]) <- rev(start(DSB_prox_nuc) - 0:10 * nuc_dist)
        ranges(ideal_nucs[12:21]) <- start(DSB_prox_nuc) + 1:10 * nuc_dist
        # add numbers
        if(DSB_prox_nuc_dist > 37){
          # if distance of DSB-proximal nucleosome to DSB is > 37 bp (cut at > 1/4 nucleosome width)
          # nucleosomes are labeled ..., -2, -1, 1, 2, ...
          ideal_nucs$idx <- c(-11:-1, 1:10)
        } else {
          # if distance is smaller (DSB within nucleosome)
          # nucleosomes are labeled ..., -2, -1, 0, 1, 2, ...
          # (the DSB-proximal nucleosome is labeled "0")
          ideal_nucs$idx <- -10:10
        }
        
        # The DSB-proximal nucleosome might be a bit out of register.
        # Slide the "ideal nucleosomes" GRanges a bit up and down and find
        # optimal shift to minimize the distance to the "real nucleosomes".
        optimal_shift <- 0
        lowest_summed_dist <- 10000
        for(nt_shift in (-41:41)){
          shifted_ideal_nucs <- shift(x = ideal_nucs, shift = nt_shift)
          tmp_hits <- distanceToNearest(x = peaks, subject = shifted_ideal_nucs)
          summed_dist <- sum(mcols(tmp_hits)$distance)
          if(summed_dist < lowest_summed_dist){
            optimal_shift <- nt_shift
            lowest_summed_dist <- summed_dist
          }
        }
        
        # assign nucleosome numbers to MNase-seq coverage peaks
        shifted_ideal_nucs <- shift(x = ideal_nucs, shift = optimal_shift)
        tmp_hits <- distanceToNearest(x = peaks, subject = shifted_ideal_nucs)
        peaks$idx <- shifted_ideal_nucs$idx[subjectHits(tmp_hits)]
        peaks_0h <- peaks
      
      }else{
        
        # assign nucleosome numbers based on t = 0 nucleosomes 
        tmp_hits <- distanceToNearest(x = peaks, subject = peaks_0h)
        peaks$idx <- peaks_0h$idx[subjectHits(tmp_hits)]
        
      }
      
      # save identified nucleosomes to "nuc_pos_...h"
      assign(x = paste0("nuc_pos_", t, "h"), value = c(get(x = paste0("nuc_pos_", t, "h")), peaks))
  
      # plotting ----------------------------------------------------------------
      coverage_117_roi <- subsetByIntersect(subject = get(paste0("coverage_117_", t, "h")), query = roi[r])
      coverage_117_roi$score <- coverage_117_roi$score / max(coverage_117_roi$score)  # scale to 0...1
      coverage_fft <- coverage_40_roi
      coverage_fft$score <- fft / max(fft)  # scale to 0...1

      my_plot(raw_GPos = coverage_117_roi, fft_GPos = coverage_fft, nuc_GPos = peaks, DSB_GPos = DSBs[r])
    }
    
    dev.off()  # close plotting device
    file_name = paste0(plot_dir, "/", sub(pattern = ":", replacement = "_", x = as.character(SrfIcs[r])), ".pdf")
    GS_embed_fonts(input = "tmp.pdf", output = file_name)  # embed fonts into plot
    par(mfrow = c(1, 1))  # reset mfrow plotting parameter

  }
  
  # save to file
  cat("\n\nSaving data...")
  assign(x = paste0(dash_to_underscore(strain), "_0h_nucleosome_positions"), value = nuc_pos_0h)
  assign(x = paste0(dash_to_underscore(strain), "_1h_nucleosome_positions"), value = nuc_pos_1h)
  assign(x = paste0(dash_to_underscore(strain), "_2h_nucleosome_positions"), value = nuc_pos_2h)
  assign(x = paste0(dash_to_underscore(strain), "_4h_nucleosome_positions"), value = nuc_pos_4h)
  save(list = paste0(dash_to_underscore(strain), "_", c(0, 1, 2, 4), "h_nucleosome_positions"), 
       file = paste0(save_dir, "/", strain, "_nucleosome_positions.RData"))
  
}
