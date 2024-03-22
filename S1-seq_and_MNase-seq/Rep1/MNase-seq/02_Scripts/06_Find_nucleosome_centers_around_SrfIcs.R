# info --------------------------------------------------------------------
# purpose: identify nucleosome positions next to DSBs
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/22/24
# last modified: 03/22/24

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


# file base paths ---------------------------------------------------------
BAM_dir <- "/media/robert/Elements/Deep_sequencing_data/24-03-20-MNase-seq/BAM"
save_dir <- "03_Processed_data/Nucleosome_positions"
dir.create(path = save_dir, showWarnings = FALSE)


# define DSB regions ------------------------------------------------------
DSBs <- SrfIcs[-c(9, 17)]  # exclude duplicated genome regions
roi <- DSB_regions(DSBs = DSBs, region_width = 6000)

# process all samples =====================================================

for(strain in c("LSY4518-13B", "LSY5415", "LSY5934", "LSY5935")){
  
  plot_dir <- paste0("04_Plots/Nucleosome_centers/", strain)
  dir.create(path = plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for(t in c(0, 1, 2, 4)){
    
    samples <- paste0(strain, "_", t, "h")
    
    for(sample in samples){
      
      cat("\n\nProcessing ", sample, "...", sep = "")
      
      cat("\nReading and processing BAM file...")

      # read alignments in DSB regions
      params <- ScanBamParam(which = roi, what = scanBamWhat())  # only read alignments in region of interest
      tmp <- readGAlignmentPairs(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = params)
      total_algns <- sum(idxstatsBam(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"))$mapped)  # count total mapped alignments for normalization
      
      # process alignments
      tmp <- GRanges(tmp)  # convert to GRanges object (necessary for size manipulations)
      tmp <- tmp[width(tmp) <= 250]  # remove alignments larger 250 bp
      tmp_untrimmed <- tmp  # for plotting
      tmp[width(tmp) > 40] <- resize(x = tmp[width(tmp) > 40], width = 40, fix = "center")  # trim to 40 bp size
      
      tmp <- GRanges(coverage(tmp))
      tmp$score <- tmp$score / total_algns * 1e6  # convert to RPM
      tmp <- as_nt_resolved_GRanges(GRanges = tmp)  # make nt-resolved
      tmp <- subsetByIntersect(subject = tmp, query = roi)
      
      tmp_untrimmed <- GRanges(coverage(tmp_untrimmed))  # for plotting
      tmp_untrimmed$score <- tmp_untrimmed$score / total_algns * 1e6  # convert to RPM
      tmp_untrimmed <- as_nt_resolved_GRanges(GRanges = tmp_untrimmed)
      tmp_untrimmed <- subsetByIntersect(subject = tmp_untrimmed, query = roi)
      
      nuc_pos <- GRanges()  # to collect identified nucleosome centers
      
      for(r in 1:length(roi)){
        
        DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(SrfIcs[r]))
        cat("\nIdentifying nucleosome positions at", DSB_locus, "...")
        
        tmp_roi <- subsetByIntersect(subject = tmp, query = roi[r])
        tmp_untrimmed_roi <- subsetByIntersect(subject = tmp_untrimmed, query = roi[r])
        
        # filter noise using Fast Fourier Transform
        fft <- filterFFT(data = tmp_roi$score, pcKeepComp = 0.01) 
        # detect peaks
        peaks <- peakDetection(fft, threshold="25%", score=TRUE, min.cov = 0.1, width = 147)
        # width = 147 is necessary to get also the width score ("fuzziness")
        
        peaks <- resize(x = peaks, width = 1, fix = "center")  # to use start(peaks) as indexes on tmp_roi
        colnames(mcols(peaks)) <- c("nuc_score", "nuc_score_w", "nuc_score_h")
        
        out <- granges(subsetByIntersect(subject = tmp_untrimmed, query = roi[r])[start(peaks)])
        mcols(out) <- mcols(peaks)
        
        nuc_pos <- c(nuc_pos, out)
        
        # plotting
        tmp_untrimmed_roi$score <- tmp_untrimmed_roi$score / max(tmp_untrimmed_roi$score)
        fft_Gpos <- tmp_untrimmed_roi
        fft_Gpos$score <- fft / max(fft)
        nuc_Gpos <- out
        nuc_Gpos$score <- fft_Gpos$score[start(peaks)]
        mcols(nuc_Gpos) <- nuc_Gpos$score
        colnames(mcols(nuc_Gpos)) <- "score"
        
        DT_untrimmed <- DataTrack(range = tmp_untrimmed_roi, type = "polygon", col = NA, fill.mountain = rep("gray", 2), name = "MNase-seq", ylim = c(0, 1))
        DT_fft <- DataTrack(range = fft_Gpos, type = "l", col = JFly_colors[1], name = "MNase-seq", ylim = c(0, 1))
        DT_nuc_pos <- DataTrack(range = nuc_Gpos, type = "h", col = JFly_colors[8], name = "MNase-seq", ylim = c(0, 1))
        OT <- OverlayTrack(trackList = list(DT_untrimmed, DT_fft, DT_nuc_pos))
        
        pdf(file = "tmp.pdf", width = 3, height = 1.5)
          plotTracks(trackList = list(OT), from = start(roi[r]), to = end(roi[r]), chromosome = seqnames(roi[r]), showTitle = FALSE, showAxis = FALSE, margin = 0)
        dev.off()
        GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "_", t ,"h.pdf"))
        
      }
      
      assign(x = paste0(dash_to_underscore(sample), "_nucleosome_positions"), value = nuc_pos)
      
    }
    
  }
  
  # save data
  cat("\n\nSaving data...")
  save(list = paste0(dash_to_underscore(samples), "_nucleosome_positions"), file = paste0(save_dir, "/", strain, "_nucleosome_positions.RData"))
  
}

