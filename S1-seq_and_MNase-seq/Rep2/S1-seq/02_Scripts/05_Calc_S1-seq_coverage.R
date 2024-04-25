# info --------------------------------------------------------------------
# purpose: filter alignments, calculate S1-seq coverage, and extract mapping stats
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 04/23/24
# last modified: 04/23/24

# load libraries ----------------------------------------------------------
library(GenomicAlignments)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")


# function definitions ====================================================

# remove alignments with 5' soft-clipping ---------------------------------
# argument: GAlignments object
# result: GAlignments object without 5' soft-clipped alignments
remove_5prime_S <- function(GAlignments, return_filtered_out_GAlignments = FALSE, verbose = FALSE){
  if(verbose){
    cat("\nExtracting CIGAR ops...")
  }
  cigar_ops_list <- explodeCigarOps(cigar = cigar(GAlignments))
  if(verbose){
    cat("\nFinding alignments with 5' soft clipping...")
  }
  has_5prime_S <- unlist(lapply(X = cigar_ops_list, FUN = function(x) {x[1] == "S"}))
  if(verbose){
    cat("\nFiltering out alignments with 5' soft clipping...")
  }
  if(return_filtered_out_GAlignments){
    out <- list(remaining = GAlignments[!has_5prime_S], filtered_out = GAlignments[has_5prime_S])
  }else{
    out <- GAlignments[!has_5prime_S]  
  }
  return(out)
}


# calculate strand-specific coverage --------------------------------------
# argument: GAlignments object
# result: GRanges object with coverage (mcol "score") for each strand
stranded_coverage <- function(GAlignments, rpm_normalization = TRUE, verbose = FALSE){
  if(verbose){
    cat("\nSplitting alignments according to strand...")
  }
  GAlignments_by_strand <- split(x = GAlignments, f = strand(GAlignments), drop = FALSE)  # keeps "+" and "-" levels, but drops unused "*" level
  if(verbose){
    cat("\nCalculating coverage per strand...")
  }
  fw <- GRanges(coverage(x = GAlignments_by_strand$`+`))
  rev <- GRanges(coverage(x = GAlignments_by_strand$`-`))
  strand(fw) <- "+"
  strand(rev) <- "-"
  if(verbose){
    cat("\nConcatenating and sorting...")
  }
  out <- sort(x = c(fw, rev), ignore.strand = TRUE)
  if(rpm_normalization){
    out$score <- out$score / length(GAlignments) * 1e6
  }
  return(out)
}


# process all samples =====================================================

# file base paths
BAM_dir <- "03_Processed_data/BAM"
# BAM_dir <- "/media/robert/Elements/Deep_sequencing_data/24-03-13-S1-seq/BAM"
save_dir <- "03_Processed_data/S1-seq_coverage"
dir.create(path = save_dir, showWarnings = FALSE)
plot_dir <- "04_Plots/Mapping_edit_distance_distributions"
dir.create(path = plot_dir, showWarnings = FALSE)

# initialize data.frames for result recording
mapping_stats <- data.frame()

# iterate through samples
samples <- paste0(rep(c("LSY4518-13B", "LSY5415"), rep(4, 2)), "_", c(0, 1, 2, 4), "h")

for(sample in samples){
  
  cat("\n\nProcessing ", sample, "...", sep = "")
  
  # read BAM file -----------------------------------------------------------
  cat("\nReading BAM file...")
  tmp <- readGAlignments(file = paste0(BAM_dir, "/", sample, "/", sample, ".bam"), param = ScanBamParam(tag = c("AS", "XS", "NM")))
  # AS: alignment score (max. read length * match bonus [--ma, default: 2])
  # XS: alignment score for 2nd best alignment
  # NM: edit distance
  mapped <- length(tmp)  # for mapping statistics
  cat(mapped, "alignments read.")

  # remove alignments with 5' soft clipping ---------------------------------
  cat("\nRemoving alignments with 5' soft clipping...")
  tmp <- remove_5prime_S(GAlignments = tmp)
  alignments_wo_5prime_S <- length(tmp)
  cat(" kept ", alignments_wo_5prime_S, " alignments (",  round(100 * alignments_wo_5prime_S / mapped), "%).", sep = "")
  
  # remove non-unique alignments --------------------------------------------
  cat("\nRemoving non-unique alignments...")
  idx <- which(!is.na(mcols(tmp)$XS))  # unique alignments have XS == NA
  idx <- idx[mcols(tmp[idx])$XS == mcols(tmp[idx])$AS]  # find non-unique alignments (where primary and alternative alignment have same alignment score; XS == AS)
  tmp <- tmp[-idx]
  unique_alignments <- length(tmp)
  cat(" kept ", unique_alignments, " alignments (",  round(100 * unique_alignments / mapped), "%).", sep = "")
  
  # plot edit distance distribution -----------------------------------------
  cat("\nPlotting edit distance distribution...")
  pdf(file = "tmp.pdf", width=3, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(1.6, 0.6, 3.6, 2.1), las = 1, tcl = -0.3, mgp = c(2.5, 0.6, 0))
    # plot fractions of alignments with edit distance = 0, 1, 2, ..., 5, >5
    h <- hist(x = mcols(tmp)$NM, breaks = 0:max(mcols(tmp)$NM), plot = FALSE)  # use hist function to calculate densities (fractions)
    bp <- barplot(height = c(h$density[1:6], sum(h$density[7:length(h$density)])), ylim = c(0, 1),
                  names.arg = c(as.character(0:5), ">5"), ylab = "Fraction of Alignments", xlab = "Edit Distance")  # plot density values as bar plot
    # add edit distance <= 3 threshold
    th <- mean(bp[4:5, ])
    abline(v = th, col = "red")
    text(x = 0.95 * th, y = 1, adj = c(1,1), xpd = TRUE, col = "red", labels = paste0(round(x = 100 * ecdf(mcols(tmp)$NM)(3), digits = 2), "%"))
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", sample, ".pdf"))
  
  # calculate S1-seq coverage ----------------------------------------------- 
  # save RPM normalized coverage
  # also save un-normalized coverage e.g. for sequence impact analysis
  cat("\nCalculating strand-specific S1-seq coverage...")
  tmp <- resize(x = GRanges(tmp), width = 1)  # only keep 5' ends
  tmp_coverage <- stranded_coverage(GAlignments = tmp, rpm_normalization = TRUE)
  assign(x = paste0(dash_to_underscore(sample), "_S1_seq"), value = tmp_coverage)
  tmp_coverage <- stranded_coverage(GAlignments = tmp, rpm_normalization = FALSE)
  assign(x = paste0(dash_to_underscore(sample), "_S1_seq_unnormalized"), value = tmp_coverage)

  # record mapping statistics data ------------------------------------------
  mapping_stats <- rbind(mapping_stats,
                         data.frame(sample = sample,
                                    reads_mapped = mapped,
                                    wo_5primeS = alignments_wo_5prime_S,
                                    unique = unique_alignments))

}

# save data ---------------------------------------------------------------
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY4518-13B", x = samples)]), "_S1_seq"), file = paste0(save_dir, "/LSY4518-13B_S1-seq.RData"))
save(list = paste0(dash_to_underscore(samples[grepl(pattern = "LSY5415", x = samples)]), "_S1_seq"), file = paste0(save_dir, "/LSY5415_S1-seq.RData"))
save(list = paste0(dash_to_underscore(samples), "_S1_seq_unnormalized"), file = paste0(save_dir, "/S1-seq_unnormalized.RData"))

write.table(x = mapping_stats, file = paste0(plot_dir, "/Mapping_stats.txt"), row.names = FALSE)
