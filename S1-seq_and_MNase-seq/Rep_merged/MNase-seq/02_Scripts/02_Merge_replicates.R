# info --------------------------------------------------------------------
# purpose: merge replicates
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/06/24
# last modified: 05/06/24

# load libraries ----------------------------------------------------------
library(GenomicRanges)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")

rep1_dir <- "../../Rep1/MNase-seq/03_Processed_data/MNase-seq_coverage/"
rep2_dir <- "../../Rep2/MNase-seq/03_Processed_data/MNase-seq_coverage/"
out_dir <- "03_Processed_data/MNase-seq_coverage/"
dir.create(path = out_dir, showWarnings = FALSE)


# function definitions ----------------------------------------------------

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
  cat("\nDone.\n")
  return(out)
}


# merge all MNase-seq replicates ==========================================
# Samples have the same name for both Rep1 and Rep2.
# So, we need to save them under individual names.

# LSY4518-13B -------------------------------------------------------------
load(file = paste0(rep1_dir, "LSY4518-13B_MNase-seq.RData"))
rep1_0 <- LSY4518_13B_0h_MNase_seq
rep1_1 <- LSY4518_13B_1h_MNase_seq
rep1_2 <- LSY4518_13B_2h_MNase_seq
rep1_4 <- LSY4518_13B_4h_MNase_seq

load(file = paste0(rep2_dir, "LSY4518-13B_MNase-seq.RData"))
rep2_0 <- LSY4518_13B_0h_MNase_seq
rep2_1 <- LSY4518_13B_1h_MNase_seq
rep2_2 <- LSY4518_13B_2h_MNase_seq
rep2_4 <- LSY4518_13B_4h_MNase_seq

LSY4518_13B_0h_MNase_seq <- average_GRanges(rep1 = rep1_0, rep2 = rep2_0, verbose = TRUE)
LSY4518_13B_1h_MNase_seq <- average_GRanges(rep1 = rep1_1, rep2 = rep2_1, verbose = TRUE)
LSY4518_13B_2h_MNase_seq <- average_GRanges(rep1 = rep1_2, rep2 = rep2_2, verbose = TRUE)
LSY4518_13B_4h_MNase_seq <- average_GRanges(rep1 = rep1_4, rep2 = rep2_4, verbose = TRUE)

save(list = paste0("LSY4518_13B_", c(0, 1, 2, 4), "h_MNase_seq"), file = paste0(out_dir, "LSY4518-13B_MNase-seq.RData"))


# LSY5415 -------------------------------------------------------------
load(file = paste0(rep1_dir, "LSY5415_MNase-seq.RData"))
rep1_0 <- LSY5415_0h_MNase_seq
rep1_1 <- LSY5415_1h_MNase_seq
rep1_2 <- LSY5415_2h_MNase_seq
rep1_4 <- LSY5415_4h_MNase_seq

load(file = paste0(rep2_dir, "LSY5415_MNase-seq.RData"))
rep2_0 <- LSY5415_0h_MNase_seq
rep2_1 <- LSY5415_1h_MNase_seq
rep2_2 <- LSY5415_2h_MNase_seq
rep2_4 <- LSY5415_4h_MNase_seq

LSY5415_0h_MNase_seq <- average_GRanges(rep1 = rep1_0, rep2 = rep2_0, verbose = TRUE)
LSY5415_1h_MNase_seq <- average_GRanges(rep1 = rep1_1, rep2 = rep2_1, verbose = TRUE)
LSY5415_2h_MNase_seq <- average_GRanges(rep1 = rep1_2, rep2 = rep2_2, verbose = TRUE)
LSY5415_4h_MNase_seq <- average_GRanges(rep1 = rep1_4, rep2 = rep2_4, verbose = TRUE)

save(list = paste0("LSY5415_", c(0, 1, 2, 4), "h_MNase_seq"), file = paste0(out_dir, "LSY5415_MNase-seq.RData"))
