# info --------------------------------------------------------------------
# purpose: test which trimming results in nucleosome peaks of desired sharpness
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


# plotting ================================================================

# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(Gviz)
options(ucscChromosomeNames=FALSE)  # for using custom chromosome names (e.g. "micron")

# read helper functions and files -----------------------------------------
source(file = "../../../Src/Gviz_functions.R")

# load MNase-seq data -----------------------------------------------------
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_147bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_127bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_117bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_107bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_97bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_87bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_74bp.RData")

# use also old data for comparison
load(file = "/home/robert/Research/Manuscripts/My_manuscripts/20-04-17-MRX_nicking_manuscript/Data/MNase-seq/LSY4377-12B_4377-15A_merged/LSY4377-12B_LSY4377-15A_merged_trimmed.RData")

# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/MNase-seq_coverage_after_trimming"
dir.create(path = plot_dir, showWarnings = FALSE, recursive = TRUE)

# plotting colors
MNase_seq_col_1 <- gray(level = 0.6)
MNase_seq_border_col_1 <- NA

# set global Gviz parameters for plotting
options(Gviz.scheme="default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$rotation.title <- 90
scheme$GdObject$fontcolor.title <- "black"
scheme$GdObject$background.title <- "white"
scheme$GdObject$cex.title <- 0.75
scheme$GdObject$col.axis <- "black"
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new scheme
# scheme$GdObject$showAxis <- FALSE
# scheme$GdObject$showTitle <- FALSE

# MNase-seq DataTracks ----------------------------------------------------
DT_147 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_147bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "147 bp")
DT_127 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_127bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "127 bp")
DT_117 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_117bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "117 bp")
DT_107 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_107bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "107 bp")
DT_97 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_97bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "97 bp")
DT_87 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_87bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "87 bp")
DT_74 <- DataTrack(range = LSY4518_13B_0h_MNase_seq_74bp, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "74 bp")
DT_old <- DataTrack(range = LSY4377_12B_0_merged_MNase_seq, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "Old")

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and modified
load(file = "../../../Src/S_cerevisiae_genome_features.RData")
# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT_genes <- AnnotationTrack(range = all_features, name = NULL, showFeatureId = TRUE, cex = 0.67, featureAnnotation = "id",
                            arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "gray")

# plot around all SrfIcs --------------------------------------------------
for (n in 1:length(SrfIcs)){
  
  # define name for data and figure saving
  DSB_locus <- SrfIcs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  roi <- DSB_regions(DSBs = SrfIcs[n], region_width = 4000)
  
  # plot
  pdf(file = "tmp.pdf", width = 3, height = 5)
  plotTracks(trackList = list(AT_genes, DT_old, DT_147, DT_127, DT_117, DT_107, DT_97, DT_87, DT_74), sizes = c(0.05, rep(0.95/8, 8)),
             from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 0)
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, ".pdf"))
  
}
