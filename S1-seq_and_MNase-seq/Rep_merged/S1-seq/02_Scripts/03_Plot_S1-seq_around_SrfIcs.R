# info --------------------------------------------------------------------
# purpose: generate S1-seq coverage plots around genomic SrfIcs in compact form
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 06/01/25
# last modified: 06/01/25

# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(Gviz)
options(ucscChromosomeNames=FALSE)  # for using custom chromosome names (e.g. "micron")

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")
source(file = "../../../Src/Gviz_functions.R")

# load and process S1-seq data --------------------------------------------
SrfIcs <- SrfIcs[-c(9, 17)]  # exclude SrfIcs in duplicated regions
roi <- DSB_regions(DSBs = SrfIcs, region_width = 4000)
# restrict data processing to regions around DSBs
# chose a bigger region size than used for plotting to allow smoothing

process_S1_seq <- function(GRanges, roi, negative_score = FALSE){
  tmp <- subsetByIntersect(subject = GRanges, query = roi)  # only keep S1-seq coverage in DSB regions
  tmp <- sort(as_nt_resolved_GRanges(tmp), ignore.strand = TRUE)  # nt resolution, and sort
  tmp <- GRanges_zero_to_NA(tmp)
  # tmp <- as_smoothed_GRanges(GRanges = tmp, hanning_window_size = 51)  # Hanning smoothing
  if(negative_score){
    tmp$score <- -tmp$score
  }
  return(tmp)
}

load(file = "03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_0h_S1_seq, roi = roi)
LSY4518_13B_1h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_1h_S1_seq, roi = roi)
LSY4518_13B_2h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_2h_S1_seq, roi = roi)
LSY4518_13B_4h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_4h_S1_seq, roi = roi)

load(file = "03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h_S1_seq <- process_S1_seq(GRanges = LSY5415_0h_S1_seq, roi = roi, negative_score = TRUE)
LSY5415_1h_S1_seq <- process_S1_seq(GRanges = LSY5415_1h_S1_seq, roi = roi, negative_score = TRUE)
LSY5415_2h_S1_seq <- process_S1_seq(GRanges = LSY5415_2h_S1_seq, roi = roi, negative_score = TRUE)
LSY5415_4h_S1_seq <- process_S1_seq(GRanges = LSY5415_4h_S1_seq, roi = roi, negative_score = TRUE)
# make scores negative for plotting below WT data

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and modified
load(file = "../../../Src/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT_genes <- AnnotationTrack(range = all_features, name = NULL, showFeatureId = FALSE, cex = 0.67, featureAnnotation = "id",
                            arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "gray")

# # to add SrfIcs location to Annotation track
# AT_SrfIcs <- AnnotationTrack(range = SrfIcs, name = NULL, showFeatureID = FALSE,
#                              shape = "arrow", fill = "red", col = NA, min.height = 0.01, min.width = 0.01)
# # merge gene and SrfIcs annotation tracks
# OT_AT <- OverlayTrack(trackList = list(AT_genes, AT_SrfIcs))
# # change "AT_genes" to "OT_AT" in plotTracks command below!

# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/S1-seq_Plots"
dir.create(path = plot_dir, showWarnings = FALSE, recursive = TRUE)

# set global Gviz parameters for plotting
options(Gviz.scheme="default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$showAxis <- FALSE
scheme$GdObject$showTitle <- FALSE
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new scheme

# plotting colors
fw_col_1 <- JFly_colors[2]
rev_col_1 <- JFly_colors[7]
fw_col_2 <- JFly_colors[8]
rev_col_2 <- JFly_colors[3]

# to calc nice ylim for data
calc_ylim <- function(GRanges, roi, symmetric = FALSE){
  tmp <- range(subsetByIntersect(subject = GRanges, query = roi)$score, na.rm = TRUE)
  if(symmetric){
    tmp <- max(abs(tmp))
    tmp <- c(-tmp, tmp)
  }
  return(tmp)
}

# iterate through all SrfIcs
n <- 16
for (n in 1:length(SrfIcs)){
  
  # define name for data and figure saving
  DSB_locus <- SrfIcs[n]
  DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
  cat("\nPlotting", DSB_locus, "...")
  
  roi <- DSB_regions(DSBs = SrfIcs[n], region_width = 1500)
  
  # calculate S1-seq ylim
  all_S1_seq <- c(LSY4518_13B_0h_S1_seq, LSY4518_13B_1h_S1_seq, LSY4518_13B_2h_S1_seq, LSY4518_13B_4h_S1_seq,
                  LSY5415_0h_S1_seq, LSY5415_1h_S1_seq, LSY5415_2h_S1_seq, LSY5415_4h_S1_seq)
  ylim_S1_seq <- calc_ylim(GRanges = all_S1_seq, roi = roi, symmetric = TRUE)
  
  # generate S1-seq DataTracks and merge WT and fun30 into OverlayTracks
  DT_S1_seq_WT_0h <- DataTrack_fw_rev(GRanges = LSY4518_13B_0h_S1_seq, fw_col = fw_col_1, rev_col = rev_col_1, type = "h", name = "0 h", ylim = ylim_S1_seq)
  DT_S1_seq_WT_1h <- DataTrack_fw_rev(GRanges = LSY4518_13B_1h_S1_seq, fw_col = fw_col_1, rev_col = rev_col_1, type = "h", name = "1 h", ylim = ylim_S1_seq)
  DT_S1_seq_WT_2h <- DataTrack_fw_rev(GRanges = LSY4518_13B_2h_S1_seq, fw_col = fw_col_1, rev_col = rev_col_1, type = "h", name = "2 h", ylim = ylim_S1_seq)
  DT_S1_seq_WT_4h <- DataTrack_fw_rev(GRanges = LSY4518_13B_4h_S1_seq, fw_col = fw_col_1, rev_col = rev_col_1, type = "h", name = "4 h", ylim = ylim_S1_seq)
  
  DT_S1_seq_fun30_0h <- DataTrack_fw_rev(GRanges = LSY5415_0h_S1_seq, fw_col = fw_col_2, rev_col = rev_col_2, type = "h", name = "0 h", ylim = ylim_S1_seq)
  DT_S1_seq_fun30_1h <- DataTrack_fw_rev(GRanges = LSY5415_1h_S1_seq, fw_col = fw_col_2, rev_col = rev_col_2, type = "h", name = "1 h", ylim = ylim_S1_seq)
  DT_S1_seq_fun30_2h <- DataTrack_fw_rev(GRanges = LSY5415_2h_S1_seq, fw_col = fw_col_2, rev_col = rev_col_2, type = "h", name = "2 h", ylim = ylim_S1_seq)
  DT_S1_seq_fun30_4h <- DataTrack_fw_rev(GRanges = LSY5415_4h_S1_seq, fw_col = fw_col_2, rev_col = rev_col_2, type = "h", name = "4 h", ylim = ylim_S1_seq)
  
  OT_S1_seq_0h <- OverlayTrack(trackList = list(DT_S1_seq_WT_0h, DT_S1_seq_fun30_0h))
  OT_S1_seq_1h <- OverlayTrack(trackList = list(DT_S1_seq_WT_1h, DT_S1_seq_fun30_1h))
  OT_S1_seq_2h <- OverlayTrack(trackList = list(DT_S1_seq_WT_2h, DT_S1_seq_fun30_2h))
  OT_S1_seq_4h <- OverlayTrack(trackList = list(DT_S1_seq_WT_4h, DT_S1_seq_fun30_4h))
  
  # plot
  dir.create(path = paste0(plot_dir, "/", DSB_locus), showWarnings = FALSE)
  pdf(file = "tmp.pdf", width = 3, height = 4)
  plotTracks(trackList = list(OT_S1_seq_0h, OT_S1_seq_1h,  OT_S1_seq_2h, OT_S1_seq_4h, AT_genes), sizes = c(rep(0.93/4, 4), 0.07),
             from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 0, showTitle = FALSE)
  # scheme$GdObject$showAxis <- FALSE was defined above. 
  # But without "showTitle = FALSE", a gray title area is plotted anyways. Bug?
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "/Plot.pdf"))
  
  # save axis ranges
  tmp <- data.frame(SrfIcs = DSB_locus, region_size = width(roi),
                    S1_seq_min = ylim_S1_seq[1], S1_seq_max = ylim_S1_seq[2])
  write.table(x = tmp, file = paste0(plot_dir, "/", DSB_locus, "/axis_ranges.txt"), row.names = FALSE)

}

subsetByIntersect(subject = all_features, query = roi)
