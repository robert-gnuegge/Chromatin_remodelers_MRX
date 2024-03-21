# info --------------------------------------------------------------------
# purpose: generate S1-seq and MNase-seq coverage plots around genomic SrfIcs in compact form
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/21/24
# last modified: 03/21/24

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
roi <- DSB_regions(DSBs = SrfIcs, region_width = 4000)  # restrict data processing to regions around DSBs
# chose a bigger region size than used for plotting to allow smoothing

process_S1_seq <- function(GRanges, roi, negative_score = FALSE){
  tmp <- subsetByIntersect(subject = GRanges, query = roi)  # only keep S1-seq coverage in DSB regions
  tmp <- sort(as_nt_resolved_GRanges(tmp), ignore.strand = TRUE)  # nt resolution, and sort
  tmp <- GRanges_zero_to_NA(tmp)
  # tmp <- as_smoothed_GRanges(GRanges = tmp, hanning_window_size = 51)  # Hanning smoothing
  # !!! if smoothing, use a larger roi first, as the smoothing requires neighboring positions
  if(negative_score){
    tmp$score <- -tmp$score
  }
  return(tmp)
}

load(file = "../S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq.RData")
LSY4518_13B_0h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_0h_S1_seq, roi = roi)
LSY4518_13B_1h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_1h_S1_seq, roi = roi)
LSY4518_13B_2h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_2h_S1_seq, roi = roi)
LSY4518_13B_4h_S1_seq <- process_S1_seq(GRanges = LSY4518_13B_4h_S1_seq, roi = roi)

load(file = "../S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq.RData")
LSY5415_0h_S1_seq <- process_S1_seq(GRanges = LSY5415_0h_S1_seq, roi = roi, negative_score = TRUE)
LSY5415_1h_S1_seq <- process_S1_seq(GRanges = LSY5415_1h_S1_seq, roi = roi, negative_score = TRUE)
LSY5415_2h_S1_seq <- process_S1_seq(GRanges = LSY5415_2h_S1_seq, roi = roi, negative_score = TRUE)
LSY5415_4h_S1_seq <- process_S1_seq(GRanges = LSY5415_4h_S1_seq, roi = roi, negative_score = TRUE)
# make scores negative for plotting below WT data

# load and process MNase-seq data -----------------------------------------
process_MNase_seq <- function(GRanges, roi, negative_score = FALSE){
  tmp <- subsetByIntersect(subject = GRanges, query = roi)  # only keep S1-seq coverage in DSB regions
  # tmp$score <- runmed(x = tmp$score, k = 31)  # smoothing
  if(negative_score){
    tmp$score <- -tmp$score
  }
  return(tmp)
}

load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/LSY4518-13B_MNase-seq_trimmed.RData")
LSY4518_13B_0h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY4518_13B_0h_MNase_seq_trimmed, roi = roi)
LSY4518_13B_1h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY4518_13B_1h_MNase_seq_trimmed, roi = roi)
LSY4518_13B_2h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY4518_13B_2h_MNase_seq_trimmed, roi = roi)
LSY4518_13B_4h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY4518_13B_4h_MNase_seq_trimmed, roi = roi)

load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq_trimmed.RData")
LSY5415_0h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY5415_0h_MNase_seq_trimmed, roi = roi, negative_score = TRUE)
LSY5415_1h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY5415_1h_MNase_seq_trimmed, roi = roi, negative_score = TRUE)
LSY5415_2h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY5415_2h_MNase_seq_trimmed, roi = roi, negative_score = TRUE)
LSY5415_4h_MNase_seq_trimmed <- process_MNase_seq(GRanges = LSY5415_4h_MNase_seq_trimmed, roi = roi, negative_score = TRUE)
# make scores negative for plotting below WT data

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and modified
load(file = "../../../Src/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"
all_features <- all_features[!(all_features$id == "") & all_features$type == "ORF"]

AT_genes <- AnnotationTrack(range = all_features, name = NULL, showFeatureId = TRUE, cex = 0.67, featureAnnotation = "id",
                            arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "gray")

AT_SrfIcs <- AnnotationTrack(range = SrfIcs, name = NULL, showFeatureID = FALSE,
                             shape = "arrow", fill = "red", col = NA, min.height = 0.01, min.width = 0.01)

# merge gene and SrfIcs annotation tracks
OT_AT <- OverlayTrack(trackList = list(AT_genes, AT_SrfIcs))


# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/S1-seq_and_MNase-seq_Plots"
dir.create(path = plot_dir, showWarnings = FALSE, recursive = TRUE)

# set global Gviz parameters for plotting
options(Gviz.scheme="default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$showAxis <- FALSE
scheme$GdObject$showTitle <- FALSE
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new schem

# plotting colors
fw_col_1 <- JFly_colors[2]
rev_col_1 <- JFly_colors[7]
fw_col_2 <- JFly_colors[8]
rev_col_2 <- JFly_colors[3]
MNase_seq_col_1 <- gray(level = 0.6)
MNase_seq_border_col_1 <- NA
MNase_seq_col_2 <- gray(level = 0.8)
MNase_seq_border_col_2 <- NA

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
# Y_axis_ranges <- data.frame()  # initialize for collecting
n <- 1
# for (n in 1:length(SrfIcs)){
  
# define name for data and figure saving
DSB_locus <- SrfIcs[n]
DSB_locus <- gsub(pattern = ":", replacement = "_", x = as.character(DSB_locus))
cat("\nPlotting", DSB_locus, "...")

roi <- DSB_regions(DSBs = SrfIcs[n], region_width = 2000)

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

# calculate MNase-seq ylim
all_MNase_seq <- c(LSY4518_13B_0h_MNase_seq_trimmed, LSY4518_13B_1h_MNase_seq_trimmed, LSY4518_13B_2h_MNase_seq_trimmed, LSY4518_13B_4h_MNase_seq_trimmed,
                   LSY5415_0h_MNase_seq_trimmed, LSY5415_1h_MNase_seq_trimmed, LSY5415_2h_MNase_seq_trimmed, LSY5415_4h_MNase_seq_trimmed)
ylim_MNase_seq <- calc_ylim(GRanges = all_MNase_seq, roi = roi, symmetric = TRUE)

# generate MNase-seq DataTracks and merge WT and fun30 into OverlayTracks
DT_MNase_seq_WT_0h <- DataTrack(range = LSY4518_13B_0h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "0 h", ylim = ylim_MNase_seq)
DT_MNase_seq_WT_1h <- DataTrack(range = LSY4518_13B_1h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "1 h", ylim = ylim_MNase_seq)
DT_MNase_seq_WT_2h <- DataTrack(range = LSY4518_13B_2h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "2 h", ylim = ylim_MNase_seq)
DT_MNase_seq_WT_4h <- DataTrack(range = LSY4518_13B_4h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "4 h", ylim = ylim_MNase_seq)

DT_MNase_seq_fun30_0h <- DataTrack(range = LSY5415_0h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_2, fill.mountain = rep(MNase_seq_col_2, 2), name = "0 h", ylim = ylim_MNase_seq)
DT_MNase_seq_fun30_1h <- DataTrack(range = LSY5415_1h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_2, fill.mountain = rep(MNase_seq_col_2, 2), name = "1 h", ylim = ylim_MNase_seq)
DT_MNase_seq_fun30_2h <- DataTrack(range = LSY5415_2h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_2, fill.mountain = rep(MNase_seq_col_2, 2), name = "2 h", ylim = ylim_MNase_seq)
DT_MNase_seq_fun30_4h <- DataTrack(range = LSY5415_4h_MNase_seq_trimmed, type = "polygon", col = MNase_seq_border_col_2, fill.mountain = rep(MNase_seq_col_2, 2), name = "4 h", ylim = ylim_MNase_seq)

OT_MNase_seq_0h <- OverlayTrack(trackList = list(DT_MNase_seq_WT_0h, DT_MNase_seq_fun30_0h))
OT_MNase_seq_1h <- OverlayTrack(trackList = list(DT_MNase_seq_WT_1h, DT_MNase_seq_fun30_1h))
OT_MNase_seq_2h <- OverlayTrack(trackList = list(DT_MNase_seq_WT_2h, DT_MNase_seq_fun30_2h))
OT_MNase_seq_4h <- OverlayTrack(trackList = list(DT_MNase_seq_WT_4h, DT_MNase_seq_fun30_4h))

# merge S1-seq and MNase-seq OverlayTracks
OT_0h <- OverlayTrack(trackList = list(OT_MNase_seq_0h, OT_S1_seq_0h))
OT_1h <- OverlayTrack(trackList = list(OT_MNase_seq_1h, OT_S1_seq_1h))
OT_2h <- OverlayTrack(trackList = list(OT_MNase_seq_2h, OT_S1_seq_2h))
OT_4h <- OverlayTrack(trackList = list(OT_MNase_seq_4h, OT_S1_seq_4h))

# plot
dir.create(path = paste0(plot_dir, "/", DSB_locus), showWarnings = FALSE)
pdf(file = "tmp.pdf", width = 3, height = 4.5)
plotTracks(trackList = list(AT_genes, OT_0h, OT_1h,  OT_2h, OT_4h), sizes = c(0.05, rep(0.95/4, 4)),
           from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 0)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/", DSB_locus, "/Plot.pdf"))

# save axis ranges
tmp <- data.frame(SrfIcs = DSB_locus, region_size = width(roi),
                  S1_seq_min = ylim_S1_seq[1], S1_seq_max = ylim_S1_seq[2],
                  MNase_seq_min = ylim_MNase_seq[1], MNase_seq_max = ylim_MNase_seq[2])
write.table(x = tmp, file = paste0(plot_dir, "/", DSB_locus, "/axis_ranges.txt"), row.names = FALSE)


# }
