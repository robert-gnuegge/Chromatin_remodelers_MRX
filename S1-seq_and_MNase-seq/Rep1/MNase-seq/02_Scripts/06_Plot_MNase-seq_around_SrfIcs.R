# info --------------------------------------------------------------------
# purpose: plot MNase-seq coverage around SrfIcs after read trimming
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/27/24
# last modified: 03/27/24

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


# load MNase-seq data -----------------------------------------------------
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_147bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_127bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_117bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_107bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_97bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_87bp.RData")
load(file = "../MNase-seq/03_Processed_data/MNase-seq_coverage/MNase-seq_74bp.RData")

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
