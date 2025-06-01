# info --------------------------------------------------------------------
# purpose: investigate MNase-seq coverage at select regions
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/08/24
# last modified: 05/08/24

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

# load and process MNase-seq data -----------------------------------------
load(file = "03_Processed_data/MNase-seq_coverage/LSY4518-13B_MNase-seq.RData")
load(file = "03_Processed_data/MNase-seq_coverage/LSY5415_MNase-seq.RData")
# make LSY5415 scores negative to plot below LSY4518-13B
# LSY5415_0h_MNase_seq$score <- -LSY5415_0h_MNase_seq$score
# LSY5415_1h_MNase_seq$score <- -LSY5415_1h_MNase_seq$score
# LSY5415_2h_MNase_seq$score <- -LSY5415_2h_MNase_seq$score
# LSY5415_4h_MNase_seq$score <- -LSY5415_4h_MNase_seq$score

# AnnotationTrack ---------------------------------------------------------
# data downloaded from yeastmine.yeastgenome.org and modified
load(file = "../../../Src/S_cerevisiae_genome_features.RData")

# adjust for use with AnnotationTrack
names(mcols(all_features))[1] <- "id"

AT_genes <- AnnotationTrack(range = all_features[!(all_features$id == "") & all_features$type %in% c("ORF", "tRNA gene", "centromere", "ARS")], 
                            name = NULL, showFeatureId = FALSE, cex = 0.67, featureAnnotation = "id",
                            arrowHeadMaxWidth = 10, fill = "white", col = "gray", fontcolor.item = "gray")


# plotting ----------------------------------------------------------------
plot_dir <- "04_Plots/Fun30_impact_on_CEN_chromatin/"
dir.create(path = plot_dir, showWarnings = FALSE, recursive = TRUE)

# set global Gviz parameters for plotting
options(Gviz.scheme="default")
scheme <- getScheme()  # copy current scheme
scheme$GdObject$showAxis <- FALSE
scheme$GdObject$showTitle <- FALSE
addScheme(scheme, "MyScheme")  # define new scheme
options(Gviz.scheme = "MyScheme")  # set new scheme

# plotting colors
MNase_seq_col_1 <- NA
MNase_seq_border_col_1 <- "black"
MNase_seq_col_2 <- NA
MNase_seq_border_col_2 <- "gray"

# to calc nice ylim for data
calc_ylim <- function(GRanges, roi, symmetric = FALSE){
  tmp <- range(subsetByIntersect(subject = GRanges, query = roi)$score, na.rm = TRUE)
  if(symmetric){
    tmp <- max(abs(tmp))
    tmp <- c(-tmp, tmp)
  }
  return(tmp)
}

# plot coverage around centromeres
rois <- resize(x = all_features[all_features$type == "centromere"], width = 2000, fix = "center")

# calculate MNase-seq ylim
all_MNase_seq <- c(LSY4518_13B_0h_MNase_seq, LSY5415_0h_MNase_seq)

ylim_MNase_seq <- calc_ylim(GRanges = all_MNase_seq, roi = rois, symmetric = FALSE)

# generate MNase-seq DataTracks and merge WT and fun30 into OverlayTracks
DT_MNase_seq_WT_0h <- DataTrack(range = LSY4518_13B_0h_MNase_seq, type = "polygon", col = MNase_seq_border_col_1, fill.mountain = rep(MNase_seq_col_1, 2), name = "0 h", ylim = ylim_MNase_seq)
DT_MNase_seq_fun30_0h <- DataTrack(range = LSY5415_0h_MNase_seq, type = "polygon", col = MNase_seq_border_col_2, fill.mountain = rep(MNase_seq_col_2, 2), name = "0 h", ylim = ylim_MNase_seq)

OT_MNase_seq_0h <- OverlayTrack(trackList = list(DT_MNase_seq_WT_0h, DT_MNase_seq_fun30_0h))

# plot
pdf(file = "tmp.pdf", width = 3, height = 3)
n <- 2 
roi <- rois[n]
plotTracks(trackList = list(AT_genes, OT_MNase_seq_0h), #sizes = c(0.1, 0.9),
           from = start(roi), to = end(roi), chromosome = seqnames(roi), margin = 0, showTitle = FALSE)
# scheme$GdObject$showAxis <- FALSE was defined above. 
# But without "showTitle = FALSE", a gray title area is plotted anyways. Bug?
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/CEN", n, ".pdf"))


# # save axis ranges
# tmp <- data.frame(SrfIcs = DSB_locus, region_size = width(roi),
#                   S1_seq_min = ylim_S1_seq[1], S1_seq_max = ylim_S1_seq[2],
#                   MNase_seq_min = ylim_MNase_seq[1], MNase_seq_max = ylim_MNase_seq[2])
# write.table(x = tmp, file = paste0(plot_dir, "/", DSB_locus, "/axis_ranges.txt"), row.names = FALSE)
