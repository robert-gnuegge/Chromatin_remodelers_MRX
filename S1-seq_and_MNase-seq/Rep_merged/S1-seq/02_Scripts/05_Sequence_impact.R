# info --------------------------------------------------------------------
# purpose: analyze nt sequence impact on MRX nicking
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 06/15/25
# last modified: 06/15/25

# load libraries ----------------------------------------------------------
library(GenomicRanges)
library(BSgenome)
library(ggplot2)
library(ggseqlogo)

# read helper functions and files -----------------------------------------
source(file = "../../../Src/JFly_colors.R")
source(file = "../../../Src/Misc_helper_functions.R")
source(file = "../../../Src/Genomic_helper_functions.R")
source(file = "../../../Src/S_cerevisiae_SrfI_cut_sites.R")

# read modified S. cerevisiae genome
genome <- import(con = "../../../Reference_genome/03_Processed_data/S288C_R64-4-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr.fasta")

# function definitions ----------------------------------------------------

# find resection extend
# argument: GRanges
# result: GRanges
crop_to_resection_extend <- function(GRanges, DSBs, max_DSB_dist = 2500, threshold = 0.95){
  out <- DSBs
  out$DSB_id <- as.character(DSBs)
  for(n in 1:length(out)){
    # find resection tract end on - strand
    roi <- DSBs[n]
    strand(roi) <- "-"
    roi <- resize(x = roi, width = max_DSB_dist, fix = "start")
    tmp <- subsetByIntersect(subject = GRanges, query = roi)
    tmp <- rev(tmp)  # for increasing distance_to_DSB values
    idx <- which.max(cumsum(tmp$score) >= threshold * max(cumsum(tmp$score)))
    start(out[n]) <- start(tmp[idx])
    # find tract end on + strand
    roi <- DSBs[n]
    strand(roi) <- "+"
    roi <- resize(x = roi, width = max_DSB_dist, fix = "start")
    tmp <- subsetByIntersect(subject = GRanges, query = roi)
    idx <- which.max(cumsum(tmp$score) >= threshold * max(cumsum(tmp$score)))
    end(out[n]) <- end(tmp[idx])
  }
  subsetByIntersect(subject = GRanges, query = out)
}

# calculate position weight matrices
# argument: GRanges, even integer, genome (DNAStringSet)
# result: matrix
calc_PWM <- function(GRanges, width, genome, calc_beakground_PWM = FALSE){
  if(calc_beakground_PWM){
    nick_pos <- GRanges
  }else{
    nick_pos <- GRanges[rep(1:length(GRanges), GRanges$score)]
    # score counts how often a nick occurred at a specific genome position
    # let's repeat score times to retrieve the sequence context score times (see below)
  }
  if(width %% 2 == 1){
    width <- width + 1  
    warning("width must be an even integer. Changing width to ", width)
  }
  nick_context <- flank(x = nick_pos, width = 0.5 * width, both = TRUE)
  nt_sequences <- getSeq(x = genome, names = nick_context)
  PWM <- consensusMatrix(x = nt_sequences, as.prob = TRUE)[1:4, ]
  colnames(PWM) <- c(-(0.5 * width):-1, 1:(0.5 * width))
  return(PWM)
}


# process data for plotting ===============================================

# read unnormalized S1-seq coverage data
load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY4518-13B_S1-seq_unnormalized.RData")
load(file = "../../Rep_merged/S1-seq/03_Processed_data/S1-seq_coverage/LSY5415_S1-seq_unnormalized.RData")

# retain GRanges in regions around DSBs where resection was detected
# other regions are excluded to prevent pattern "dilution" by background noise
WT_1h <- crop_to_resection_extend(GRanges = LSY4518_13B_1h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
WT_2h <- crop_to_resection_extend(GRanges = LSY4518_13B_2h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
WT_4h <- crop_to_resection_extend(GRanges = LSY4518_13B_4h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])

fun30_1h <- crop_to_resection_extend(GRanges = LSY5415_1h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
fun30_2h <- crop_to_resection_extend(GRanges = LSY5415_2h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])
fun30_4h <- crop_to_resection_extend(GRanges = LSY5415_4h_S1_seq_unnormalized, DSBs = SrfIcs[- c(9, 17)])

# remove SrfIcs sequence regions
# they give a lot of S1-seq signal, but derive from unprocessed DSBs and not from MRX nicking
DSBs <- DSB_regions(DSBs = SrfIcs, region_width = 8)
unique(getSeq(x = genome, names = DSBs))  # they are indeed all the SrfIcs sequence

WT_1h <- WT_1h[-subjectHits(findOverlaps(query = DSBs, subject = WT_1h))]
WT_2h <- WT_2h[-subjectHits(findOverlaps(query = DSBs, subject = WT_2h))]
WT_4h <- WT_4h[-subjectHits(findOverlaps(query = DSBs, subject = WT_4h))]

fun30_1h <- fun30_1h[-subjectHits(findOverlaps(query = DSBs, subject = fun30_1h))]
fun30_2h <- fun30_2h[-subjectHits(findOverlaps(query = DSBs, subject = fun30_2h))]
fun30_4h <- fun30_4h[-subjectHits(findOverlaps(query = DSBs, subject = fun30_4h))]

# calculate WT PWM and background PWM
PWM_WT <- calc_PWM(GRanges = c(WT_1h, WT_2h, WT_4h), width = 100, genome = genome)
PWM_WT_bg <- calc_PWM(GRanges = c(WT_1h, WT_2h, WT_4h), width = 100, genome = genome, calc_beakground_PWM = TRUE)

# calculate fun30 PWM and background PWM
PWM_fun30 <- calc_PWM(GRanges = c(fun30_1h, fun30_2h, fun30_4h), width = 100, genome = genome)
PWM_fun30_bg <- calc_PWM(GRanges = c(fun30_1h, fun30_2h, fun30_4h), width = 100, genome = genome, calc_beakground_PWM = TRUE)


# plotting ================================================================
plot_dir <- "04_Plots/Sequence_impact"
dir.create(path = plot_dir, recursive = TRUE, showWarnings = FALSE)

# plot sequence logo ------------------------------------------------------

col_scheme = make_col_scheme(chars = c("A", "T", "G", "C"), cols = JFly_colors[1:4])  # color assignment

# plotting function
plot_seq_logo <- function(PWM){
  PWM <- PWM[, which(colnames(PWM) == "-25"):which(colnames(PWM) == "25")]
  x.ticks <- c(-25, -20, -15, -10, -5, 5, 10, 15, 20, 25)
  ggplot() + geom_logo(PWM, method = 'prob', col_scheme = col_scheme) + theme_logo() +
    scale_x_continuous(breaks = match(x = x.ticks, table = as.numeric(colnames(PWM))), labels = x.ticks) + 
    labs(x = "", y = "") +
    # theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    # theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
    geom_vline(xintercept = 25.5, linetype = "solid", col = "gray")
}

plot_seq_logo(PWM = PWM_WT)
ggsave(filename = "tmp.pdf", device = "pdf", width = 5, height = 2, units = "in")
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_logo_WT.pdf"))

plot_seq_logo(PWM = PWM_WT_bg)
ggsave(filename = "tmp.pdf", device = "pdf", width = 5, height = 2, units = "in")
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_logo_WT_bg.pdf"))

plot_seq_logo(PWM = PWM_fun30)
ggsave(filename = "tmp.pdf", device = "pdf", width = 5, height = 2, units = "in")
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_logo_fun30.pdf"))

plot_seq_logo(PWM = PWM_fun30_bg)
ggsave(filename = "tmp.pdf", device = "pdf", width = 5, height = 2, units = "in")
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "/Sequence_logo_fun30_bg.pdf"))
