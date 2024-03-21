# info --------------------------------------------------------------------
# purpose: helper functions for plotting genomic data with Gviz
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date started: 03/21/24
# last modified: 03/21/24

# make DataTrack with differently colored fw and rev scores ---------------
# argument: GRanges object, color definitions, ... (e.g. type, name)
# result: OverlayTrack
DataTrack_fw_rev <- function(GRanges, fw_col, rev_col, ...){
  fw <- DataTrack(range = GRanges[strand(GRanges) == "+"], col = fw_col, ...)
  rev <- DataTrack(range = GRanges[strand(GRanges) == "-"], col = rev_col, ...)
  OverlayTrack(trackList = list(fw, rev))
}

# set 0 to NA in GRanges object -------------------------------------------
# to avoid plotting of line at y = 0
# argument: GRanges object
# result: GRanges object
GRanges_zero_to_NA <- function(GRanges){
  out <- GRanges
  out$score[out$score == 0] <- NA
  return(out)
}