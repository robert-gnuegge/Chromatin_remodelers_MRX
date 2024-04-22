# info --------------------------------------------------------------------
# purpose: plot DSB formation data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/29/24
# last modified: 03/29/24

# read helper functions and files -----------------------------------------
source(file = "../../Src/JFly_colors.R")
source(file = "../../Src/Misc_helper_functions.R")


# read data ---------------------------------------------------------------
Cut_fraction <- read.table(file = "03_Processed_data/SrfIcs_cutting.txt", header = TRUE)

# plotting function -------------------------------------------------------
my_plot <- function(data, strains = NULL, colors, file_name, y_range = NULL, jitter = 0.2){
  
  if(is.null(y_range)){
    y_range <- range(c(data$mean + data$sd, data$mean - data$sd, 0, 1), na.rm = TRUE)
  }
  
  if(is.null(strains)){
    strains <- unique(data$strain)
  }
  
  # print to PDF
  pdf(file = "tmp.pdf", width=2.5, height=2.25)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(3.8, 2, 4.0, 2.0), tcl = -0.25, mgp = c(2.25, 0.4, 0), las = 1)
  
  # start empty plot
  plot(x = NA, y = NA, ylim = y_range, xlim = range(data$time), xlab = NA, ylab = NA)
  
  # iterate through all samples
  for (n in 1:length(strains)){
    
    strain <- strains[n]
    
    # time values with jitter
    t <- jitter(x = data$time[data$strain == strain], factor = jitter)
    
    # line plot of means
    points(x = t, 
           y = data$mean[data$strain == strain],
           pch = 20, col = colors[n], type = "o"
    )
    
    # add error bars
    arrows(x0 = t,
           y0 = data$mean[data$strain == strain] - data$sd[data$strain == strain],
           x1 = t,
           y1 = data$mean[data$strain == strain] + data$sd[data$strain == strain],
           length = 0.03, # length of arrow head
           angle = 90, # angle of arrow head
           code = 3, # to draw arrow head on both ends
           col = colors[n])
    
  }
  
  dev.off()
  
  GS_embed_fonts(input = "tmp.pdf", output = file_name)
  
}


# plotting ================================================================
plot_dir <- "04_Plots/"
plot_colors <- c("black", "gray", JFly_colors[c(4, 2)])
strains <- c("5935_1", "5935_2", "5934_1", "5934_2")

# plot legend -------------------------------------------------------------
Legend_txt <- c("WT #1", "WT #2",
                expression(italic("fun30"*Delta)~"#1"),
                expression(italic("fun30"*Delta)~"#2"))

pdf(file = "tmp.pdf", width=1.20, height=1.05)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, col = plot_colors, pch = 20)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Legend.pdf"))


# plot graphs -------------------------------------------------------------
primers_amplicons <- data.frame(primers = c("oRG46_oRG47", "oRG341_oRG342", "oRG345_oRG346", "oRG347_oRG348", "oRG353_oRG354", 
                                            "oRG355_oRG356", "oRG361_oRG362", "oRG1075_oRG1076", "oRG1079_oRG1080", "oRG1081_oRG1082", 
                                            "oRG1085_oRG1086", "oRG1091_oRG1092", "oRG1093_oRG1094", "oRG1099_oRG1100", "oRG1101_oRG1102", 
                                            "oRG1107_oRG1108", "oRG1109_oRG1110", "oRG1115_oRG1116", "oRG1119_oRG1120", "oRG1121_oRG1122"),
                                amplicon = c("Chr3_200908", "Chr2_256173", "Chr13_664938", "Chr15_27760", "Chr5_399646", 
                                             "Chr7_642694", "Chr15_370687", "Chr4_370477", "Chr5_123366", "Chr7_398052", 
                                             "Chr7_517001", "Chr10_360907", "Chr12_498747", "Chr13_129799", "Chr13_566584", 
                                             "Chr13_676597", "Chr14_301414", "Chr15_756594", "Chr15_1039563", "Chr16_431983"))

for(n in 1:nrow(primers_amplicons)){
  
  cat("\nPlotting", primers_amplicons$amplicon[n], "...")
  
  tmp <- Cut_fraction[Cut_fraction$primers == primers_amplicons$primers[n], ]
  my_plot(data = tmp, strains = strains, colors = plot_colors, file_name = paste0(plot_dir, primers_amplicons$amplicon[n], ".pdf"))
  
}

# what is the order of cutting kinetics?
Final_cut_level <- aggregate(mean ~ primers, data = Cut_fraction[Cut_fraction$time == 4, ], FUN = mean)
Final_cut_level <- Final_cut_level[order(Final_cut_level$mean, decreasing = TRUE), ]

primers_amplicons$amplicon[match(x = Final_cut_level$primers, table = primers_amplicons$primers)]
#  [1] "Chr3_200908"   "Chr12_498747"  "Chr10_360907"  "Chr14_301414"  "Chr2_256173"   "Chr15_27760"   "Chr5_123366"   "Chr15_1039563" "Chr13_664938" 
# [10] "Chr15_370687"  "Chr7_398052"   "Chr5_399646"   "Chr16_431983"  "Chr13_566584"  "Chr7_517001"   "Chr13_676597"  "Chr13_129799"  "Chr4_370477"  
# [19] "Chr15_756594"  "Chr7_642694"  


plot(x = Cut_fraction$mean[Cut_fraction$strain == "5935_1"], y = Cut_fraction$mean[Cut_fraction$strain == "5935_2"])
cor(x = Cut_fraction$mean[Cut_fraction$strain == "5935_1"], y = Cut_fraction$mean[Cut_fraction$strain == "5935_2"], method = "pearson")
lm(y ~ x, data = data.frame(x = Cut_fraction$mean[Cut_fraction$strain == "5935_1"], y = Cut_fraction$mean[Cut_fraction$strain == "5935_2"]))


plot(x = Cut_fraction$mean[Cut_fraction$strain == "5934_1"], y = Cut_fraction$mean[Cut_fraction$strain == "5934_2"])
cor(x = Cut_fraction$mean[Cut_fraction$strain == "5934_1"], y = Cut_fraction$mean[Cut_fraction$strain == "5934_2"], method = "pearson")
lm(y ~ x, data = data.frame(x = Cut_fraction$mean[Cut_fraction$strain == "5934_1"], y = Cut_fraction$mean[Cut_fraction$strain == "5934_2"]))
