# info --------------------------------------------------------------------
# purpose: plot resection data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 03/31/24
# last modified: 03/31/24

# read helper functions and files -----------------------------------------
source(file = "../../Src/JFly_colors.R")
source(file = "../../Src/Misc_helper_functions.R")


# read data ---------------------------------------------------------------
Resected_fraction <- read.table(file = "03_Processed_data/Resection.txt", header = TRUE)

# plotting function =======================================================
my_plot <- function(data, strains = NULL, colors, ylab, file_name, y_range = NULL, jitter = 0.2){
  
  if(is.null(y_range)){
    y_range <- range(c(data$mean + data$sd, data$mean - data$sd, 0, 1), na.rm = TRUE)
  }
  
  if(is.null(strains)){
    strains <- unique(data$strain)
  }
  
  # print to PDF
  pdf(file = "tmp.pdf", width=2.75, height=2.5)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.3, 1.0, 4.0, 2.0), tcl = -0.3, mgp = c(2.25, 0.6, 0), las = 1)
  
  # start empty plot
  plot(x = NA, y = NA, ylim = y_range, xlim = range(data$time), xlab = NA, ylab = ylab)
  title(xlab =  "Time (h)", line = 1.75)
  
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
plot_colors <- c("black", "gray", JFly_colors[c(4, 2, 7, 3)])
strains <- c("4518-13B_1", "4518-13B_2", "5935_1", "5935_2", "5934_1", "5934_2")

# plot legend -------------------------------------------------------------
Legend_txt <- c("WT #1", "WT #2",
                expression(italic("mre11-nd")~"#1"),
                expression(italic("mre11-nd")~"#2"),
                expression(italic("mre11-nd fun30"*Delta)~"#1"),
                expression(italic("mre11-nd fun30"*Delta)~"#2"))

pdf(file = "tmp.pdf", width=2.00, height=1.45)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, col = plot_colors, pch = 20)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Legend.pdf"))

pdf(file = "tmp.pdf", width=5.70, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, col = plot_colors, pch = 20, ncol = 3)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Legend_horiz.pdf"))

# +98 bp ----------------------------------------------------------------------
tmp <- Resected_fraction[Resected_fraction$primers == "oRG50_oRG51", ]
my_plot(data = tmp, strains = strains, colors = plot_colors, ylab = "ssDNA Fraction", file_name = paste0(plot_dir, "Resection_98_bp.pdf"))

# +640 bp ----------------------------------------------------------------------
tmp <- Resected_fraction[Resected_fraction$strain %in% strains & Resected_fraction$primers == "oRG52_oRG53", ]
my_plot(data = tmp, strains = strains, colors = plot_colors, ylab = "ssDNA Fraction", file_name = paste0(plot_dir, "Resection_640_bp.pdf"))
