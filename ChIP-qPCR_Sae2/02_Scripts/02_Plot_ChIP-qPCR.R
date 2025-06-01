# info --------------------------------------------------------------------
# purpose: plot ChIP-qPCR data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 02/25/25
# last modified: 03/02/25

# read helper functions and files -----------------------------------------
source(file = "../Src/JFly_colors.R")
source(file = "../Src/Misc_helper_functions.R")

# read data ---------------------------------------------------------------
perc_in_data <- read.table(file = "03_Processed_data/Percent_input_rep1.txt", header = TRUE)

# plotting function -------------------------------------------------------
my_plot <- function(data, strains, colors, y_range = NULL){
  
  if(is.null(y_range)){
    y_range <- range(c(data$mean + data$sd, data$mean - data$sd, 0), na.rm = TRUE)
    y_range[1] <- y_range[1] - 0.02 * y_range[2]  # increase slightly to prevent partial plotting of data points at plot border
    y_range[2] <- 1.02 * y_range[2]
  }
  
  # prepare average data for plotting
  data$sort <- match(x = data$strain, table = strains)
  data <- data[order(data$sort, data$time), ]
  heights <- matrix(data = data$mean, nrow = length(strains), byrow = TRUE)
  
  # plot barplot
  bp <- barplot(height = heights, beside = TRUE, ylim = y_range, xlab = "Time (h)", ylab = NA, col = colors)
  axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = unique(data$time), line = -0.4, lwd = 0)
  
  # add error bars (sd)
  x <- as.vector(t(bp))
  arrows(x0 = x,
         y0 = data$mean - data$sd,
         x1 = x,
         y1 = data$mean + data$sd,
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3 # to draw arrow head on both ends
  )
  
}

# plotting ================================================================
plot_dir <- "04_Plots/"

strains <- c("LSY6098", "LSY6097")
MyColors <- c("gray", JFly_colors[2])


# plot legend -------------------------------------------------------------
Legend_txt <- c("WT", expression(italic("fun30"*Delta)))

pdf(file = "tmp.pdf", width=1.1, height=0.65)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, fill = MyColors)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Legend.pdf"))

pdf(file = "tmp.pdf", width=2.05, height=0.45)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, fill = MyColors, ncol = 2)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Legend_horiz.pdf"))


# +98 bp ----------------------------------------------------------------------
avg <- perc_in_data[perc_in_data$strain %in% strains & perc_in_data$time != 2 & perc_in_data$primers == "oRG50_oRG51", ]

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.5, 4.1, 2.1), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

my_plot(data = avg, y_range = NULL, strains = strains, colors = MyColors)
title(ylab = "% Input", line = 2.75)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Percent_input_98_bp.pdf"))

# +640 bp ---------------------------------------------------------------------
avg <- perc_in_data[perc_in_data$strain %in% strains & perc_in_data$time != 2 & perc_in_data$primers == "oRG52_oRG53", ]

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0, 4.1, 2.1), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

my_plot(data = avg, y_range = NULL, strains = strains, colors = MyColors)
title(ylab = "% Input", line = 3.25)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Percent_input_647_bp.pdf"))


# ADH1 --------------------------------------------------------------------
avg <- perc_in_data[perc_in_data$strain %in% strains & perc_in_data$time != 2 & perc_in_data$primers == "oRG54_oRG55", ]

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0, 4.1, 2.1), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

my_plot(data = avg, y_range = NULL, strains = strains, colors = MyColors)
title(ylab = "% Input", line = 3.25)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Percent_input_ADH1.pdf"))


# Sae2/Mre11 ratio ==========================================================

perc_in_Mre11 <- read.table(file = "../ChIP-qPCR_Mre11/03_Processed_data/percent_input_averages.txt", header = TRUE)

# +98 bp ----------------------------------------------------------------------
sae2 <- perc_in_data[perc_in_data$strain %in% strains & perc_in_data$time != 2 & perc_in_data$primers == "oRG50_oRG51", ]
sae2 <- sae2[order(sae2$strain, decreasing = TRUE),]
mre11 <- perc_in_Mre11[perc_in_Mre11$strain %in% c("4518-13B", "5415") & perc_in_Mre11$sample == "IP" & perc_in_Mre11$time != 2 & perc_in_Mre11$primers == "oRG50_oRG51", ]

avg <- data.frame(strain = c(rep("WT", 3), rep("fun30", 3)),
                  time = c(0, 1, 4), 
                  ratio = sae2$mean / mre11$mean)

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.5, 4, 2.1), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

heights <- matrix(data = avg$ratio[avg$time %in% c(1, 4)], nrow = 2, byrow = TRUE)

# plot barplot
bp <- barplot(height = heights, beside = TRUE, xlab = "Time (h)", ylab = NA, col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(1, 4), line = -0.4, lwd = 0)

title(ylab = "Sae2/Mre11 (AU)", line = 2.75)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Sae2_over_Mre11_98_bp.pdf"))


# +640 bp ----------------------------------------------------------------------
sae2 <- perc_in_data[perc_in_data$strain %in% strains & perc_in_data$time != 2 & perc_in_data$primers == "oRG52_oRG53", ]
sae2 <- sae2[order(sae2$strain, decreasing = TRUE),]
mre11 <- perc_in_Mre11[perc_in_Mre11$strain %in% c("4518-13B", "5415") & perc_in_Mre11$sample == "IP" & perc_in_Mre11$time != 2 & perc_in_Mre11$primers == "oRG52_oRG53", ]

avg <- data.frame(strain = c(rep("WT", 3), rep("fun30", 3)),
                  time = c(0, 1, 4), 
                  ratio = sae2$mean / mre11$mean)

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0, 4, 2.1), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

heights <- matrix(data = avg$ratio[avg$time %in% c(1, 4)], nrow = 2, byrow = TRUE)

# plot barplot
bp <- barplot(height = heights, beside = TRUE, xlab = "Time (h)", ylab = NA, col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(1, 4), line = -0.4, lwd = 0)

title(ylab = "Sae2/Mre11 (AU)", line = 3.25)

dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Sae2_over_Mre11_647_bp.pdf"))
