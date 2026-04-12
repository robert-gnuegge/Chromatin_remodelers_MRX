# info --------------------------------------------------------------------
# purpose: plot ChIP-qPCR data
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# created: 05/11/25
# last modified: 02/27/26

# read helper functions and files -----------------------------------------
source(file = "../Src/JFly_colors.R")
source(file = "../Src/Misc_helper_functions.R")


# read data ---------------------------------------------------------------
avg_data <- read.table(file = "03_Processed_data/percent_input_averages.txt", header = TRUE)
avg_data <- avg_data[avg_data$time != 2, ]

raw_data <- read.table(file = "03_Processed_data/percent_input_replicates.txt", header = TRUE)
raw_data <- raw_data[raw_data$time != 2, ]

# plotting function -------------------------------------------------------
my_plot <- function(avg_data, raw_data, strains, colors, file_name, y_range = NULL, width=2.5, height=2.5){
  
  if(is.null(y_range)){
    y_range <- range(c(avg_data$mean + avg_data$sd, avg_data$mean - avg_data$sd, raw_data$mean, 0), na.rm = TRUE)
    y_range[1] <- y_range[1] - 0.02 * y_range[2]  # increase slightly to prevent partial plotting of data points at plot border
    y_range[2] <- 1.02 * y_range[2]
  }

  # print to PDF
  pdf(file = "tmp.pdf", width=width, height=height)
  par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.5, 4.1, 2.1), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)
  
  # prepare average data for plotting
  avg_data$sort <- match(x = avg_data$strain, table = strains)
  avg_data <- avg_data[order(avg_data$sort, avg_data$time), ]
  heights <- matrix(data = avg_data$mean, nrow = length(strains), byrow = TRUE)
  
  # plot barplot
  bp <- barplot(height = heights, beside = TRUE, ylim = y_range, xlab = "Time (h)", ylab = NA, col = colors)
  axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 4), line = -0.4, lwd = 0)
  title(ylab = "% Input", line = 2.75)
  
  # add individual data.points
  x <- as.vector(t(bp))
  x_shift <- min(diff(sort(x))) / 5
  
  raw_data$sort <- match(x = raw_data$strain, table = strains)
  raw_data <- raw_data[order(raw_data$sort, raw_data$time), ]

  points(x = x - x_shift, y = raw_data$mean[raw_data$replicate == 1], pch = 20, cex = 0.67, col = gray(level = 0.4))
  points(x = x + x_shift, y = raw_data$mean[raw_data$replicate == 2], pch = 20, cex = 0.67, col = gray(level = 0.4))
  points(x = x[1:3], y = raw_data$mean[raw_data$replicate == 3], pch = 20, cex = 0.67, col = gray(level = 0.4))
  
  # add error bars (sd)
  arrows(x0 = x,
         y0 = avg_data$mean - avg_data$sd,
         x1 = x,
         y1 = avg_data$mean + avg_data$sd,
         length = 0.03, # length of arrow head
         angle = 90, # angle of arrow head
         code = 3 # to draw arrow head on both ends
  )
  
  dev.off()
  GS_embed_fonts(input = "tmp.pdf", output = file_name)
  
}


# plotting ================================================================
plot_dir <- "04_Plots/"

strains <- c("LSY6098", "LSY6097")
MyColors <- c("gray", JFly_colors[2])


# plot legend -------------------------------------------------------------
Legend_txt <- c(expression(italic("FUN30")), expression(italic("fun30"*Delta)))

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
avg <- avg_data[avg_data$strain %in% strains & avg_data$primers == "oRG50_oRG51", ]
raw <- raw_data[raw_data$strain %in% strains & raw_data$primers == "oRG50_oRG51", ]

wilcox.test(mean ~ strain, data = raw, subset = (raw$time == 1))
wilcox.test(mean ~ strain, data = raw, subset = (raw$time == 4))

my_plot(avg_data = avg, raw_data = raw, y_range = c(-0.001, 0.09),
        strains = strains, colors = MyColors, file_name = paste0(plot_dir, "Percent_input_98_bp.pdf"))

# +647 bp ----------------------------------------------------------------------
avg <- avg_data[avg_data$strain %in% strains & avg_data$primers == "oRG52_oRG53", ]
raw <- raw_data[raw_data$strain %in% strains & raw_data$primers == "oRG52_oRG53", ]

wilcox.test(mean ~ strain, data = raw, subset = (raw$time == 1))
wilcox.test(mean ~ strain, data = raw, subset = (raw$time == 4))

my_plot(avg_data = avg, raw_data = raw, y_range = c(-0.001, 0.09),
        strains = strains, colors = MyColors, file_name = paste0(plot_dir, "Percent_input_647_bp.pdf"))

# ADH1 --------------------------------------------------------------------
avg <- avg_data[avg_data$strain %in% strains & avg_data$primers == "oRG54_oRG55", ]
raw <- raw_data[raw_data$strain %in% strains & raw_data$primers == "oRG54_oRG55", ]

my_plot(avg_data = avg, raw_data = raw, y_range = c(-0.001, 0.09),
        strains = strains, colors = MyColors, file_name = paste0(plot_dir, "Percent_input_ADH1.pdf"))
