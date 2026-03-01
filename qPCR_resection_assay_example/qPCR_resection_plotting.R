# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# create directory for plots
dir.create(path = "Plots/", showWarnings = FALSE)

# read in helper files and functions
source(file = "/home/robert/Research/Software/R_scripts/JFly_colors.R")
source(file = "/home/robert/Research/Software/R_scripts/Misc_helper_functions.R")


# read data ===============================================================
Cut_fraction <- read.table(file = "Moments/Cutting.txt", header = TRUE, stringsAsFactors = FALSE)
Resected_fraction <- read.table(file = "Moments/Resection.txt", header = TRUE, stringsAsFactors = FALSE)


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
plot_dir <- "Plots/"
dir.create(path = plot_dir, showWarnings = FALSE)
strains <- c("6072", "6073", "6074", "6075")
MyColors <- c("gray", "black", JFlyColors[c(2, 4)])

# plot legend -------------------------------------------------------------
Legend_txt <- c(expression("2"*mu),
                expression("2"*mu*"-"*italic("SAE2")),
                expression(italic("fun30"*Delta)~"2"*mu),
                expression(italic("fun30"*Delta)~"2"*mu*"-"*italic("SAE2")))

pdf(file = "tmp.pdf", width=1.7, height=1.05)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Legend_txt, col = MyColors, pch = 20)
dev.off()
GS_embed_fonts(input = "tmp.pdf", output = paste0(plot_dir, "Legend.pdf"))

# cutting -----------------------------------------------------------------
tmp <- Cut_fraction[Cut_fraction$strain %in% strains & Cut_fraction$primers == "oRG46_oRG47", ]
my_plot(data = tmp, strains = strains, colors = MyColors, ylab = "Cut Fraction", file_name = paste0(plot_dir, "Cut_fraction.pdf"))

# +98 bp ----------------------------------------------------------------------
tmp <- Resected_fraction[Resected_fraction$strain %in% strains & Resected_fraction$primers == "oRG50_oRG51", ]
my_plot(data = tmp, strains = strains, colors = MyColors, ylab = "ssDNA Fraction", file_name = paste0(plot_dir, "Resection_98_bp.pdf"))

# +640 bp ----------------------------------------------------------------------
tmp <- Resected_fraction[Resected_fraction$strain %in% strains & Resected_fraction$primers == "oRG52_oRG53", ]
my_plot(data = tmp, strains = strains, colors = MyColors, ylab = "ssDNA Fraction", file_name = paste0(plot_dir, "Resection_640_bp.pdf"))
