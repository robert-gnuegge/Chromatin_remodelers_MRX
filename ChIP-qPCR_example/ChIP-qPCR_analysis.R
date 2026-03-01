# set working directory to this file's location
wd.path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd.path)

# create directories for plots and distribution data
dir.create(path = "Plots/", showWarnings = FALSE)
dir.create(path = "Moments/", showWarnings = FALSE)

# read in helper files and functions
source(file = "/home/robert/Research/Software/R_scripts/JFly_colors.R")
source(file = "/home/robert/Research/Software/R_scripts/Misc_functions.R")
source(file = "/home/robert/Research/Software/R_scripts/ChIP-qPCR.R")

# read primer efficiencies
efficiencies <- read.table(file = "/home/robert/Research/Resources/qPCR_primer_efficiencies/primer_efficiencies.txt", header = TRUE)

# read, check, and summarize Cq values ========================================

# define plate layout ---------------------------------------------------------
qPCR <- data.frame(plate = c(rep(1, 18 * 14), rep(2, 18 * 12 + 9)),
                   well = c(paste0(rep(LETTERS[2:5], rep(9, 4)), formatC(x = 2:10, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[2:5], rep(9, 4)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[6:9], rep(9, 4)), formatC(x = 2:10, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[6:9], rep(9, 4)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[10:13], rep(9, 4)), formatC(x = 2:10, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[10:13], rep(9, 4)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[14:15], rep(9, 2)), formatC(x = 2:10, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[14:15], rep(9, 2)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[2:5], rep(9, 4)), formatC(x = 2:10, width = 2, format = "d", flag = "0")), # plate 2
                            paste0(rep(LETTERS[2:5], rep(9, 4)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[6:9], rep(9, 4)), formatC(x = 2:10, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[6:9], rep(9, 4)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[10:13], rep(9, 4)), formatC(x = 2:10, width = 2, format = "d", flag = "0")),
                            paste0(rep(LETTERS[10:13], rep(9, 4)), formatC(x = 11:19, width = 2, format = "d", flag = "0")),
                            paste0("N", formatC(x = 2:10, width = 2, format = "d", flag = "0"))),
                   strain = c(rep(c("4518-13B", "5415", "5495-30D"), rep(72, 3)), rep("4518-13B", 36),
                              rep(c("5692-15B", "5495-30D (04/12)", "5692-15B (04/12)"), rep(72, 3)), rep("NTC", 9)),
                   sample = c(rep(c("Input", "IP", "Input", "IP", "Input", "IP"), rep(36, 6)), rep("No Ab", 36),
                              rep(c("Input", "IP", "Input", "IP", "Input", "IP"), rep(36, 6)), rep("NTC", 9)),
                   time = c(rep(rep(c(0, 1, 2, 4), rep(9, 4)), 7),
                            rep(rep(c(0, 1, 2, 4), rep(9, 4)), 6), rep(0, 9)),
                   primers = c(rep(rep(c("oRG54_oRG55", "oRG50_oRG51", "oRG52_oRG53"), rep(3, 3)), 28),
                               rep(rep(c("oRG54_oRG55", "oRG50_oRG51", "oRG52_oRG53"), rep(3, 3)), 25)),
                   Cq = 0)

head(qPCR)
tail(qPCR, n = 40)

# read Cq Results file --------------------------------------------------------
file_path <- grep(pattern = "Cq", x = list.files(path = "qPCR_data", full.names = TRUE, recursive = TRUE), value = TRUE)

# save Cq values in qPCR data frame
for (plate in 1:2){
  tmp <- read.table(file = file_path[plate], sep = "\t", header = TRUE)
  for (well in qPCR$well[qPCR$plate == plate]){
    qPCR$Cq[qPCR$well == well & qPCR$plate == plate] <- tmp$Cq[tmp$Well == well & tmp$Fluor == "SYBR"]
  }  
}

head(qPCR)
tail(qPCR, n = 20)

# check NTCs ------------------------------------------------------------------

# min NTC Cq
min(qPCR$Cq[qPCR$sample == "NTC"], na.rm = TRUE)
# 33.6317

# max non-NTC Cq
max(qPCR$Cq[qPCR$sample != "NTC"], na.rm = TRUE)
# 36.19696

# calc Cq mean, sd, and cv ----------------------------------------------------
Cq.modes <- aggregate(Cq ~ strain + sample + time + primers, data = qPCR,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          cv = sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                                          )
                      )

# sort
Cq.modes <- Cq.modes[order(Cq.modes$strain, Cq.modes$sample, Cq.modes$primers, Cq.modes$time),]

# mean, sd, and cv are in a single column (as a matrix)
# convert into separate columns
Cq.modes <- data.frame(Cq.modes[, 1:(ncol(Cq.modes)-1)], as.data.frame(Cq.modes$Cq))

# what is the largest cv?
max(Cq.modes$cv, na.rm = TRUE)
# 0.04986006

# which are the samples with large cv?
Cq.modes[Cq.modes$cv > 0.025, ]
#               strain sample time     primers     mean        sd         cv
# 80          4518-13B  No Ab    1 oRG52_oRG53 34.59452 1.7248846 0.04986006
# 143 5495-30D (04/12)     IP    2 oRG54_oRG55 33.92950 1.2559145 0.03701541
# 79  5692-15B (04/12)     IP    1 oRG52_oRG53 26.39536 0.7019615 0.02659412
# 120              NTC    NTC    0 oRG54_oRG55 38.63764 1.3220484 0.03421659

# what is the ADH1 fold spread in Input samples?
tmp <- Cq.modes$mean[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample == "Input"]
efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(tmp)))
# 11.66166

# what is the ADH1 fold spread within each strain
aggregate(mean ~ strain, data = Cq.modes[Cq.modes$primers == "oRG54_oRG55" & Cq.modes$sample == "Input", ],
          FUN = function(x) efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"]^(diff(range(x))))
#             strain     mean
# 1         4518-13B 2.505153
# 2             5415 1.790438
# 3         5495-30D 4.524464
# 4 5495-30D (04/12) 6.495585
# 5         5692-15B 5.880518
# 6 5692-15B (04/12) 2.909435


# calc percent input enrichment ===========================================

# initialize data.frame to collect results
percent_input <- data.frame()

for (strain in c("4518-13B", "5415", "5495-30D", "5692-15B", "5495-30D (04/12)", "5692-15B (04/12)")){
  for (primers in c("oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55")){
    for (time in c(0, 1, 2, 4)){
      tmp <- Calc_percent_input_err_propagated(Cq_Input_mean_and_sd = Cq.modes[Cq.modes$strain == strain & Cq.modes$sample == "Input" & Cq.modes$primers == primers & Cq.modes$time == time, c("mean", "sd")],
                                               Cq_IP_mean_and_sd = Cq.modes[Cq.modes$strain == strain & Cq.modes$sample == "IP" & Cq.modes$primers == primers & Cq.modes$time == time, c("mean", "sd")],
                                               Input_dil = 0.01,
                                               E = efficiencies$efficiency[efficiencies$primers == primers])
      percent_input <- rbind(percent_input,
                             data.frame(strain = strain, sample = "IP", primers = primers, time = time, tmp))
    }
  }
}

# add -Ab ctrl
for (strain in c("4518-13B")){
  for (primers in c("oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55")){
    for (time in c(0, 1, 2, 4)){
      tmp <- Calc_percent_input_err_propagated(Cq_Input_mean_and_sd = Cq.modes[Cq.modes$strain == strain & Cq.modes$sample == "Input" & Cq.modes$primers == primers & Cq.modes$time == time, c("mean", "sd")],
                                               Cq_IP_mean_and_sd = Cq.modes[Cq.modes$strain == strain & Cq.modes$sample == "No Ab" & Cq.modes$primers == primers & Cq.modes$time == time, c("mean", "sd")],
                                               Input_dil = 0.01,
                                               E = efficiencies$efficiency[efficiencies$primers == primers])
      percent_input <- rbind(percent_input,
                             data.frame(strain = strain, sample = "No Ab", primers = primers, time = time, tmp))
    }
  }
}

# write data to file
write.table(x = percent_input, file = "Moments/Percent_input.txt", row.names = FALSE, col.names = TRUE)


# plotting ================================================================
strains <- c("4518-13B", "5415", "5495-30D", "5692-15B")  # lexO-SrfI WT, lexO-SrfI fun30, lexO-HO WT, lexO-HO fun30 
MyColors <- c(JFlyColors[c(7, 3)], "black", "gray", "white") # white = No Ab


# plot legend =============================================================
Leg_txt <- c(expression(italic("lexO-SrfI")), 
             expression(italic("lexO-SrfI fun30")*Delta), 
             expression(italic("lexO-HO")), 
             expression(italic("lexO-HO fun30")*Delta),
             "No Ab")

pdf(file = "tmp.pdf", width=1.85, height=1.25)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Leg_txt, fill = MyColors)
dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots/Legend.pdf"), RemoveInputFile = TRUE)


# plot percent input ======================================================

# calc axis ranges
y_range <- 1.01 * range(c(percent_input$mean + percent_input$sd, percent_input$mean - percent_input$sd, 0), na.rm = TRUE)

# 98 bp -------------------------------------------------------------------
tmp <- percent_input[percent_input$primers == "oRG50_oRG51", ]
tmp <- tmp[tmp$strain %in% strains, ]
heights <- matrix(data = tmp$mean, nrow = 5, byrow = TRUE)

# print to PDF
pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.9, 4.0, 2.0), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

bp <- barplot(height = heights, beside = TRUE, ylim = y_range, ylab = NA, xlab = "Time (h)", col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 2, 4), line = -0.3, lwd = 0)
title(ylab = "% Input", line = 2.25)

# add error bars
x <- as.vector(t(bp))
arrows(x0 = x,
       y0 = tmp$mean - tmp$sd,
       x1 = x,
       y1 = tmp$mean + tmp$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
       )

dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots/Percent_input_98bp.pdf"), RemoveInputFile = TRUE)


# 647 bp ------------------------------------------------------------------
tmp <- percent_input[percent_input$primers == "oRG52_oRG53", ]
tmp <- tmp[tmp$strain %in% strains, ]
heights <- matrix(data = tmp$mean, nrow = 5, byrow = TRUE)

# print to PDF
pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.9, 4.0, 2.0), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

bp <- barplot(height = heights, beside = TRUE, ylim = y_range, ylab = NA, xlab = "Time (h)", col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 2, 4), line = -0.3, lwd = 0)
title(ylab = "% Input", line = 2.25)

# add error bars
x <- as.vector(t(bp))
arrows(x0 = x,
       y0 = tmp$mean - tmp$sd,
       x1 = x,
       y1 = tmp$mean + tmp$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
)

dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots/Percent_input_647bp.pdf"), RemoveInputFile = TRUE)


# ADH1 ctrl ---------------------------------------------------------------
tmp <- percent_input[percent_input$primers == "oRG54_oRG55", ]
tmp <- tmp[tmp$strain %in% strains, ]
heights <- matrix(data = tmp$mean, nrow = 5, byrow = TRUE)


# print to PDF
pdf(file = "tmp.pdf", width=3, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.9, 4.0, 2.0), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

bp <- barplot(height = heights, beside = TRUE, ylim = y_range, ylab = NA, xlab = "Time (h)", col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 2, 4), line = -0.3, lwd = 0)
title(ylab = "% Input", line = 2.25)

# add error bars
x <- as.vector(t(bp))
arrows(x0 = x,
       y0 = tmp$mean - tmp$sd,
       x1 = x,
       y1 = tmp$mean + tmp$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
)

dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots/Percent_input_ADH1.pdf"), RemoveInputFile = TRUE)


# plot percent input for 04/12 samples ====================================
dir.create(path = "Plots_04-12/", showWarnings = FALSE)
strains <- c("5495-30D (04/12)", "5692-15B (04/12)")  # lexO-HO WT, lexO-HO fun30 
MyColors <- gray.colors(n = 3)

# plot legend =============================================================
Leg_txt <- c(expression(italic("lexO-HO")), 
             expression(italic("lexO-HO fun30")*Delta),
             "No Ab")

pdf(file = "tmp.pdf", width=1.8, height=0.85)
par(cex = 1, mar = rep(0, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, xjust=0.5, yjust=0.5, legend = Leg_txt, fill = MyColors)
dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots_04-12/Legend.pdf"), RemoveInputFile = TRUE)


# calc axis ranges
percent_input_04_12 <- percent_input[percent_input$strain %in% c("5495-30D (04/12)", "5692-15B (04/12)"), ]
percent_input_04_12 <- rbind(percent_input_04_12, percent_input[percent_input$sample == "No Ab", ])

y_range <- 1.01 * range(c(percent_input_04_12$mean + percent_input_04_12$sd, percent_input_04_12$mean - percent_input_04_12$sd, 0), na.rm = TRUE)

# 98 bp -------------------------------------------------------------------
tmp <- percent_input_04_12[percent_input_04_12$primers == "oRG50_oRG51", ]
heights <- matrix(data = tmp$mean, nrow = 3, byrow = TRUE)

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.9, 4.0, 2.0), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

bp <- barplot(height = heights, beside = TRUE, ylim = y_range, ylab = NA, xlab = "Time (h)", col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 2, 4), line = -0.3, lwd = 0)
title(ylab = "% Input", line = 2.25)

# add error bars
x <- as.vector(t(bp))
arrows(x0 = x,
       y0 = tmp$mean - tmp$sd,
       x1 = x,
       y1 = tmp$mean + tmp$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
)

dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots_04-12/98bp.pdf"), RemoveInputFile = TRUE)


# 647 bp ------------------------------------------------------------------
tmp <- percent_input_04_12[percent_input_04_12$primers == "oRG52_oRG53", ]
heights <- matrix(data = tmp$mean, nrow = 3, byrow = TRUE)

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.9, 4.0, 2.0), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

bp <- barplot(height = heights, beside = TRUE, ylim = y_range, ylab = NA, xlab = "Time (h)", col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 2, 4), line = -0.3, lwd = 0)
title(ylab = "% Input", line = 2.25)

# add error bars
x <- as.vector(t(bp))
arrows(x0 = x,
       y0 = tmp$mean - tmp$sd,
       x1 = x,
       y1 = tmp$mean + tmp$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
)

dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots_04-12/647bp.pdf"), RemoveInputFile = TRUE)


# ADH1 ctrl ---------------------------------------------------------------
tmp <- percent_input_04_12[percent_input_04_12$primers == "oRG54_oRG55", ]
heights <- matrix(data = tmp$mean, nrow = 3, byrow = TRUE)

# print to PDF
pdf(file = "tmp.pdf", width=2.5, height=2.5)
par(cex = 1, mar = c(5.1, 4.1, 4.1, 2.1) - c(2.5, 0.9, 4.0, 2.0), tcl = -0.3, mgp = c(1.5, 0.6, 0), las = 1)

bp <- barplot(height = heights, beside = TRUE, ylim = y_range, ylab = NA, xlab = "Time (h)", col = MyColors)
axis(side = 1, at = apply(X = bp, MARGIN = 2, FUN = mean), labels = c(0, 1, 2, 4), line = -0.3, lwd = 0)
title(ylab = "% Input", line = 2.25)

# add error bars
x <- as.vector(t(bp))
arrows(x0 = x,
       y0 = tmp$mean - tmp$sd,
       x1 = x,
       y1 = tmp$mean + tmp$sd,
       length = 0.03, # length of arrow head
       angle = 90, # angle of arrow head
       code = 3 # to draw arrow head on both ends
)

dev.off()
GSEmbedFonts(input = "tmp.pdf", output = paste0("Plots_04-12/ADH1.pdf"), RemoveInputFile = TRUE)



# calc fold enrichment over control =======================================

library(errors)

# create data frame with Cq means and assigned errors (as attribute)
Cq_err <- Cq.modes[, c("strain", "sample", "time", "primers", "mean")]
errors(Cq_err$mean) <- Cq.modes$sd

# define funtion
Calc_fold_over_ctrl <- function(Cq_Input_target, Cq_IP_target, Cq_Input_ctrl, Cq_IP_ctrl, Input_dil, E_target, E_ctrl){
  fold_over_ctrl <- E_target^(Cq_Input_target - log(1 / Input_dil, E_target) - Cq_IP_target) / E_ctrl^(Cq_Input_ctrl - log(1 / Input_dil, E_ctrl) - Cq_IP_ctrl)
  return(fold_over_ctrl)
}

# initialize data.frame to collect results
fold_enrich <- data.frame()

for (strain in c("4518-13B", "5415", "5495-30D", "5692-15B", "5495-30D (04/12)", "5692-15B (04/12)")){
  for (primers in c("oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55")){
    for (time in c(0, 1, 2, 4)){
      tmp <- Calc_fold_over_ctrl(Cq_Input_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == primers & Cq_err$time == time], 
                                 Cq_Input_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time], 
                                 Cq_IP_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == primers & Cq_err$time == time], 
                                 Cq_IP_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "IP" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time],
                                 Input_dil = 0.01, 
                                 E_target = efficiencies$efficiency[efficiencies$primers == primers],
                                 E_ctrl = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"])
      fold_enrich <- rbind(fold_enrich,
                           data.frame(strain = strain, sample = "IP", primers = primers, time = time, mean = as.numeric(tmp), sd = errors(tmp)))
    }
  }
}


# add -Ab ctrl
for (strain in c("4518-13B")){
  for (primers in c("oRG50_oRG51", "oRG52_oRG53", "oRG54_oRG55")){
    for (time in c(0, 1, 2, 4)){
      tmp <- Calc_fold_over_ctrl(Cq_Input_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == primers & Cq_err$time == time], 
                                 Cq_Input_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "Input" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time], 
                                 Cq_IP_target = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "No Ab" & Cq_err$primers == primers & Cq_err$time == time], 
                                 Cq_IP_ctrl = Cq_err$mean[Cq_err$strain == strain & Cq_err$sample == "No Ab" & Cq_err$primers == "oRG54_oRG55" & Cq_err$time == time],
                                 Input_dil = 0.01, 
                                 E_target = efficiencies$efficiency[efficiencies$primers == primers],
                                 E_ctrl = efficiencies$efficiency[efficiencies$primers == "oRG54_oRG55"])
      fold_enrich <- rbind(fold_enrich,
                           data.frame(strain = strain, sample = "No Ab", primers = primers, time = time, mean = as.numeric(tmp), sd = errors(tmp)))
    }
  }
}

# write data to file
write.table(x = fold_enrich, file = "Moments/Fold_over_ctrl.txt", row.names = FALSE, col.names = TRUE)
