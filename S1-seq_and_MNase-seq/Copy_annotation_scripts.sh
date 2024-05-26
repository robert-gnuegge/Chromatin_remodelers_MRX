#!/bin/bash

declare -a arr=("chrII_256173" "chrIII_200753" "chrIV_370477" "chrV_123366" "chrV_399646" "chrVII_398052" "chrVII_517001" "chrVII_642694" "chrX_360907" "chrXII_498747" "chrXIII_129799" "chrXIII_566584" "chrXIII_664938" "chrXIII_676597" "chrXIV_301414" "chrXV_1039563" "chrXV_370687" "chrXV_756594" "chrXVI_431983")

for S in "${arr[@]}"
do

	cp "/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/S1-seq_and_MNase-seq/Rep1/Plots/04_Plots/MRE11/S1-seq_and_MNase-seq_Plots/$S/Annotated_Plot.tex" \
	"/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/S1-seq_and_MNase-seq/Rep_merged/Analyze_and_plot/04_Plots/MRE11/S1-seq_and_MNase-seq_Plots/$S/"

done