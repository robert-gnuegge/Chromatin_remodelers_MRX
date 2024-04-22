#!/bin/bash

cd /home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/S1-seq_and_MNase-seq/Rep2/S1-seq/01_Raw_data
mkdir -p Merged_FastQ_files

declare -a arr=("LSY4518-13B_0h" "LSY4518-13B_1h" "LSY4518-13B_2h" "LSY4518-13B_4h"
	"LSY5415_0h" "LSY5415_1h" "LSY5415_2h" "LSY5415_4h")

for OUT in "${arr[@]}"
do

	printf "Merging reads for $OUT ...\n"

	mkdir ./Merged_FastQ_files/$OUT
	find ./BaseSpace_Download/ -print | grep fastq | grep $OUT | grep "R1" | xargs cat > ./Merged_FastQ_files/$OUT/$OUT"_R1.fastq.gz"
	find ./BaseSpace_Download/ -print | grep fastq | grep $OUT | grep "R2" | xargs cat > ./Merged_FastQ_files/$OUT/$OUT"_R2.fastq.gz"

done