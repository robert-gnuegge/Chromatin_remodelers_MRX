#!/bin/bash

declare -a arr=( #"LSY4518-13B_0h" "LSY4518-13B_1h" "LSY4518-13B_2h" "LSY4518-13B_4h")
	# "LSY5415_0h" "LSY5415_1h" "LSY5415_2h" "LSY5415_4h")
	# "LSY5934_0h" "LSY5934_1h" "LSY5934_2h" "LSY5934_4h")
	# "LSY5935_0h" "LSY5935_1h" "LSY5935_2h" 
	"LSY5935_4h")

BASEDIR="/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/S1-seq_and_MNase-seq/Rep2/MNase-seq"

for SAMPLE in "${arr[@]}"
do

	OUTDIR="$BASEDIR"/03_Processed_data/BAM/"$SAMPLE"
	mkdir -p $OUTDIR
	cd $OUTDIR

	bash "$BASEDIR"/02_Scripts/03_From_fastq_to_bam.sh $SAMPLE

done