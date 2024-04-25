#!/bin/bash

TARGET_DIR="/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/S1-seq_and_MNase-seq/Rep2/MNase-seq/01_Raw_data"
OUT_DIR=$TARGET_DIR"/BaseSpace_Download"
mkdir -p $OUT_DIR

# bs auth	
bs list projects
read -p "Specify the ID of the project to be downloaded: " ID
bs download project -i $ID -o $OUT_DIR --extension=fastq.gz