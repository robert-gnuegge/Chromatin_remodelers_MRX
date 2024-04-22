#!/bin/bash

# info --------------------------------------------------------------------
# purpose: derive deduplicated alignments for NGS reads (FASTQ -> BAM)
# author: Robert Gnuegge (robert.gnuegge@gmail.com)
# date started: 04/22/24
# date last modified: 04/22/24

SAMPLE=$1
printf "\nProcessing $SAMPLE..."

READDIR="/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/S1-seq_and_MNase-seq/Rep2/S1-seq/01_Raw_data/Merged_FastQ_files"
# READDIR="/media/robert/Elements/Deep_sequencing_data/24-04-22-S1-seq/Merged_FastQ_files"
# READDIR="/media/robert/One Touch/Deep_sequencing_data/Robert_Gnuegge/24-04-22-S1-seq/Merged_FastQ_files"
IDX="S288C_R64-4-1_W303_SNPs_MATa_hocs2SrfIcs_hml_hmr_idx"
export BOWTIE2_INDEXES="/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/Reference_genome/03_Processed_data"
ADAPTERS="/home/robert/Research/Manuscripts/My_manuscripts/23-09-Chromatin_remodeler_Impact_on_MRX_nicking/Data_Analyses_and_Figures/Src/2nd_Adpt_linkers.fa"

# To keep time
SECONDS=0
Elapsed () {
	TIME="[$(printf "%02d" $(($SECONDS /3600))):$(printf "%02d" $(($SECONDS / 60))):$(printf "%02d" $(($SECONDS % 60)))]"
}

Elapsed
printf "\n$TIME Read2 length filtering..."
fastp --in1 "$READDIR"/"$SAMPLE"/"$SAMPLE"_R1.fastq.gz --in2 "$READDIR"/"$SAMPLE"/"$SAMPLE"_R2.fastq.gz \
--out1 tmp_R1.fastq.gz --out2 tmp_R2.fastq.gz --length_required=12 \
--disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g \
-h fastp_1.html 2> fastp_1.log # --reads_to_process 10000

Elapsed
printf "\n$TIME UMI extraction..."
fastp --in1 tmp_R1.fastq.gz --in2 tmp_R2.fastq.gz --out1 tmp2_R1.fastq.gz --out2 tmp2_R2.fastq.gz \
--umi --umi_loc=read2 --umi_len=12 \
--disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering \
-h fastp_2.html 2> fastp_2.log

rm tmp_R1.fastq.gz tmp_R2.fastq.gz tmp2_R2.fastq.gz

Elapsed
printf "\n$TIME Preprocessing and mapping Read1; converting and sorting alignments..."
fastp --in1 tmp2_R1.fastq.gz --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_fasta $ADAPTERS -h fastp_3.html 2>fastp_3.log --stdout \
| bowtie2 --local --very-sensitive -p 4 -x $IDX - 2> Bowtie2.log \
| samtools view -bu - \
| samtools sort - -o tmp.bam

rm tmp2_R1.fastq.gz

Elapsed
printf "\n$TIME Indexing..."
samtools index tmp.bam

samtools idxstats tmp.bam > "$SAMPLE"_alignment_summary.txt
Elapsed
printf "\n$TIME Removing alignments to chrM and 2 micron..."
samtools view -b tmp.bam chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI > tmp_2.bam

rm tmp.bam
rm tmp.bam.bai

Elapsed
printf "\n$TIME Indexing..."
samtools index tmp_2.bam

Elapsed
printf "\n$TIME Deduplicating and sorting..."
umi_tools dedup  --stdin=tmp_2.bam --umi-separator=":" --no-sort-output --log=umi_tools.log --method=unique \
| samtools sort - -o "$SAMPLE".bam

rm tmp_2.bam
rm tmp_2.bam.bai

Elapsed
printf "\n$TIME Indexing..."
samtools index "$SAMPLE".bam

Elapsed
printf "\n$TIME Done.\n\n"

# Notes -------------------------------------------------------------------
# * Read2 has to be length-filtered, as truncated UMIs break umi_tools dedup.
# * UMIs are moved from Read2 to the Read1 headers. Subsequently, only Read1 is processed further.
# * Alignments to mitochondrial and 2 micron DNAs are not needed for further analyses.
#   We record the number of alignments to each cellular DNA (using samtools idxstats).
#   We then exclude alignments to mitochondrial and 2 micron DNAs from further processing to reduce run time.