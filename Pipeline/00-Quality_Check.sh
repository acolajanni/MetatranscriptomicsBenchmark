#! /bin/sh

################################################################################
# > september 2023
# > Script : 00_QC.sh
# > Function : Quality check
# @ COLAJANNI Antonin
################################################################################


PATH_RES="${1:-/shared/projects/microbiome_translocation/data/Douek_cell2021/}"
RAW_extension="${2:-}"


### Run FastQC ###

# On raw reads
mkdir -p ${PATH_RES}"raw_reads/fastqc/"
files=${PATH_RES}"raw_reads/*.fastq"${RAW_extension}
# fastqc -o ${PATH_RES}"raw_reads/fastqc/" -t 20 $files
# multiqc ${PATH_RES}"raw_reads/fastqc/" -o ${PATH_RES}"raw_reads/" -n raw_report

# On trimmed reads
files=${PATH_RES}Trimmed_reads/*.trimmed.fastq${RAW_extension}  
mkdir -p ${PATH_RES}"Trimmed_reads/fastqc/"
fastqc -o ${PATH_RES}"Trimmed_reads/fastqc"/ -t 20 $files
multiqc ${PATH_RES}"Trimmed_reads/fastqc/" -o ${PATH_RES}"Trimmed_reads/" -n trimmed_report
##################



mkdir -p ${PATH_RES}fastqc_bam/
mkdir -p ${PATH_RES}fastqc_bam/multiqc/
files=${PATH_RES}Bowtie2_mapping_hg38/${ID}*/*.bam


# fastqc -o ${PATH_RES}fastqc_bam/ -t 20 $files
# multiqc ${PATH_RES}fastqc_bam/ -o ${PATH_RES}fastqc_bam/multiqc/ -n mapping_bwt2hg38_report


# PATH_results=/shared/projects/microbiome_translocation/results/VRC_IVI/wholeblood/Benchmark/krakenuniq_4filters/
# files=${PATH_results}/*_RNA/*classified#.fastq

# mkdir -p $PATH_results/fastqc/


#fastqc -o ${PATH_results}fastqc/ -t 20 --nogroup $files 
#multiqc ${PATH_results}fastqc/ -o ${PATH_RES}fastqc/ -n mapping_Kraken