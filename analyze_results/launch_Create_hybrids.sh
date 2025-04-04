#!/bin/bash

#############################
#SBATCH -J make_hybrid
#############################

PROJECT_DIR=/shared/projects/microbiome_translocation/
SCRIPT_DIR=${PROJECT_DIR}fastq_scripts/
cd $PROJECT_DIR


RESULTS_PATH=${PROJECT_DIR}results/Douek_Cleveland/


ID=${PROJECT_DIR}data/Douek_Cleveland/sra_list_RNA.txt
ID=${PROJECT_DIR}data/Douek_Cleveland/sra_list_RNA_error.txt


Classif_contigs=${RESULTS_PATH}Contigs/
Classif_kuniq=${RESULTS_PATH}kraken/microbialDB/kuniq/
Classif_kraken2=${RESULTS_PATH}kraken/nt/k2uniq/


sbatch -A microbiome_translocation --mem=8G --cpus-per-task=2 -J blastKU -t 01:00:00 ${PROJECT_DIR}fastq_scripts/12-Hybrid_classif.sh \
    ${RESULTS_PATH} $Classif_contigs $Classif_kuniq $ID hybrid_Blast-kuniq

sbatch -A microbiome_translocation --mem=8G --cpus-per-task=2 -J blastK2 -t 01:00:00 ${PROJECT_DIR}fastq_scripts/12-Hybrid_classif.sh \
    ${RESULTS_PATH} $Classif_contigs $Classif_kraken2 $ID hybrid_Blast-kraken2

sbatch -A microbiome_translocation --mem=8G --cpus-per-task=2 -J ku-k2 -t 01:00:00 ${PROJECT_DIR}fastq_scripts/12-Hybrid_classif.sh \
    ${RESULTS_PATH} $Classif_kuniq $Classif_kraken2 $ID hybrid_kuniq-kraken2

sbatch -A microbiome_translocation --mem=8G --cpus-per-task=2 -J k2-ku -t 01:00:00 ${PROJECT_DIR}fastq_scripts/12-Hybrid_classif.sh \
    ${RESULTS_PATH} $Classif_kraken2 $Classif_kuniq $ID hybrid_kraken2-kuniq