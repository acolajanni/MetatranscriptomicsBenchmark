#!/bin/bash
#############################
#SBATCH -J ContToReads
#############################

module purge
module load seqkit
module load samtools

# MAIN_PATH="${1:-/shared/projects/microbiome_translocation/data/Simulation/A0}"
# PATH_RESULTS="${2-/shared/projects/microbiome_translocation/results/Simulation/Contigs/A0/}"

#MAIN_PATH="${2-/home/acolajanni/Documents/work/data/Simulation/A12/}"
PATH_RES="${1:-/home/acolajanni/Documents/work/results/Simulation/Contigs/}"
SRA_ID="${2:-A12}"
TAX_lvl="${3:-strain}"

PATH_RESULTS=${PATH_RES}${SRA_ID}/


echo $SRA_ID

mkdir -p ${PATH_RESULTS}/ContigsToReads/
PATH_classif=${PATH_RESULTS}/ContigsToReads/

PATH_quant=${PATH_RESULTS}Quantification/${TAX_lvl}/

cd $PATH_classif

# Index/sort SAM file produced by Salmon
samtools view -@ 20 -u ${PATH_quant}/${TAX_lvl}_quant_${SRA_ID}.sam > ${PATH_quant}/${TAX_lvl}_quant_${SRA_ID}.bam
samtools sort -@ 20 ${PATH_quant}/${TAX_lvl}_quant_${SRA_ID}.bam > ${PATH_quant}/${TAX_lvl}_quant_${SRA_ID}.sorted
samtools index -@ 20 ${PATH_quant}/${TAX_lvl}_quant_${SRA_ID}.sorted    

BAM_file=${PATH_quant}/${TAX_lvl}_quant_${SRA_ID}.sorted 



# seqkit grep --pattern-file ${pattern_file} \
#     ${PATH_RESULTS}trinity_out_dir.Trinity.fasta > \
#     ${PATH_classif}/${TAX_lvl}_contigs.fasta

#samtools view -@ 20 $BAM_file ${PATH_classif}/${TAX_lvl}_contigs.fasta | awk '{print $1 "\t" $3}' > ${PATH_classif}${TAX_lvl}_readsToContigs.txt


# ==> Refiltrer parmis les outputs en fonction des contigs ID
samtools view -@ 20 $BAM_file | awk '{print $1 "\t" $3}' > ${PATH_classif}${TAX_lvl}_readsToContigs.txt