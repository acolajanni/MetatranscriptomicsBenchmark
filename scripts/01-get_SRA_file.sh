#!/bin/bash

#############################
#SBATCH -J get_sra
#SBATCH --array=0-54
#############################

################################################################################
# > september 2023
# > Script : get SRA file.sh
# > Function : download sra file then the correpsonding fastq
# @ COLAJANNI Antonin
################################################################################

module load sra-tools/2.11.0


# First argument: path to data, if no argument go to douek_cell_2021
OUTPUT_DIR="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"


declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}
echo $SRA_ID


# Fetch sra files
prefetch $SRA_ID --output-directory $OUTPUT_DIR


mkdir -p ${OUTPUT_DIR}raw_reads/


fastq-dump ${OUTPUT_DIR}${SRA_ID}/${SRA_ID}.sra --gzip --clip --skip-technical --split-files --split-spot --outdir $OUTPUT_DIR"raw_reads/" 

    
if [ $SEQ_TYPE = "PE" ] ; then 
    gzip ${OUTPUT_DIR}raw_reads/${SRA_ID}_2.fastq &
    gzip ${OUTPUT_DIR}raw_reads/${SRA_ID}_1.fastq
else
    gzip ${OUTPUT_DIR}raw_reads/${SRA_ID}.fastq
fi

echo " --- "



remove downloaded files
rm -r  ${OUTPUT_DIR}${SRA_ID}/