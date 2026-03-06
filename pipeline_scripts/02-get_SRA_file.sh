#!/bin/bash

#############################
#SBATCH -J get_sra
#SBATCH --array=0-9
#############################

################################################################################
# > september 2023
# > Script : get SRA file.sh
# > Function : download sra file then the correpsonding fastq
# @ COLAJANNI Antonin
################################################################################

module load sra-tools/2.11.0


# First argument: path to data, if no argument go to douek_cell_2021
OUTPUT_DIR="${1:-~/data/Douek_cell2021/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-~/data/Douek_cell2021/sra_subset_in_manuscript.txt}"


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

# fasterq-dump ${OUTPUT_DIR}${SRA_ID}/${SRA_ID}.sra --skip-technical --outdir $OUTPUT_DIR"raw_reads/" --threads 8   # SRA file are named like SRA_id.sra / at the end we got fastq file
    
# if [ $SEQ_TYPE = "PE" ] ; then 
#     gzip ${OUTPUT_DIR}raw_reads/${SRA_ID}_2.fastq &
#     gzip ${OUTPUT_DIR}raw_reads/${SRA_ID}_1.fastq
# else
#     gzip ${OUTPUT_DIR}raw_reads/${SRA_ID}.fastq
# fi

echo " --- "



# remove downloaded files
rm -r  ${OUTPUT_DIR}${SRA_ID}/