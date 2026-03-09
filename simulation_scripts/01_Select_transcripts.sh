#!/bin/bash

#############################
#SBATCH -J simu
#SBATCH --array=0-9
#############################

module load r/4.2.3
module load seqkit



PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Simulation/with_replacement/}"
n_transcript="${2:-200}"
replicate="${3:-10}"
PROJECT_NAME="${4:-with_replacement}"
ONLY_HUMAN="${5:-FALSE}"

MAIN_PATH=/shared/projects/microbiome_translocation/
SCRIPT_DIR=${MAIN_PATH}fastq_scripts/Simulation/
cd $SCRIPT_DIR


if [[ $PROJECT_NAME = "family" ]] ; then
    family=TRUE
    PATH_DATABASE=${MAIN_PATH}database_clean/transcript_database/family/${SLURM_ARRAY_TASK_ID}/dedup_transcripts_ID/

else 
    family=FALSE
    PATH_DATABASE=${MAIN_PATH}database_clean/transcript_database/${SLURM_ARRAY_TASK_ID}/dedup_transcripts_ID/
fi





##### 1: Create R scripts to select the transcripts (Build simulated metagenome)
Rscript --no-save --no-restore ${SCRIPT_DIR}00_Build_simulated_metagenome_v3.r \
    $n_transcript $replicate $SLURM_ARRAY_TASK_ID \
    Simulation/${PROJECT_NAME}/ TRUE ${ONLY_HUMAN} $family


##### 2: take transcripts ID and extracting them into separate fasta file
replicate_dir=${PATH_DATA}transcripts/${SLURM_ARRAY_TASK_ID}/


# Build the transcripts
for phylum_file in ${replicate_dir}*_${n_transcript}.txt ; do

    phylum="${phylum_file%_*}" 
    phylum="${phylum##*/}"

    # Seqkit de-duplicate given ID
    n_id=$(sort $phylum_file | uniq -c | wc -l)

    # If there are exactly n_transcript unique transcripts, proceed without duplicates    
    if [ "$n_id" -eq "$n_transcript" ]; then
        seqkit grep --pattern-file $phylum_file ${PATH_DATABASE}${phylum}_transcripts.fasta > ${replicate_dir}/${phylum}_${n_transcript}.fasta
    else 
        # if duplicate ID, needs to be written 1 by 1
        rm -f ${replicate_dir}/${phylum}_${n_transcript}.fasta
        while IFS= read -r id; do
            # Step 1: Run seqkit grep as usual
            seqkit grep -p "$id" ${PATH_DATABASE}${phylum}_transcripts.fasta >> ${replicate_dir}/${phylum}_${n_transcript}.fasta

        done < $phylum_file

        # Step 2: Process the output file to add suffixes to duplicate transcript IDs
        awk '/^>/ {
                header=$0; 
                count[header]++; 
                if (count[header] > 1) {
                    print header "_" count[header] - 1;
                } else { 
                    print header;
                }
        } !/^>/ {print}' ${replicate_dir}/${phylum}_${n_transcript}.fasta > ${replicate_dir}/temp_${phylum}_${n_transcript}.fasta

        # Step 3: rename the file
        mv ${replicate_dir}/temp_${phylum}_${n_transcript}.fasta ${replicate_dir}/${phylum}_${n_transcript}.fasta
    
    fi

done








