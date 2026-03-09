#!/bin/bash

#############################
#SBATCH -J simu
#SBATCH --array=0-49
#############################

module load r/4.2.3
module load seqkit



PATH_DATA="${1:-~//data/Simulation/realist/}"
PROJECT_NAME="${2:-realist}"
ONLY_HUMAN="${3:-FALSE}"


MAIN_PATH=~//
SCRIPT_DIR=${MAIN_PATH}fastq_scripts/Simulation/
cd $SCRIPT_DIR

PATH_DATABASE=${MAIN_PATH}database_clean/transcript_database/realist/${SLURM_ARRAY_TASK_ID}/dedup_transcripts_ID/

##### 1: Create R scripts to select the transcripts (Build simulated metagenome)
Rscript --no-save --no-restore ${SCRIPT_DIR}realist/00_Build_simulated_metagenome_realist.r 


##### 2: take transcripts ID and extracting them into separate fasta file
replicate_dir=${PATH_DATA}transcripts/${SLURM_ARRAY_TASK_ID}/


# Build the transcripts
for genus_file in ${replicate_dir}*.txt ; do

    # Error in this line ??
    genus="${genus_file%_*}" 
    genus="${genus##*/}"

    # Seqkit de-duplicate given ID
    n_id=$(sort $genus_file | uniq -c | wc -l)
    n_transcript=$(sort $genus_file | wc -l)



    if [ "$genus" == "Homo" ]; then
        trancripts_file=${MAIN_PATH}database_clean/transcript_database/realist/Homo_sapiens_rna_transcripts.fasta

    else 
        trancripts_file=${PATH_DATABASE}${genus}_transcripts.fasta
    fi

    # If there are exactly n_transcript unique transcripts, proceed without duplicates    
    if [ "$n_id" -eq "$n_transcript" ]; then
        seqkit grep --pattern-file $genus_file $trancripts_file > ${replicate_dir}/${genus}_${n_transcript}.fasta
    else 
        # if duplicate ID, needs to be written 1 by 1
        rm -f ${replicate_dir}/${genus}_${n_transcript}.fasta
        while IFS= read -r id; do
            # Step 1: Run seqkit grep as usual
            seqkit grep -p "$id" $trancripts_file >> ${replicate_dir}/${genus}_${n_transcript}.fasta

        done < $genus_file

        # Step 2: Process the output file to add suffixes to duplicate transcript IDs
        awk '/^>/ {
                header=$0; 
                count[header]++; 
                if (count[header] > 1) {
                    print header "_" count[header] - 1;
                } else { 
                    print header;
                }
        } !/^>/ {print}' ${replicate_dir}/${genus}_${n_transcript}.fasta > ${replicate_dir}/temp_${genus}_${n_transcript}.fasta

        # Step 3: rename the file
        mv ${replicate_dir}/temp_${genus}_${n_transcript}.fasta ${replicate_dir}/${genus}_${n_transcript}.fasta
    
    fi




done








