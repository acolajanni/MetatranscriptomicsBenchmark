#!/bin/bash

#############################
#SBATCH -J launch_simu
#SBATCH --array=0-9
#############################

module purge
module load r/4.2.3


PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Simulation/with_replacement/}"
n_transcript="${2:-200}"
read_per_transcript="${3:-5}"


PATH_MAIN=/shared/projects/microbiome_translocation/
cd $PATH_MAIN

echo ${PATH_DATA}transcripts/${SLURM_ARRAY_TASK_ID}/

Rscript --no-save --no-restore \
    ${PATH_MAIN}fastq_scripts/Simulation/02_Simulate_reads_v2.r \
        ${PATH_DATA} \
        ${PATH_DATA}transcripts/${SLURM_ARRAY_TASK_ID}/ \
        ${n_transcript} \
        ${read_per_transcript} \
        ${PATH_MAIN}outputFile_${sra_id}.Rout 2>&1


# For each phylum, merge reads into one file
read_path=${PATH_DATA}raw_reads_${n_transcript}/${SLURM_ARRAY_TASK_ID}/
for file in ${read_path}*/; do 

    phylum=$(basename "$file")
    echo $phylum

    phylum_path=${read_path}/${phylum}/

    cat ${phylum_path}/tmp_merge/${phylum}_${read_per_transcript}rt_*_R1.fastq.gz > ${phylum_path}/${phylum}_${read_per_transcript}rt_R1.fastq.gz 
    cat ${phylum_path}/tmp_merge/${phylum}_${read_per_transcript}rt_*_R2.fastq.gz > ${phylum_path}/${phylum}_${read_per_transcript}rt_R2.fastq.gz 

    #rm ${phylum_path}/tmp_merge/${phylum}_${read_per_transcript}rt_*_R*.fastq.gz
done



