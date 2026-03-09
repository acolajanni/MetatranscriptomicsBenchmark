#!/bin/bash

#############################
#SBATCH -J launch_simu
#SBATCH --array=0-49
#############################

module purge
module load r/4.2.3



PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Simulation/realist/}"
n_transcript="${2:-200}"
read_per_transcript="${3:-100}"
transcript_dir="${4:-/shared/projects/microbiome_translocation/data/Simulation/realist/transcripts/}"
array="${5:-TRUE}"


if [[ "$array" == "TRUE" ]] ; then
    transcript_dir=${transcript_dir}/${SLURM_ARRAY_TASK_ID}/
fi


PATH_MAIN=/shared/projects/microbiome_translocation/
cd $PATH_MAIN

echo ${PATH_DATA}transcripts/${SLURM_ARRAY_TASK_ID}/

Rscript --no-save --no-restore \
    ${PATH_MAIN}fastq_scripts/Simulation/realist/01_Simulate_reads_realist.r \
        $SLURM_ARRAY_TASK_ID \
        ${PATH_MAIN}outputFile_${sra_id}.Rout 2>&1


# For each phylum, merge reads into one file
read_path=${PATH_DATA}raw_reads/${SLURM_ARRAY_TASK_ID}/
for file in ${read_path}*/; do 

    phylum=$(basename "$file")
    echo $phylum

    phylum_path=${read_path}/${phylum}/

    cat ${phylum_path}/tmp_merge/${phylum}_${read_per_transcript}rt_*_R1.fastq.gz > ${phylum_path}/${phylum}_${read_per_transcript}rt_R1.fastq.gz 
    cat ${phylum_path}/tmp_merge/${phylum}_${read_per_transcript}rt_*_R2.fastq.gz > ${phylum_path}/${phylum}_${read_per_transcript}rt_R2.fastq.gz 

    #rm ${phylum_path}/tmp_merge/${phylum}_${read_per_transcript}rt_*_R*.fastq.gz
done



