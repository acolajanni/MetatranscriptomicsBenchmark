#!/bin/bash

#############################
#SBATCH -J simu
#SBATCH --array=0-49
#############################

module purge
module load seqkit
module load r/4.2.3


PATH_DATA="${1:-~/data/Simulation/realist/}"

read_per_transcript="100"


PATH_MAIN=~/


replicate_ID=$(( SLURM_ARRAY_TASK_ID))
echo $replicate_ID

compression=.gz
extension=.fastq

SIMULATION_DIR=${PATH_DATA}raw_reads/
mkdir -p ${PATH_DATA}raw_reads/${replicate_ID}
OUTPUT_DIR=${PATH_DATA}raw_reads/${replicate_ID}/

R1_out=${OUTPUT_DIR}${replicate_ID}_unmapped_1${extension}
R2_out=${OUTPUT_DIR}${replicate_ID}_unmapped_2${extension}

rm $R1_out
touch $R1_out

rm $R2_out
touch $R2_out

# For each phylum, gather the previously producted reads except for human if we don't want them
for phylum_path in ${SIMULATION_DIR}${replicate_ID}/* ; do
    phylum="${phylum_path##*/}"
    echo $phylum


    R1=${phylum_path}/${phylum}_${read_per_transcript}rt_R1${extension}${compression}
    R2=${phylum_path}/${phylum}_${read_per_transcript}rt_R2${extension}${compression}

    echo $R1

    # carefull Seqkit removes the compression

    # ADD /1 and /2 at the end of read header to prevent potential problems downstream
    seqkit replace --threads 4 -p '(.+)' -r '${1}/1' $R1 >> $R1_out
    seqkit replace --threads 4 -p '(.+)' -r '${1}/2' $R2 >> $R2_out     

done






