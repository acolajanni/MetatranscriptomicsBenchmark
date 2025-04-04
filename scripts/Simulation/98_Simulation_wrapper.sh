#!/bin/bash


#############################
#SBATCH -J simu
#############################

module load r/4.2.3


PATH_MAIN=/shared/projects/microbiome_translocation/
SCRIPT_DIR=${PATH_MAIN}fastq_scripts/Simulation/
PROJECT_NAME=with_replacement


PATH_DATA=${PATH_MAIN}data/Simulation/${PROJECT_NAME}/




### 0: Reunite all the transcripts from each phylum into one file and gather information on them
sbatch -A microbiome_translocation --mem=16G --cpus-per-task=12 ${SCRIPT_DIR}00_Build_transcript_database.sh



Rscript --no-save --no-restore ${SCRIPT_DIR}00_Get_transcripts_info.r $SLURM_ARRAY_TASK_ID


### 1: Randomly select and isolate n_transcript
# 10-90 / 50-50 and 90-10 /// 45*37 phylum = 1665 transcripts - 185/1665=10/90 - 14985/1665=90/10
sbatch -A microbiome_translocation --mem=16G --cpus-per-task=4 --exclude=cpu-node-17 ${SCRIPT_DIR}01_Select_transcripts.sh \
     ${PATH_DATA} 45 10 $PROJECT_NAME FALSE

sbatch -A microbiome_translocation --mem=16G --cpus-per-task=4 --exclude=cpu-node-17 ${SCRIPT_DIR}01_Select_transcripts.sh \
     ${PATH_DATA} 185 10 $PROJECT_NAME TRUE

sbatch -A microbiome_translocation --mem=16G --cpus-per-task=4 --exclude=cpu-node-17 ${SCRIPT_DIR}01_Select_transcripts.sh \
     ${PATH_DATA} 1665 10 $PROJECT_NAME TRUE

sbatch -A microbiome_translocation --mem=16G --cpus-per-task=4 --exclude=cpu-node-17 ${SCRIPT_DIR}01_Select_transcripts.sh \
     ${PATH_DATA} 14985 10 $PROJECT_NAME TRUE



## 2: Simulate reads from these transcripts (5 and 10 reads/transcripts) (1 job per replicate set: 10)
sbatch -A microbiome_translocation --mem=16G --cpus-per-task=10 --exclude=cpu-node-17 ${SCRIPT_DIR}02_launch_reads_simulation.sh \
    $PATH_DATA 45 100 $PROJECT_NAME
 
sbatch -A microbiome_translocation --mem=16G --cpus-per-task=10 --exclude=cpu-node-17 ${SCRIPT_DIR}02_launch_reads_simulation.sh \
    $PATH_DATA 185 100 $PROJECT_NAME

sbatch -A microbiome_translocation --mem=32G --cpus-per-task=10 --exclude=cpu-node-17 ${SCRIPT_DIR}02_launch_reads_simulation.sh \
    $PATH_DATA 1665 100 $PROJECT_NAME

sbatch -A microbiome_translocation --mem=32G --cpus-per-task=10 --exclude=cpu-node-17 ${SCRIPT_DIR}02_launch_reads_simulation.sh \
    $PATH_DATA 14985 100 $PROJECT_NAME


sbatch -A microbiome_translocation --mem=16G --cpus-per-task=4 --exclude=cpu-node-17 ${SCRIPT_DIR}03_Build_Sample_OR_human.sh \
    ${PATH_DATA}condition_list.tsv \
    ${PATH_DATA} unmapped
