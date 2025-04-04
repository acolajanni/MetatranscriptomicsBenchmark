#!/bin/bash

#############################
#SBATCH -J readClassif
#SBATCH --array=0-3
#############################

# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################
module purge
module load r/4.2.3
module load use.own

MAIN_PATH=/shared/projects/microbiome_translocation/



task_array=(I L H K)
letter=${task_array[$SLURM_ARRAY_TASK_ID]}


echo $letter 

SCRIPT_DIR=${MAIN_PATH}/fastq_scripts/analyze_results
cd $SCRIPT_DIR


RESULTS_PATH=${MAIN_PATH}/results/${dataset}/

# Rscript --no-save --no-restore ${SCRIPT_DIR}/analyze_aligned_reads.R $letter

Rscript --no-save --no-restore ${SCRIPT_DIR}/analyze_aligned_reads.R $letter