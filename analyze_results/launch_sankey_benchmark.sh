#!/bin/bash

#############################
#SBATCH -J sankey
#SBATCH --array=0-1
############################

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


PROJECT_DIR=/shared/projects/microbiome_translocation/
SCRIPT_DIR=${PROJECT_DIR}fastq_scripts/
cd $PROJECT_DIR


RESULTS_PATH=${PROJECT_DIR}results/Douek_Cleveland/
ID=${PROJECT_DIR}data/Douek_Cleveland/sra_list_RNA.txt


declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

SRA_IDs=(responder non_responder)

SRA=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}
#SRA=all



Rscript --no-save --no-restore \
        /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/99-Sankey_benchmark_real.R \
        $SRA Douek_Cleveland

        