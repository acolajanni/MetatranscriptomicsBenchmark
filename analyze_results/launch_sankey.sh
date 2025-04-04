#!/bin/bash

#############################
#SBATCH -J metrics
#SBATCH --array=0-3
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
ID="${1:-/shared/projects/microbiome_translocation/data/Simulation/family/Simu_condition.txt}"
dataset="${2:-Simulation}"
result_dir="${3:/shared/projects/microbiome_translocation/results/Simulation/family/}"

declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

compare_assembled_reads_only=FALSE

#task_array=(A C E M N O P Q R)
task_array=(M N O P R)
task_array=(A B C D)
#task_array=(M N O P Q R)
letter=${task_array[$SLURM_ARRAY_TASK_ID]}
SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}

echo $letter 

result_dir=/shared/projects/microbiome_translocation/results/Simulation/family/

cd /shared/projects/microbiome_translocation/
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/99-Get_sankey.R ${SRA_ID} $dataset \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1


# cd /shared/projects/microbiome_translocation/
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/99-Get_sankey.R ${SRA_ID} Douek_Cleveland \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1


# Rscript --no-save --no-restore /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/99-Get_sankey.R all Douek_Cleveland \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1


# echo "Kraken2"

# # Kraken2:
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} normal kraken2 none aligned $compare_assembled_reads_only \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"
# echo "Krakenuniq"

# ### KrakenUniq:
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} normal kuniq none aligned $compare_assembled_reads_only \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"

# # # ### Contig based:
# echo "Contigs"

# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} normal Blast none aligned $compare_assembled_reads_only \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"




## hybrid: 
echo "hybrid ku"
Rscript --no-save --no-restore \
        /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
        ${letter} ${result_dir} hybrid Blast kuniq aligned \
        /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
echo "---------------------------------------------"
echo "hybrid k2"

Rscript --no-save --no-restore \
        /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
        ${letter} ${result_dir} hybrid Blast kraken2 aligned \
        /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
echo "---------------------------------------------"

# echo "hybrid ku-k2"
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} hybrid kuniq kraken2 aligned $compare_assembled_reads_only \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"

# echo "hybrid k2-ku"
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} hybrid kraken2 kuniq aligned $compare_assembled_reads_only \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"



# ### human kuniq first: 
# echo "hybrid ku - human kuniq"
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} hybrid Blast kuniq human_kuniq \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"
# echo "hybrid k2- human kuniq"

# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} hybrid Blast kraken2 human_kuniq \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"

# echo "hybrid ku-k2- human kuniq"
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} hybrid kuniq kraken2 human_kuniq \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"

# echo "hybrid k2-ku - human kuniq"
# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} hybrid kraken2 kuniq human_kuniq \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"




# # Kraken2:
# echo "k2 - human kuniq"

# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} normal kraken2 none human_kuniq \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"

# # # ### Contig based:
# echo "Contigs - human kuniq"

# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/Simulation/12_Analyze_classif_SalmonKraken.r \
#         ${letter} ${result_dir} normal Blast none human_kuniq \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
# echo "---------------------------------------------"









# Rscript --no-save --no-restore \
#         /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/Standardize_classification.r \
#         ${SRA_ID} $dataset \
#         /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1