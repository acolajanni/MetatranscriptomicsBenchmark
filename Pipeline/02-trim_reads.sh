#! /bin/bash

################################################################################
# > september 2023
# > Script : trim_reads.sh
# > Function : Trim fastq file <ith trimmomatic + fastqc
# @ COLAJANNI Antonin
################################################################################

#############################
#SBATCH -J SimuTRI
#SBATCH --array=0-54%10
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

module load trimmomatic/0.39


PATH_RES="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
SEQ_TYPE="${2:-PE}"
RAW_extension="${3:-}"
ID="${4:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"
adapter_path="${5:-/shared/projects/microbiome_translocation/data/TruSeq3-PE.fa}" 
INPUT_FOLDERNAME="${6:-raw_reads}"

#PATH_RES=/shared/projects/microbiome_translocation/data/Douek_Cleveland/TEST_RAM/
#FOLDERS=${PATH_RES}/${ID}*/

DATA_DIR=/shared/projects/microbiome_translocation/data/
mkdir -p $(echo $PATH_RES"Trimmed_reads/")
OUTPUT_PATH=$(echo $PATH_RES"Trimmed_reads/")

# Storing folder names in an array 
declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

# Number of job in the job Array
length=${#SRA_IDs[@]}
# length divided by n, where n is the number of task asked
n=54
#n_operation=$(( (length / n) + 1 ))
n_operation=$(( (length / n) )) 
index_iter=$((n_operation * SLURM_ARRAY_TASK_ID))

SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}



echo $SRA_ID
    
    
if [ ${SEQ_TYPE} = "PE" ] ; then 

        echo - Paired-end -

        READ1_path=${PATH_RES}${INPUT_FOLDERNAME}/${SRA_ID}_1.fastq${RAW_extension}
        READ2_path=${PATH_RES}${INPUT_FOLDERNAME}/${SRA_ID}_2.fastq${RAW_extension}

        output_path_R1=${OUTPUT_PATH}${SRA_ID}_1
        output_path_R2=${OUTPUT_PATH}${SRA_ID}_2

        mkdir -p $(echo $PATH_RES"Trimmed_reads/Orphaned_reads/")
        #adapter_path=${DATA_DIR}TruSeq3-PE.fa

    
        trimmomatic PE -threads 20 -phred33  $READ1_path $READ2_path \
               ${output_path_R1}.trimmed.fastq.gz ${OUTPUT_PATH}Orphaned_reads/${SRA_ID}_1_orphanedReads.fastq.gz \
               ${output_path_R2}.trimmed.fastq.gz ${OUTPUT_PATH}Orphaned_reads/${SRA_ID}_2_orphanedReads.fastq.gz \
               ILLUMINACLIP:$adapter_path:2:40:15 \
               LEADING:2 TRAILING:2 \
               SLIDINGWINDOW:4:15 \
               MINLEN:36
                
                

        
elif [ ${SEQ_TYPE} = "SE" ] ; then 

        echo - Single-end -

        READ_input=${PATH_RES}/raw_reads/${SRA_ID}".fastq"${RAW_extension}
        READ_output=${OUTPUT_PATH}${SRA_ID}".trimmed.fastq"${RAW_extension}

        adapter_path=${DATA_DIR}TruSeq2-SE.fa

        trimmomatic SE -threads 20 -phred33 \
                $READ_input $READ_output \
                ILLUMINACLIP:$adapter_path:2:40:15 \
                LEADING:2 TRAILING:2 \
                SLIDINGWINDOW:4:15 \
                MINLEN:36

        #gzip $READ_output


fi

echo " --- "
